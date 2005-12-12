//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/OhmmsBlas.h"

namespace qmcplusplus {

    DiracDeterminantBase::DiracDeterminantBase(int first): 
      NP(0), FirstIndex(first) {}

    ///default destructor
    DiracDeterminantBase::~DiracDeterminantBase() {}
  
    ///reset the size: with the number of particles and number of orbtials
    void DiracDeterminantBase::resize(int nel, int morb) {
      int norb=morb;
      if(norb <= 0) norb = nel; // for morb == -1 (default)
      psiM.resize(nel,norb);
      dpsiM.resize(nel,norb);
      d2psiM.resize(nel,norb);
      psiM_temp.resize(nel,norb);
      dpsiM_temp.resize(nel,norb);
      d2psiM_temp.resize(nel,norb);
      psiMinv.resize(nel,norb);
      psiV.resize(norb);
      WorkSpace.resize(nel);
      Pivot.resize(nel);
      LastIndex = FirstIndex + nel;
      NumPtcls=nel;
      NumOrbitals=norb;
    }

    DiracDeterminantBase::ValueType 
    DiracDeterminantBase::registerData(ParticleSet& P, PooledData<RealType>& buf) {

      if(NP == 0) {//first time, allocate once
	//int norb = cols();
	dpsiV.resize(NumOrbitals);
	d2psiV.resize(NumOrbitals);
	workV1.resize(NumOrbitals);
	workV2.resize(NumOrbitals);
	NP=P.getTotalNum();
	myG.resize(NP);
	myL.resize(NP);
	myG_temp.resize(NP);
	myL_temp.resize(NP);
	FirstAddressOfG = &myG[0][0];
	LastAddressOfG = FirstAddressOfG + NP*DIM;
	FirstAddressOfdV = &(dpsiM(0,0)[0]); //(*dpsiM.begin())[0]);
	LastAddressOfdV = FirstAddressOfdV + NumPtcls*NumOrbitals*DIM;
      }

      //allocate once but each walker calls this
      myG=0.0;
      myL=0.0;

      ValueType x=evaluate(P,myG,myL); 

      P.G += myG;
      P.L += myL;

      //add the data: determinant, inverse, gradient and laplacians
      buf.add(psiM.begin(),psiM.end());
      buf.add(FirstAddressOfdV,LastAddressOfdV);
      buf.add(d2psiM.begin(),d2psiM.end());
      buf.add(myL.begin(), myL.end());
      buf.add(FirstAddressOfG,LastAddressOfG);
      buf.add(CurrentDet);

      return CurrentDet;
    }

    DiracDeterminantBase::ValueType 
      DiracDeterminantBase::updateBuffer(ParticleSet& P, PooledData<RealType>& buf) {

      myG=0.0;
      myL=0.0;
      ValueType x=evaluate(P,myG,myL); 
      P.G += myG;
      P.L += myL;
      buf.put(psiM.begin(),psiM.end());
      buf.put(FirstAddressOfdV,LastAddressOfdV);
      buf.put(d2psiM.begin(),d2psiM.end());
      buf.put(myL.begin(), myL.end());
      buf.put(FirstAddressOfG,LastAddressOfG);
      buf.put(CurrentDet);
      return CurrentDet;
    }

    void DiracDeterminantBase::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {

      buf.get(psiM.begin(),psiM.end());
      buf.get(FirstAddressOfdV,LastAddressOfdV);
      buf.get(d2psiM.begin(),d2psiM.end());
      buf.get(myL.begin(), myL.end());
      buf.get(FirstAddressOfG,LastAddressOfG);
      buf.get(CurrentDet);

      //re-evaluate it for testing
      //Phi.evaluate(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM);
      //CurrentDet = Invert(psiM.data(),NumPtcls,NumOrbitals);
      //need extra copy for gradient/laplacian calculations without updating it
      psiM_temp = psiM;
      dpsiM_temp = dpsiM;
      d2psiM_temp = d2psiM;
    }

    DiracDeterminantBase::ValueType 
    DiracDeterminantBase::ratio(ParticleSet& P, int iat) {
      //evaluateSingle(P, iat, psiV);
      evaluateSingle(P, iat);
#ifdef DIRAC_USE_BLAS
      return curRatio = BLAS::dot(NumOrbitals,psiM[iat-FirstIndex],&psiV[0]);
#else
      return DetRatio(psiM, psiV.begin(),iat-FirstIndex);
#endif
    }

    DiracDeterminantBase::ValueType 
    DiracDeterminantBase::ratio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG, 
		    ParticleSet::ParticleLaplacian_t& dL) {

      //evaluateSPO(P, iat, psiV, dpsiV, d2psiV);
      evaluateSingleAll(P, iat);

      WorkingIndex = iat-FirstIndex;
#ifdef DIRAC_USE_BLAS
      curRatio = BLAS::dot(NumOrbitals,psiM_temp[WorkingIndex],&psiV[0]);
#else
      curRatio= DetRatio(psiM_temp, psiV.begin(),WorkingIndex);
#endif

      //update psiM_temp with the row substituted
      DetUpdate(psiM_temp,psiV,workV1,workV2,WorkingIndex,curRatio);

      //update dpsiM_temp and d2psiM_temp 
      for(int j=0; j<NumOrbitals; j++) {
	dpsiM_temp(WorkingIndex,j)=dpsiV[j];
	d2psiM_temp(WorkingIndex,j)=d2psiV[j];
      }

      int kat=FirstIndex;
      for(int i=0; i<NumPtcls; i++,kat++) {
	PosType rv =psiM_temp(i,0)*dpsiM_temp(i,0);
	ValueType lap=psiM_temp(i,0)*d2psiM_temp(i,0);
        for(int j=1; j<NumOrbitals; j++) {
	  rv += psiM_temp(i,j)*dpsiM_temp(i,j);
	  lap += psiM_temp(i,j)*d2psiM_temp(i,j);
	}
	lap -= dot(rv,rv);
	dG[kat] += rv - myG[kat];  myG_temp[kat]=rv;
	dL[kat] += lap -myL[kat];  myL_temp[kat]=lap;
      }

      return curRatio;
    }

    DiracDeterminantBase::ValueType 
    DiracDeterminantBase::logRatio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG, 
		    ParticleSet::ParticleLaplacian_t& dL) {
      ValueType r=ratio(P,iat,dG,dL);
      SignValue = (r<0.0)? -1.0: 1.0;
      return log(abs(r));
    }


    /** move was accepted, update the real container
     */
    void DiracDeterminantBase::update(ParticleSet& P, int iat) {
      CurrentDet *= curRatio;
      myG = myG_temp;
      myL = myL_temp;
      psiM = psiM_temp;
      for(int j=0; j<NumOrbitals; j++) {
	dpsiM(WorkingIndex,j)=dpsiV[j];
	d2psiM(WorkingIndex,j)=d2psiV[j];
      }
      curRatio=1.0;
    }

    /** move was rejected. copy the real container to the temporary to move on
     */
    void DiracDeterminantBase::restore(int iat) {
      psiM_temp = psiM;
      for(int j=0; j<NumOrbitals; j++) {
	dpsiM_temp(WorkingIndex,j)=dpsiM(WorkingIndex,j);
	d2psiM_temp(WorkingIndex,j)=d2psiM(WorkingIndex,j);
      }
      curRatio=1.0;
    }

    
    void DiracDeterminantBase::update(ParticleSet& P, 
		ParticleSet::ParticleGradient_t& dG, 
		ParticleSet::ParticleLaplacian_t& dL,
		int iat) {

      DetUpdate(psiM,psiV,workV1,workV2,WorkingIndex,curRatio);
      for(int j=0; j<NumOrbitals; j++) {
	dpsiM(WorkingIndex,j)=dpsiV[j];
	d2psiM(WorkingIndex,j)=d2psiV[j];
      }

      int kat=FirstIndex;
      for(int i=0; i<NumPtcls; i++,kat++) {
	PosType rv =psiM(i,0)*dpsiM(i,0);
	ValueType lap=psiM(i,0)*d2psiM(i,0);
	for(int j=1; j<NumOrbitals; j++) {
	  rv += psiM(i,j)*dpsiM(i,j);
	  lap += psiM(i,j)*d2psiM(i,j);
	}
	lap -= dot(rv,rv);
	dG[kat] += rv - myG[kat]; myG[kat]=rv;
	dL[kat] += lap -myL[kat]; myL[kat]=lap;
      }

      //not very useful
      CurrentDet *= curRatio;
      curRatio=1.0;
    }

    DiracDeterminantBase::ValueType DiracDeterminantBase::evaluate(ParticleSet& P, PooledData<RealType>& buf) {

      buf.put(psiM.begin(),psiM.end());
      buf.put(FirstAddressOfdV,LastAddressOfdV);
      buf.put(d2psiM.begin(),d2psiM.end());
      buf.put(myL.begin(), myL.end());
      buf.put(FirstAddressOfG,LastAddressOfG);
      buf.put(CurrentDet);

      return CurrentDet;
    }

    DiracDeterminantBase::ValueType
    DiracDeterminantBase::evaluate(ParticleSet& P, 
	     ParticleSet::ParticleGradient_t& G, 
	     ParticleSet::ParticleLaplacian_t& L){

      //evaluateSPO(P, FirstIndex, LastIndex, psiM,dpsiM, d2psiM);
      evaluateAll(P);

      if(NumPtcls==1) {
        CurrentDet=psiM(0,0);
        ValueType y=1.0/CurrentDet;
        psiM(0,0)=y;
        PosType rv = y*dpsiM(0,0);
        G(FirstIndex) += rv;
        L(FirstIndex) += y*d2psiM(0,0) - dot(rv,rv);
      } else {
        CurrentDet = Invert(psiM.data(),NumPtcls,NumOrbitals, WorkSpace.data(), Pivot.data());
        //CurrentDet = Invert(psiM.data(),NumPtcls,NumOrbitals);
        int iat = FirstIndex; //the index of the particle with respect to P
        for(int i=0; i<NumPtcls; i++, iat++) {
          PosType rv = psiM(i,0)*dpsiM(i,0);
          ValueType lap=psiM(i,0)*d2psiM(i,0);
          for(int j=1; j<NumOrbitals; j++) {
            rv += psiM(i,j)*dpsiM(i,j);
            lap += psiM(i,j)*d2psiM(i,j);
          }
          G(iat) += rv;
          L(iat) += lap - dot(rv,rv);
        }
      }
      return CurrentDet;
    }

  void DiracDeterminantBase::resizeByWalkers(int nwalkers) {
    //if(psiM_v.size() < nwalkers) {
    //  psiM_v.resize(nwalkers);
    //  dpsiM_v.resize(nwalkers);
    //  d2psiM_v.resize(nwalkers);
    //  for(int iw=0; iw<nwalkers; iw++) psiM_v[iw].resize(NumPtcls, NumOrbitals);
    //  for(int iw=0; iw<nwalkers; iw++) dpsiM_v[iw].resize(NumPtcls, NumOrbitals);
    //  for(int iw=0; iw<nwalkers; iw++) d2psiM_v[iw].resize(NumPtcls, NumOrbitals);
    //}
    app_error() << "::resizeByWalkers obsolete" << endl;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
