//////////////////////////////////////////////////////////////////
// (c) Copyright 2006- by Jeongnim Kim
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCWaveFunctions/AGPDeterminant.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/OhmmsBlas.h"

namespace qmcplusplus {

  AGPDeterminant::AGPDeterminant(BasisSetType* bs): BasisSet(bs), NumPtcls(0){
  }
  AGPDeterminant::~AGPDeterminant() {}

  void AGPDeterminant::resize(int nup, int ndown) {

    if(NumPtcls == 0) { //use Numptcls to ensure resizing only once
      BasisSize=BasisSet->size();
      Coeff.resize(BasisSize,BasisSize);

      Nup=nup;
      Ndown=ndown;
      NumPtcls=nup+ndown;

      psiM.resize(nup,nup);
      dpsiM.resize(nup,nup);
      d2psiM.resize(nup,nup);
      psiM_temp.resize(nup,nup);
      dpsiM_temp.resize(nup,nup);
      d2psiM_temp.resize(nup,nup);
      psiMinv.resize(nup,nup);
      psiV.resize(nup);
      WorkSpace.resize(nup);
      Pivot.resize(nup);
      phiU.resize(nup,BasisSize);
      phiD.resize(ndown,BasisSize);
    }
  }

  /** Calculate the value of the Dirac determinant for particles
   *@param P input configuration containing N particles
   *@param G a vector containing N gradients
   *@param L a vector containing N laplacians
   *@return the value of the determinant
   *
   *\f$ (first,first+nel). \f$  Add the gradient and laplacian 
   *contribution of the determinant to G(radient) and L(aplacian)
   *for local energy calculations.
   */ 
  AGPDeterminant::ValueType
  AGPDeterminant::evaluate(ParticleSet& P, 
      ParticleSet::ParticleGradient_t& G, 
      ParticleSet::ParticleLaplacian_t& L){

    BasisSet->evaluate(P);

    /* evaluate psi_up(iat)= \sum_{j} C_{ij} \phi_j^{u}(r_{iat}) 
     * psi_down(iat-Nup) =  \sum_{j} C_{ij} \phi_j^{d}(r_{iat})
     */
    //rewrite it by gemm
    for(int u=0; u<Nup; u++) {
      BLAS::gemv(BasisSize,BasisSize,Coeff.data(),BasisSet->y(u),phiU[u]);
    }

    for(int d=0,jat=Nup; d<Ndown; d++,jat++) {
      BLAS::gemv(BasisSize,BasisSize,Coeff.data(),BasisSet->y(jat),phiD[d]);
    }

    for(int u=0; u<Nup; u++) {
      for(int d=0, jat=Nup; d<Ndown; d++,jat++) {
        psiM(d,u) = BLAS::dot(BasisSize,phiU[u],BasisSet->y(jat));
      }
    }

    CurrentDet = Invert(psiM.data(),Nup,Nup,WorkSpace.data(),Pivot.data());

    for(int iat=0; iat<Nup; iat++) {
      GradType rv;
      ValueType lap=0;
      for(int d=0, jat=Nup; d<Ndown; d++,jat++) {
        ValueType dfac=psiM(iat,d);
        for(int k=0; k<BasisSize; k++) {
          rv += dfac*BasisSet->dY(iat,k)*phiD(d,k);
          lap += dfac*BasisSet->d2Y(iat,k)*phiD(d,k);
        }
      }
      G(iat) += rv;
      L(iat) += lap-dot(rv,rv);
    }

    for(int jat=Nup,d=0; jat<NumPtcls; jat++,d++) {
      GradType rv;
      ValueType lap=0;
      for(int u=0; u<Nup; u++) {
        ValueType dfac=psiM(d,u); //transpose
        for(int k=0; k<BasisSize; k++) {
          rv += dfac*BasisSet->dY(jat,k)*phiU(u,k);
          lap += dfac*BasisSet->d2Y(jat,k)*phiU(u,k);
        }
      }
      G(jat) += rv;
      L(jat) += lap-dot(rv,rv);
    }

    return CurrentDet;
  }

  AGPDeterminant::ValueType 
  AGPDeterminant::registerData(ParticleSet& P, PooledData<RealType>& buf) {
//
//      if(NP == 0) {//first time, allocate once
//	//int norb = cols();
//	dpsiV.resize(NumOrbitals);
//	d2psiV.resize(NumOrbitals);
//	workV1.resize(NumOrbitals);
//	workV2.resize(NumOrbitals);
//	NP=P.getTotalNum();
//	myG.resize(NP);
//	myL.resize(NP);
//	myG_temp.resize(NP);
//	myL_temp.resize(NP);
//	FirstAddressOfG = &myG[0][0];
//	LastAddressOfG = FirstAddressOfG + NP*DIM;
//	FirstAddressOfdV = &(dpsiM(0,0)[0]); //(*dpsiM.begin())[0]);
//	LastAddressOfdV = FirstAddressOfdV + NumPtcls*NumOrbitals*DIM;
//      }
//
//      //allocate once but each walker calls this
//      myG=0.0;
//      myL=0.0;
//
//      ValueType x=evaluate(P,myG,myL); 
//
//      P.G += myG;
//      P.L += myL;
//
//      //add the data: determinant, inverse, gradient and laplacians
//      buf.add(psiM.begin(),psiM.end());
//      buf.add(FirstAddressOfdV,LastAddressOfdV);
//      buf.add(d2psiM.begin(),d2psiM.end());
//      buf.add(myL.begin(), myL.end());
//      buf.add(FirstAddressOfG,LastAddressOfG);
//      buf.add(CurrentDet);
//
      return CurrentDet;
    }

  AGPDeterminant::ValueType 
  AGPDeterminant::updateBuffer(ParticleSet& P, PooledData<RealType>& buf) {
//      myG=0.0;
//      myL=0.0;
//      ValueType x=evaluate(P,myG,myL); 
//      P.G += myG;
//      P.L += myL;
//      buf.put(psiM.begin(),psiM.end());
//      buf.put(FirstAddressOfdV,LastAddressOfdV);
//      buf.put(d2psiM.begin(),d2psiM.end());
//      buf.put(myL.begin(), myL.end());
//      buf.put(FirstAddressOfG,LastAddressOfG);
//      buf.put(CurrentDet);

      return CurrentDet;
    }

    void AGPDeterminant::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {

//      buf.get(psiM.begin(),psiM.end());
//      buf.get(FirstAddressOfdV,LastAddressOfdV);
//      buf.get(d2psiM.begin(),d2psiM.end());
//      buf.get(myL.begin(), myL.end());
//      buf.get(FirstAddressOfG,LastAddressOfG);
//      buf.get(CurrentDet);
//
//      //re-evaluate it for testing
//      //Phi.evaluate(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM);
//      //CurrentDet = Invert(psiM.data(),NumPtcls,NumOrbitals);
//      //need extra copy for gradient/laplacian calculations without updating it
//      psiM_temp = psiM;
//      dpsiM_temp = dpsiM;
//      d2psiM_temp = d2psiM;
    }


    /** return the ratio only for the  iat-th partcle move
     * @param P current configuration
     * @param iat the particle thas is being moved
     */
    AGPDeterminant::ValueType AGPDeterminant::ratio(ParticleSet& P, int iat) {
      BasisSet->evaluate(P,iat);
      return 1.0;

//      ValueType* yptr=BasisSet->y(0);
//      if(iat<Nup) {
//        for(int i=0; i<Ndown; i++) {
//          psiV[i]=BLAS::dot(BasisSize,phi_down(i),yptr);
//        }
//        //add the non-paired phi
//      } else {
//        for(int i=0; i<Ndown; i++) {
//          psiV[i]=BLAS::dot(BasisSize,phi_up(i),yptr);
//        }
//        iat-=Nup;
//      }
//
//#ifdef DIRAC_USE_BLAS
//      return curRatio = BLAS::dot(Nup,psiM[iat],&psiV[0]);
//#else
//      return DetRatio(psiM, psiV.begin(),iat);
//#endif
    }

    /** return the ratio
     * @param P current configuration
     * @param iat particle whose position is moved
     * @param dG differential Gradients
     * @param dL differential Laplacians
     *
     * Data member *_temp contain the data assuming that the move is accepted
     * and are used to evaluate differential Gradients and Laplacians.
     */
    AGPDeterminant::ValueType AGPDeterminant::ratio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG, 
		    ParticleSet::ParticleLaplacian_t& dL) {

//      Phi.evaluate(P, iat, psiV, dpsiV, d2psiV);
//      WorkingIndex = iat-FirstIndex;
//
//#ifdef DIRAC_USE_BLAS
//      curRatio = BLAS::dot(NumOrbitals,psiM_temp[WorkingIndex],&psiV[0]);
//#else
//      curRatio= DetRatio(psiM_temp, psiV.begin(),WorkingIndex);
//#endif
//
//      //update psiM_temp with the row substituted
//      DetUpdate(psiM_temp,psiV,workV1,workV2,WorkingIndex,curRatio);
//
//      //update dpsiM_temp and d2psiM_temp 
//      for(int j=0; j<NumOrbitals; j++) {
//	dpsiM_temp(WorkingIndex,j)=dpsiV[j];
//	d2psiM_temp(WorkingIndex,j)=d2psiV[j];
//      }
//
//      int kat=FirstIndex;
//      for(int i=0; i<NumPtcls; i++,kat++) {
//	PosType rv =psiM_temp(i,0)*dpsiM_temp(i,0);
//	ValueType lap=psiM_temp(i,0)*d2psiM_temp(i,0);
//        for(int j=1; j<NumOrbitals; j++) {
//	  rv += psiM_temp(i,j)*dpsiM_temp(i,j);
//	  lap += psiM_temp(i,j)*d2psiM_temp(i,j);
//	}
//	lap -= dot(rv,rv);
//	dG[kat] += rv - myG[kat];  myG_temp[kat]=rv;
//	dL[kat] += lap -myL[kat];  myL_temp[kat]=lap;
//      }
//
      return curRatio;
    }


    /** move was accepted, update the real container
     */
    void AGPDeterminant::update(ParticleSet& P, int iat) {
//      CurrentDet *= curRatio;
//      myG = myG_temp;
//      myL = myL_temp;
//      psiM = psiM_temp;
//      for(int j=0; j<NumOrbitals; j++) {
//	dpsiM(WorkingIndex,j)=dpsiV[j];
//	d2psiM(WorkingIndex,j)=d2psiV[j];
//      }
      curRatio=1.0;
    }

    /** move was rejected. copy the real container to the temporary to move on
     */
    void AGPDeterminant::restore(int iat) {
//      psiM_temp = psiM;
//      for(int j=0; j<NumOrbitals; j++) {
//	dpsiM_temp(WorkingIndex,j)=dpsiM(WorkingIndex,j);
//	d2psiM_temp(WorkingIndex,j)=d2psiM(WorkingIndex,j);
//      }
//      curRatio=1.0;
    }

    
    void AGPDeterminant::update(ParticleSet& P, 
		ParticleSet::ParticleGradient_t& dG, 
		ParticleSet::ParticleLaplacian_t& dL,
		int iat) {

//      DetUpdate(psiM,psiV,workV1,workV2,WorkingIndex,curRatio);
//      for(int j=0; j<NumOrbitals; j++) {
//	dpsiM(WorkingIndex,j)=dpsiV[j];
//	d2psiM(WorkingIndex,j)=d2psiV[j];
//      }
//
//      int kat=FirstIndex;
//      for(int i=0; i<NumPtcls; i++,kat++) {
//	PosType rv =psiM(i,0)*dpsiM(i,0);
//	ValueType lap=psiM(i,0)*d2psiM(i,0);
//	for(int j=1; j<NumOrbitals; j++) {
//	  rv += psiM(i,j)*dpsiM(i,j);
//	  lap += psiM(i,j)*d2psiM(i,j);
//	}
//	lap -= dot(rv,rv);
//	dG[kat] += rv - myG[kat]; myG[kat]=rv;
//	dL[kat] += lap -myL[kat]; myL[kat]=lap;
//      }
//
//      //not very useful
//      CurrentDet *= curRatio;
//      curRatio=1.0;
//    }
//
//    ValueType evaluate(ParticleSet& P, PooledData<RealType>& buf) {
//
//      buf.put(psiM.begin(),psiM.end());
//      buf.put(FirstAddressOfdV,LastAddressOfdV);
//      buf.put(d2psiM.begin(),d2psiM.end());
//      buf.put(myL.begin(), myL.end());
//      buf.put(FirstAddressOfG,LastAddressOfG);
//      buf.put(CurrentDet);
//
    }

  AGPDeterminant::ValueType 
  AGPDeterminant::evaluate(ParticleSet& P, PooledData<RealType>& buf) {
    return CurrentDet;
  }

    void AGPDeterminant::resizeByWalkers(int nwalkers) {
      //don't know what to do here
    }


}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
