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
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-

#include "QMCWaveFunctions/Fermion/MultiDiracDeterminantBase.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/OhmmsBlas.h"
#include "Numerics/MatrixOperators.h"

namespace qmcplusplus {

  /** set the index of the first particle in the determinant and reset the size of the determinant
   *@param first index of first particle
   *@param nel number of particles in the determinant
   */
  void MultiDiracDeterminantBase::set(int first, int nel) {
    FirstIndex = first;
    resize(nel,nel);
  }



  void MultiDiracDeterminantBase::set_Multi(int first, int nel,int norb) {
    FirstIndex = first;
    resize(nel,norb);
  }


  MultiDiracDeterminantBase::RealType MultiDiracDeterminantBase::updateBuffer(ParticleSet& P, 
      PooledData<RealType>& buf, bool fromscratch) 
  {
    ///needs to be fixed to work with tranposed psiM
    myG=0.0;
    myL=0.0;

    if(fromscratch)
      LogValue=evaluateLog(P,myG,myL);
    else
    {
      if(UpdateMode == ORB_PBYP_RATIO) 
	Phi->evaluate(P, FirstIndex, LastIndex, psiM_temp,dpsiM, d2psiM);    

      if(NumPtcls==1) {
	APP_ABORT("Evaluate Log with 1 particle in MultiMultiDiracDeterminantBase is potentially dangerous");
        ValueType y=1.0/psiM_temp(0,0);
        psiM(0,0)=y;
        GradType rv = y*dpsiM(0,0);
        myG(FirstIndex) += rv;
        myL(FirstIndex) += y*d2psiM(0,0) - dot(rv,rv);
      } else {
        const ValueType* restrict yptr=psiM.data();
        const ValueType* restrict d2yptr=d2psiM.data();
        const GradType* restrict dyptr=dpsiM.data();
        for(int i=0, iat=FirstIndex; i<NumPtcls; i++, iat++) 
        {
          GradType rv;
          ValueType lap=0.0;
          for(int j=0; j<NumOrbitals; j++,yptr++) {
            rv += *yptr * *dyptr++;
            lap += *yptr * *d2yptr++;
          }
          myG(iat) += rv;
          myL(iat) += lap - dot(rv,rv);
        }
      }
    }

    P.G += myG;
    P.L += myL;

    buf.put(psiM.first_address(),psiM.last_address());
    buf.put(psiM_actual.first_address(),psiM_actual.last_address());
    buf.put(FirstAddressOfdV,LastAddressOfdV);
    buf.put(d2psiM.first_address(),d2psiM.last_address());
    buf.put(myL.first_address(), myL.last_address());
    buf.put(FirstAddressOfG,LastAddressOfG);
    buf.put(LogValue);
    buf.put(PhaseValue);

    return LogValue;
  }

  void MultiDiracDeterminantBase::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {

    buf.get(psiM.first_address(),psiM.last_address());
    buf.get(psiM_actual.first_address(),psiM_actual.last_address());
    buf.get(FirstAddressOfdV,LastAddressOfdV);
    buf.get(d2psiM.first_address(),d2psiM.last_address());
    buf.get(myL.first_address(), myL.last_address());
    buf.get(FirstAddressOfG,LastAddressOfG);
    buf.get(LogValue);
    buf.get(PhaseValue);

    //re-evaluate it for testing
    //Phi.evaluate(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM);
    //CurrentDet = Invert(psiM.data(),NumPtcls,NumOrbitals);
    //need extra copy for gradient/laplacian calculations without updating it
    psiM_temp = psiM;
    dpsiM_temp = dpsiM;
    d2psiM_temp = d2psiM;
  }
  


  MultiDiracDeterminantBase::GradType 
    MultiDiracDeterminantBase::evalGrad(ParticleSet& P, int iat)
  {
    const ValueType* restrict yptr=psiM[iat-FirstIndex];
    const GradType* restrict dyptr=dpsiM[iat-FirstIndex];
    GradType rv;
    for(int j=0; j<NumOrbitals; ++j) rv += (*yptr++) *(*dyptr++);
    return rv;
  }
    MultiDiracDeterminantBase::ValueType 
      MultiDiracDeterminantBase::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
  {
    Phi->evaluate(P, iat, psiV, dpsiV, d2psiV);
    WorkingIndex = iat-FirstIndex;

    UpdateMode=ORB_PBYP_PARTIAL;
    
    const ValueType* restrict vptr = psiV.data();
    const ValueType* restrict yptr = psiM[WorkingIndex];
    const GradType* restrict dyptr = dpsiV.data();
    
    GradType rv;
    curRatio = 0.0;
    for(int j=0; j<NumOrbitals; ++j) {
      rv       += yptr[j] * dyptr[j];
      curRatio += yptr[j] * vptr[j];
    }
    grad_iat += (1.0/curRatio) * rv;
    return curRatio;
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
  MultiDiracDeterminantBase::ValueType MultiDiracDeterminantBase::ratio(ParticleSet& P, int iat,
      ParticleSet::ParticleGradient_t& dG, 
      ParticleSet::ParticleLaplacian_t& dL) 
  {
    APP_ABORT("  MultiDiracDeterminantBase::ratio with grad is distabled");
    UpdateMode=ORB_PBYP_ALL;
    Phi->evaluate(P, iat, psiV, dpsiV, d2psiV);
    WorkingIndex = iat-FirstIndex;

    //psiM_temp = psiM;
#ifdef DIRAC_USE_BLAS
    curRatio = BLAS::dot(NumOrbitals,psiM_temp[WorkingIndex],&psiV[0]);
#else
    curRatio= DetRatio(psiM_temp, psiV.begin(),WorkingIndex);
#endif

    if(abs(curRatio)<numeric_limits<RealType>::epsilon()) 
    {
      UpdateMode=ORB_PBYP_RATIO; //singularity! do not update inverse 
      return 0.0;
    }

    //update psiM_temp with the row substituted
    DetUpdate(psiM_temp,psiV,workV1,workV2,WorkingIndex,curRatio);

    //update dpsiM_temp and d2psiM_temp 
    for(int j=0; j<NumOrbitals; j++) {
      dpsiM_temp(WorkingIndex,j)=dpsiV[j];
      d2psiM_temp(WorkingIndex,j)=d2psiV[j];
    }

    int kat=FirstIndex;

    const ValueType* restrict yptr=psiM_temp.data();
    const ValueType* restrict d2yptr=d2psiM_temp.data();
    const GradType* restrict dyptr=dpsiM_temp.data();
    for(int i=0; i<NumPtcls; i++,kat++) {
      //This mimics gemm with loop optimization
      GradType rv;
      ValueType lap=0.0;
      for(int j=0; j<NumOrbitals; j++,yptr++) {
        rv += *yptr * *dyptr++;
        lap += *yptr * *d2yptr++;
      }

      //using inline dot functions
      //GradType rv=dot(psiM_temp[i],dpsiM_temp[i],NumOrbitals);
      //ValueType lap=dot(psiM_temp[i],d2psiM_temp[i],NumOrbitals);

      //Old index: This is not pretty
      //GradType rv =psiM_temp(i,0)*dpsiM_temp(i,0);
      //ValueType lap=psiM_temp(i,0)*d2psiM_temp(i,0);
      //for(int j=1; j<NumOrbitals; j++) {
      //  rv += psiM_temp(i,j)*dpsiM_temp(i,j);
      //  lap += psiM_temp(i,j)*d2psiM_temp(i,j);
      //}
      lap -= dot(rv,rv);
      dG[kat] += rv - myG[kat];  myG_temp[kat]=rv;
      dL[kat] += lap -myL[kat];  myL_temp[kat]=lap;
    }

    return curRatio;
  }


  /** return the ratio only for the  iat-th partcle move
   * @param P current configuration
   * @param iat the particle thas is being moved
   * psiM and psiM_temp shoudl be the same upon entering
   */
  MultiDiracDeterminantBase::ValueType MultiDiracDeterminantBase::ratio(ParticleSet& P, int iat)
  {
    cerr<<"Calling multi-dirac-determinat-base-ratio"<<" "<<psiMInv.size()<<" "<<psiM_temp.size()<<endl;
    UpdateMode=ORB_PBYP_ALL;
    WorkingIndex = iat-FirstIndex;
    cerr<<"Working Index is "<<WorkingIndex<<endl;
    cerr<<"pre-a"<<psiMInv.size()<<" "<<psiM_temp.size()<<endl;
    Phi->evaluate(P, iat, psiV);
    cerr<<"also pre-a"<<" "<<psiMInv.size()<<" "<<psiM_temp.size()<<endl;


	///Copying ground state data to psiM
	for (int orbital=0;orbital<NumOrbitals;orbital++)
	  for (int ptcl=0;ptcl<NumPtcls;ptcl++) 
	    psiM_temp(orbital,ptcl)=psiM_actual(orbital,ptcl); //Note: Because psiM is typically an inverse, it's typically ptcl,orbital
	//done copying
	LogValue=InvertWithLog(psiM_temp.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);


    psiMInv=psiM_temp;
    cerr<<"and pre-tranpose"<<endl;
    MatrixOperators::transpose(psiMInv);
    cerr<<"A"<<endl;
    Excitations.BuildDotProducts(psiMInv,psiM_actual);
    ValueType oldVal=1.0; //eventually store previous valuess
    int coefIndex=0;
    Excitations.CalcSingleExcitations(coefs,oldVal,coefIndex);
    //    Excitations.CalcDoubleExcitations(coefs,oldVal,coefIndex)
    cerr<<"B"<<endl;
    ValueType gs_ratio=0.0;
    for (int orbital=0;orbital<NumOrbitals;orbital++)
      gs_ratio+=psiM_temp(WorkingIndex,orbital)*psiV(orbital);

    cerr<<"ratio gs_ratio is "<<gs_ratio<<endl;


    //relying here on the fact that DetUpdate only uses the first (psiM_temp.cols()) rows of psiV!
    //update psiM_temp with the row substituted
    cerr<<"C"<<endl;

    for (int orbital=0;orbital<NumOrbitals_total;orbital++){
      psiV_old(orbital)=psiM_actual(orbital,WorkingIndex);
      psiM_actual(orbital,WorkingIndex)=psiV(orbital);
    }

    DetUpdate(psiM_temp,psiV,workV1,workV2,WorkingIndex,gs_ratio);


	///Copying ground state data to psiM
	for (int orbital=0;orbital<NumOrbitals;orbital++)
	  for (int ptcl=0;ptcl<NumPtcls;ptcl++) 
	    psiM_temp(orbital,ptcl)=psiM_actual(orbital,ptcl); //Note: Because psiM is typically an inverse, it's typically ptcl,orbital
	//done copying
	LogValue=InvertWithLog(psiM_temp.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);



    psiMInv=psiM_temp;
    MatrixOperators::transpose(psiMInv);
    cerr<<"D"<<endl;
   //currnetly psiM sending incorrectly transposed!
    Excitations.BuildDotProducts(psiMInv,psiM_actual);
    ValueType val=1.0;
    coefIndex=0;
    cerr<<"E"<<endl;
    Excitations.CalcSingleExcitations(coefs,val,coefIndex);
    //    Excitations.CalcDoubleExcitations(coefs,val,coefIndex);
    curRatio=(gs_ratio*val)/oldVal;
    cerr<<"curr Ratio is "<<curRatio<<" "<<psiMInv.size()<<" "<<psiM_temp.size()<<endl;
    return curRatio;
  }





  /** move was accepted, update the real container
  */
  void MultiDiracDeterminantBase::acceptMove(ParticleSet& P, int iat) 
  {
    cerr<<"MultiDirac accepting"<<endl;
    PhaseValue += evaluatePhase(curRatio);
    LogValue +=std::log(std::abs(curRatio));
    switch(UpdateMode)
    {
      case ORB_PBYP_RATIO:
        DetUpdate(psiM,psiV,workV1,workV2,WorkingIndex,curRatio);
        break;
      case ORB_PBYP_PARTIAL:
        //psiM = psiM_temp;
	DetUpdate(psiM,psiV,workV1,workV2,WorkingIndex,curRatio);
        std::copy(dpsiV.begin(),dpsiV.end(),dpsiM[WorkingIndex]);
        std::copy(d2psiV.begin(),d2psiV.end(),d2psiM[WorkingIndex]);

        //////////////////////////////////////
        ////THIS WILL BE REMOVED. ONLY FOR DEBUG DUE TO WAVEFUNCTIONTEST
        //myG = myG_temp;
        //myL = myL_temp;
        ///////////////////////

        break;
      default:
	cerr<<"Accepting move"<<endl;
        myG = myG_temp;
        myL = myL_temp;
        psiM = psiM_temp;
        std::copy(dpsiV.begin(),dpsiV.end(),dpsiM[WorkingIndex]);
        std::copy(d2psiV.begin(),d2psiV.end(),d2psiM[WorkingIndex]);
	cerr<<"Done Accepting move"<<endl;
        break;
    }

    curRatio=1.0;
  }

  /** move was rejected. copy the real container to the temporary to move on
  */
  void MultiDiracDeterminantBase::restore(int iat) {
    cerr<<"INto restore"<<endl;
    if(UpdateMode == ORB_PBYP_ALL) {
      cerr<<"Rejecting move"<<endl;
      psiM_temp = psiM;
      for (int orbital=0;orbital<NumOrbitals_total;orbital++)
	psiM_actual(orbital,WorkingIndex)=psiV_old(orbital);
      std::copy(dpsiM[WorkingIndex],dpsiM[WorkingIndex+1],dpsiM_temp[WorkingIndex]);
      std::copy(d2psiM[WorkingIndex],d2psiM[WorkingIndex+1],d2psiM_temp[WorkingIndex]);
      cerr<<"Done Rejecting move"<<endl;
    }
    else if (UpdateMode== ORB_PBYP_RATIO){
      for (int orbital=0;orbital<NumOrbitals_total;orbital++)
	psiM_actual(orbital,WorkingIndex)=psiV_old(orbital);
    }
    curRatio=1.0;
  }


  MultiDiracDeterminantBase::RealType
    MultiDiracDeterminantBase::evaluateLog(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& G, 
        ParticleSet::ParticleLaplacian_t& L)
    {
      cerr<<"I'm calling evaluate log"<<endl;
      //      Phi->evaluate(P, FirstIndex, LastIndex, psiM,dpsiM, d2psiM);
      Phi->evaluate(P, FirstIndex, LastIndex, psiM_actual,dpsiM, d2psiM);

      if(NumPtcls==1) 
      {
	APP_ABORT("Evaluate Log with 1 particle in MultiMultiDiracDeterminantBase is potentially dangerous");
        //CurrentDet=psiM(0,0);
        ValueType det=psiM(0,0);
        ValueType y=1.0/det;
        psiM(0,0)=y;
        GradType rv = y*dpsiM(0,0);
        G(FirstIndex) += rv;
        L(FirstIndex) += y*d2psiM(0,0) - dot(rv,rv);
        LogValue = evaluateLogAndPhase(det,PhaseValue);
      } else {
	
	cerr<<"WHEN COPYING NUMPTCLS IS "<<NumPtcls<<endl;
	///Copying ground state data to psiM
	for (int orbital=0;orbital<NumOrbitals;orbital++)
	  for (int ptcl=0;ptcl<NumPtcls;ptcl++) 
	    psiM(orbital,ptcl)=psiM_actual(orbital,ptcl); //Note: Because psiM is typically an inverse, it's typically ptcl,orbital
	//done copying
	LogValue=InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
        const ValueType* restrict yptr=psiM.data();
        const ValueType* restrict d2yptr=d2psiM.data();
        const GradType* restrict dyptr=dpsiM.data();
        for(int i=0, iat=FirstIndex; i<NumPtcls; i++, iat++) {
          GradType rv;
          ValueType lap=0.0;
          for(int j=0; j<NumOrbitals; j++,yptr++) {
            rv += *yptr * *dyptr++;
            lap += *yptr * *d2yptr++;
          }
          G(iat) += rv;
          L(iat) += lap - dot(rv,rv);
        }
      }
      psiM_temp = psiM;
      cerr<<"EVAL LOG gs_ratio is "<<std::exp(LogValue)<<endl;
      double val1=std::exp(LogValue)*std::cos(abs(PhaseValue));
      for (int orb=NumOrbitals;orb<NumOrbitals_total;orb++){
	for (int orbital=0;orbital<NumOrbitals;orbital++){
	  for (int ptcl=0;ptcl<NumPtcls;ptcl++){ 
	    psiM(orbital,ptcl)=psiM_actual(orbital,ptcl);
	  }
	}
	for (int ptcl=0;ptcl<NumPtcls;ptcl++){
	  cerr<<"reasonable "<<psiM_actual(orb,ptcl)<<endl;
	  psiM(NumOrbitals-1,ptcl)=psiM_actual(orb,ptcl);
	}
	double PhaseValuep;
	double LogValuep=InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValuep);
	cerr<<"Pre-LOGVALUE PHASEVALUE "<<LogValue<<" "<<PhaseValue<<" "<<LogValuep<<" "<<PhaseValuep<<endl;

	double val2=std::exp(LogValuep)*std::cos(abs(PhaseValuep));
	cerr<<"RATIO: "<<val2/val1<<endl;
	double val=std::exp(LogValue)*std::cos(abs(PhaseValue))+std::exp(LogValuep)*std::cos(abs(PhaseValuep));
	LogValue=std::log(abs(val));
	if (val>0)
	  PhaseValue=0.0;
	else 
	  PhaseValue=M_PI;
	cerr<<"Post-LOGVALUE PHASEVALUE "<<LogValue<<" "<<PhaseValue<<" "<<LogValuep<<" "<<PhaseValuep<<endl;
	psiM=psiM_temp;
      }
      cerr<<"LOGVALUE PHASEVALUE "<<LogValue<<" "<<PhaseValue<<endl;
      return LogValue;
    }


#include "MultiDiracDeterminantBase_help.cpp";
}
/***************************************************************************
 * $RCSfile$   $Author: kesler $
 * $Revision: 3582 $   $Date: 2009-02-20 18:07:41 -0600 (Fri, 20 Feb 2009) $
 * $Id: MultiDiracDeterminantBase.cpp 3582 2009-02-21 00:07:41Z kesler $ 
 ***************************************************************************/
