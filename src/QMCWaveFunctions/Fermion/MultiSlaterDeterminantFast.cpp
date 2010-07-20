//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#include "QMCWaveFunctions/Fermion/MultiSlaterDeterminantFast.h"
#include "QMCWaveFunctions/Fermion/MultiDiracDeterminantBase.h"
#include "ParticleBase/ParticleAttribOps.h"

namespace qmcplusplus {

  MultiSlaterDeterminantFast::MultiSlaterDeterminantFast(ParticleSet& targetPtcl, MultiDiracDeterminantBase* up, MultiDiracDeterminantBase* dn):
    RatioTimer("MultiSlaterDeterminantFast::ratio"),
    RatioGradTimer("MultiSlaterDeterminantFast::ratioGrad"),
    RatioAllTimer("MultiSlaterDeterminantFast::ratio(all)"),
    UpdateTimer("MultiSlaterDeterminantFast::updateBuffer"),
    EvaluateTimer("MultiSlaterDeterminantFast::evaluate")
  { 
    registerTimers();
    //Optimizable=true;
    Optimizable=true;
    OrbitalName="MultiSlaterDeterminantFast";

    NP = targetPtcl.getTotalNum();
    nels_up = targetPtcl.last(0)-targetPtcl.first(0);
    nels_dn = targetPtcl.last(1)-targetPtcl.first(1);
    FirstIndex_up=targetPtcl.first(0);
    FirstIndex_dn=targetPtcl.first(1);
    Dets.resize(2);
    Dets[0]=up;
    Dets[1]=dn;
    myG.resize(NP);
    myL.resize(NP);
    myG_temp.resize(NP);
    myL_temp.resize(NP);
    DetID.resize(NP);
    for(int i=0; i<targetPtcl.groups(); ++i)
      for(int j=targetPtcl.first(i); j<targetPtcl.last(i); ++j) DetID[j]=i;    
  }
  
  OrbitalBasePtr MultiSlaterDeterminantFast::makeClone(ParticleSet& tqp) const
  { 
    MultiSlaterDeterminantFast* clone = new MultiSlaterDeterminantFast(tqp,(MultiDiracDeterminantBase*) Dets[0]->makeClone(tqp),(MultiDiracDeterminantBase*) Dets[1]->makeClone(tqp)); 
    clone->Dets[0]->resetTargetParticleSet(tqp);
    clone->Dets[1]->resetTargetParticleSet(tqp);
    clone->C2node_up=C2node_up;
    clone->C2node_dn=C2node_dn;
    clone->Optimizable=Optimizable;
    clone->C=C;
    clone->myVars=myVars;
    
    return clone;
  }
  
  MultiSlaterDeterminantFast::~MultiSlaterDeterminantFast() { }

  void MultiSlaterDeterminantFast::resetTargetParticleSet(ParticleSet& P) 
  {
    for(int i=0; i<Dets.size(); i++)
      Dets[i]->resetTargetParticleSet(P);
  }

//  void MultiSlaterDeterminantFast::resize(int n1, int n2)
//  {
//  }

  void MultiSlaterDeterminantFast::testMSD(ParticleSet& P, int iat)
  {

//     APP_ABORT("Testing disabled for safety");
     app_log() <<"Testing MSDFast. \n";

     int n = nels_up+nels_dn;
     ParticleSet::ParticleGradient_t G(n),G0(n);
     ParticleSet::ParticleLaplacian_t L(n),L0(n);
     ValueType log, log0;
     log = msd->evaluate(P,G,L);
     log0 = evaluate(P,G0,L0);

/*
     app_log() <<"Testing evaluate(P,G,L). \n";
     cout<<endl <<endl;
     cout<<"Psi: " <<log <<"   " <<log0 <<"   " <<log/log0 <<endl;

     for(int i=0; i<n; i++) {
       cout<<i  <<"\n"
           <<"  x: " <<G(i)[0]-G0(i)[0] <<"\n"
           <<"  y: " <<G(i)[1]-G0(i)[1] <<"\n"
           <<"  z: " <<G(i)[2]-G0(i)[2] <<"\n"
           <<"  d2: " <<L(i)-L0(i) <<"\n"
           <<endl;
     } 
     cout<<endl <<endl;
     APP_ABORT("end of test 1");
*/

     Walker_t::Buffer_t wbuffer;
     wbuffer.clear();
     log=registerData(P,wbuffer);

     log = msd->evaluate(P,G,L);
     log0 = evaluate(P,G0,L0);

     PosType dr;
     dr[0] = 0.1; dr[1]=0.05; dr[2] = -0.01;
     PosType newpos(P.makeMove(iat,dr));

     app_log() <<"Testing ratio(P,dG,dL). \n";
     G=0;G0=0;L=0;L0=0;
     log = msd->ratio(P,iat,G,L);
     log0 = ratio(P,iat,G0,L0);
     cout<<"Psi: " <<log <<"   " <<log0 <<"   " <<log/log0 <<endl;

     for(int i=0; i<n; i++) {
       cout<<i  <<"\n"
           <<"  x: " <<G(i)[0]-G0(i)[0] <<"  " <<G(i)[0]   <<"\n"
           <<"  y: " <<G(i)[1]-G0(i)[1] <<"  " <<G(i)[1] <<"\n"
           <<"  z: " <<G(i)[2]-G0(i)[2] <<"  " <<G(i)[2] <<"\n"
           <<"  d2: " <<L(i)-L0(i) <<"  " <<L(i) <<"\n"
           <<endl;
     }

     cout<<endl <<endl;
     APP_ABORT("After MultiSlaterDeterminantFast::testMSD()");

  }

  OrbitalBase::ValueType MultiSlaterDeterminantFast::evaluate(ParticleSet& P
      , ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
  { 
    EvaluateTimer.start();
    //for(int i=0; i<Dets.size(); i++)
    //  Dets[i]->evaluateForWalkerMove(P);
    Dets[0]->evaluateForWalkerMove(P);
    Dets[1]->evaluateForWalkerMove(P);
    // can this change over time??? I don't know yet
    ValueVector_t& detValues_up = Dets[0]->detValues; 
    ValueVector_t& detValues_dn = Dets[1]->detValues; 
    GradMatrix_t& grads_up = Dets[0]->grads;
    GradMatrix_t& grads_dn = Dets[1]->grads;
    ValueMatrix_t& lapls_up = Dets[0]->lapls;
    ValueMatrix_t& lapls_dn = Dets[1]->lapls;
    int N1 = Dets[0]->FirstIndex;
    int N2 = Dets[1]->FirstIndex;
    int NP1 = Dets[0]->NumPtcls;
    int NP2 = Dets[1]->NumPtcls;

    psiCurrent=0.0;
    myG=0.0;
    myL=0.0;
    vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
    vector<RealType>::iterator it(C.begin()),last(C.end());
    while(it != last) {
      psiCurrent += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
      for(int k=0,n=N1; k<NP1; k++,n++) {
        myG(n) += (*it)*grads_up(*upC,k)*detValues_dn[*dnC]; 
        myL(n) += (*it)*lapls_up(*upC,k)*detValues_dn[*dnC]; 
      }
      for(int k=0,n=N2; k<NP2; k++,n++) {
        myG(n) += (*it)*grads_dn(*dnC,k)*detValues_up[*upC];
        myL(n) += (*it)*lapls_dn(*dnC,k)*detValues_up[*upC]; 
      }
      it++;upC++;dnC++;
    }
    ValueType psiinv = 1.0/psiCurrent;
    myG *= psiinv;
    myL *= psiinv;
    G += myG;
    for(int i=0; i<L.size(); i++) 
      L(i) += myL[i] - dot(myG[i],myG[i]);
    EvaluateTimer.stop();
    return psiCurrent;
  }

  OrbitalBase::RealType MultiSlaterDeterminantFast::evaluateLog(ParticleSet& P
      , ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
  {
    ValueType psi = evaluate(P,G,L);
    return LogValue = evaluateLogAndPhase(psi,PhaseValue);
  }

  OrbitalBase::GradType MultiSlaterDeterminantFast::evalGrad(ParticleSet& P, int iat)
  {
    GradType grad_iat=0.0;
    if(DetID[iat] == 0) {
      ValueVector_t& detValues_up = Dets[0]->detValues;
      ValueVector_t& detValues_dn = Dets[1]->detValues;
      GradMatrix_t& grads_up = Dets[0]->grads;
      int N1 = Dets[0]->FirstIndex;
      ValueType psi=0.0;
      vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
      vector<RealType>::iterator it(C.begin()),last(C.end());
      while(it != last) {
        psi += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
        grad_iat += (*it)*grads_up(*upC,iat-N1)*detValues_dn[*dnC]; 
        it++;upC++;dnC++;
      }
      grad_iat *= 1.0/psi;
      return grad_iat;
    } else {
      ValueType psi=0.0;
      ValueVector_t& detValues_up = Dets[0]->detValues;
      ValueVector_t& detValues_dn = Dets[1]->detValues;
      GradMatrix_t& grads_dn = Dets[1]->grads;
      int N2 = Dets[1]->FirstIndex;
      vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
      vector<RealType>::iterator it(C.begin()),last(C.end());
      while(it != last) {
        psi += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
        grad_iat += (*it)*grads_dn(*dnC,iat-N2)*detValues_up[*upC];
        it++;upC++;dnC++;
      }
      grad_iat *= 1.0/psi;
      return grad_iat;
    }
  }

  OrbitalBase::ValueType MultiSlaterDeterminantFast::ratioGrad(ParticleSet& P
      , int iat, GradType& grad_iat)
  {
    UpdateMode=ORB_PBYP_PARTIAL;
    if(DetID[iat] == 0) {
      RatioGradTimer.start();
      Dets[0]->evaluateDetsAndGradsForPtclMove(P,iat);
      ValueVector_t& detValues_up = Dets[0]->new_detValues;
      ValueVector_t& detValues_dn = Dets[1]->detValues;
      GradMatrix_t& grads_up = Dets[0]->new_grads;
      int N1 = Dets[0]->FirstIndex;
      vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
      vector<RealType>::iterator it(C.begin()),last(C.end());
      ValueType psiNew=0.0;
      GradType dummy=0.0; 
      while(it != last) {
        psiNew += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
        dummy += (*it)*grads_up(*upC,iat-N1)*detValues_dn[*dnC]; 
        it++;upC++;dnC++;
      }
      grad_iat+=dummy/psiNew;
      curRatio = psiNew/psiCurrent;
      RatioGradTimer.stop();
      return curRatio; 
    } else {
      RatioGradTimer.start();
      Dets[1]->evaluateDetsAndGradsForPtclMove(P,iat);

      ValueVector_t& detValues_up = Dets[0]->detValues;
      ValueVector_t& detValues_dn = Dets[1]->new_detValues;
      GradMatrix_t& grads_dn = Dets[1]->new_grads;
      int N2 = Dets[1]->FirstIndex;
      vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
      vector<RealType>::iterator it(C.begin()),last(C.end());
      ValueType psiNew=0.0;
      GradType dummy=0.0;
      while(it != last) {
        psiNew += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
        dummy += (*it)*grads_dn(*dnC,iat-N2)*detValues_up[*upC];
        it++;upC++;dnC++;
      }
      grad_iat+=dummy/psiNew;
      curRatio = psiNew/psiCurrent;
      RatioGradTimer.stop();
      return curRatio;
    }
  }


// This routine need work, sloppy for now
  OrbitalBase::ValueType  MultiSlaterDeterminantFast::ratio(ParticleSet& P, int iat
     , ParticleSet::ParticleGradient_t& dG,ParticleSet::ParticleLaplacian_t& dL)
  {
    UpdateMode=ORB_PBYP_ALL;
    if(DetID[iat] == 0) {
      RatioAllTimer.start();

/*
      P.acceptMove(iat);
      Dets[0]->evaluateForWalkerMove(P);
      ValueVector_t& detValues_up = Dets[0]->detValues;
      ValueVector_t& detValues_dn = Dets[1]->detValues;
      GradMatrix_t& grads_up = Dets[0]->grads;
      GradMatrix_t& grads_dn = Dets[1]->grads;
      ValueMatrix_t& lapls_up = Dets[0]->lapls;
      ValueMatrix_t& lapls_dn = Dets[1]->lapls;
*/

//*
      Dets[0]->evaluateAllForPtclMove(P,iat);

      ValueVector_t& detValues_up = Dets[0]->new_detValues;
      ValueVector_t& detValues_dn = Dets[1]->detValues;
      GradMatrix_t& grads_up = Dets[0]->new_grads;
      GradMatrix_t& grads_dn = Dets[1]->grads;
      ValueMatrix_t& lapls_up = Dets[0]->new_lapls;
      ValueMatrix_t& lapls_dn = Dets[1]->lapls;
//*/

      int N1 = Dets[0]->FirstIndex;
      int N2 = Dets[1]->FirstIndex;
      int NP1 = Dets[0]->NumPtcls;
      int NP2 = Dets[1]->NumPtcls;
      ValueType psiNew=0.0;

      // myG,myL should contain current grad and lapl

      vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
      vector<RealType>::iterator it(C.begin()),last(C.end());
      myG_temp=0.0;
      myL_temp=0.0;
      while(it != last) {
        psiNew += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
        for(int k=0,n=N1; k<NP1; k++,n++) {
          myG_temp(n) += (*it)*grads_up(*upC,k)*detValues_dn[*dnC];
          myL_temp(n) += (*it)*lapls_up(*upC,k)*detValues_dn[*dnC];
        }
        for(int k=0,n=N2; k<NP2; k++,n++) {
          myG_temp(n) += (*it)*grads_dn(*dnC,k)*detValues_up[*upC];
          myL_temp(n) += (*it)*lapls_dn(*dnC,k)*detValues_up[*upC];
        }
        it++;upC++;dnC++;
      }

      ValueType psiNinv=1.0/psiNew;
      myG_temp *= psiNinv;
      myL_temp *= psiNinv;
      dG += myG_temp-myG;
      for(int i=0; i<dL.size(); i++)
        dL(i) += myL_temp[i] - myL[i] - dot(myG_temp[i],myG_temp[i]) + dot(myG[i],myG[i]);

      curRatio = psiNew/psiCurrent;
      RatioAllTimer.stop();
      return curRatio;

    } else {
      RatioAllTimer.start();

      Dets[1]->evaluateAllForPtclMove(P,iat);

      ValueVector_t& detValues_up = Dets[0]->detValues;
      ValueVector_t& detValues_dn = Dets[1]->new_detValues;
      GradMatrix_t& grads_up = Dets[0]->grads;
      GradMatrix_t& grads_dn = Dets[1]->new_grads;
      ValueMatrix_t& lapls_up = Dets[0]->lapls;
      ValueMatrix_t& lapls_dn = Dets[1]->new_lapls;

      int N1 = Dets[0]->FirstIndex;
      int N2 = Dets[1]->FirstIndex;
      int NP1 = Dets[0]->NumPtcls;
      int NP2 = Dets[1]->NumPtcls;
      ValueType psiNew=0.0;

      // myG,myL should contain current grad and lapl

      vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
      vector<RealType>::iterator it(C.begin()),last(C.end());
      myG_temp=0.0;
      myL_temp=0.0;
      while(it != last) {
        psiNew += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
        for(int k=0,n=N1; k<NP1; k++,n++) {
          myG_temp(n) += (*it)*grads_up(*upC,k)*detValues_dn[*dnC];
          myL_temp(n) += (*it)*lapls_up(*upC,k)*detValues_dn[*dnC];
        }
        for(int k=0,n=N2; k<NP2; k++,n++) {
          myG_temp(n) += (*it)*grads_dn(*dnC,k)*detValues_up[*upC];
          myL_temp(n) += (*it)*lapls_dn(*dnC,k)*detValues_up[*upC];
        }
        it++;upC++;dnC++;
      }

      ValueType psiNinv=1.0/psiNew;
      myG_temp *= psiNinv;
      myL_temp *= psiNinv;
      dG += myG_temp-myG;
      for(int i=0; i<dL.size(); i++)
        dL(i) += myL_temp[i] - myL[i] - dot(myG_temp[i],myG_temp[i]) + dot(myG[i],myG[i]);

      curRatio = psiNew/psiCurrent;
      RatioAllTimer.stop();
      return curRatio;
    }
  }

// use ci_node for this routine only
  OrbitalBase::ValueType MultiSlaterDeterminantFast::ratio(ParticleSet& P, int iat)
  {
// debug
//    testMSD(P,iat);
    UpdateMode=ORB_PBYP_RATIO;
    if(DetID[iat] == 0) {
      RatioTimer.start();
      Dets[0]->evaluateDetsForPtclMove(P,iat);
      ValueVector_t& detValues_up = Dets[0]->new_detValues;
      ValueVector_t& detValues_dn = Dets[1]->detValues;
      ValueType psiNew=0.0;
      vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
      vector<RealType>::iterator it(C.begin()),last(C.end());
      while(it != last) 
        psiNew += (*(it++))*detValues_up[*(upC++)]*detValues_dn[*(dnC++)];
      curRatio = psiNew/psiCurrent;
      RatioTimer.stop();
      return curRatio;
    } else {
      RatioTimer.start();
      Dets[1]->evaluateDetsForPtclMove(P,iat);
      ValueVector_t& detValues_up = Dets[0]->detValues;
      ValueVector_t& detValues_dn = Dets[1]->new_detValues;
      ValueType psiNew=0.0;
      vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
      vector<RealType>::iterator it(C.begin()),last(C.end());
      while(it != last) 
        psiNew += (*(it++))*detValues_up[*(upC++)]*detValues_dn[*(dnC++)];
      curRatio = psiNew/psiCurrent;
      RatioTimer.stop();
      return curRatio;
    }
  }

  void MultiSlaterDeterminantFast::acceptMove(ParticleSet& P, int iat)
  {
// this should depend on the type of update, ratio / ratioGrad 
// for now is incorrect fot ratio(P,iat,dG,dL) updates 

// update psiCurrent,myG_temp,myL_temp

    psiCurrent *= curRatio;
    curRatio=1.0;
    Dets[DetID[iat]]->acceptMove(P,iat); 
    switch(UpdateMode) 
    {
      case ORB_PBYP_ALL:
      // ratio(P,iat,dG,dL)
        myG = myG_temp;
        myL = myL_temp;
        break;
      default: 
        break;
    }
  }

  void MultiSlaterDeterminantFast::restore(int iat)
  {
    Dets[DetID[iat]]->restore(iat);
    curRatio=1.0;
  }

  void MultiSlaterDeterminantFast::update(ParticleSet& P
      , ParticleSet::ParticleGradient_t& dG, ParticleSet::ParticleLaplacian_t& dL
      , int iat)
  {
    APP_ABORT("IMPLEMENT MultiSlaterDeterminantFast::update");
  }

  OrbitalBase::RealType MultiSlaterDeterminantFast::evaluateLog(ParticleSet& P,BufferType& buf)
  {

    Dets[0]->evaluateLog(P,buf);
    Dets[1]->evaluateLog(P,buf);
    buf.put(psiCurrent);
    buf.put(myL.first_address(), myL.last_address());
    buf.put(FirstAddressOfG,LastAddressOfG);

    return LogValue = evaluateLogAndPhase(psiCurrent,PhaseValue);
  }

  OrbitalBase::RealType MultiSlaterDeterminantFast::registerData(ParticleSet& P, BufferType& buf)
  {
    Dets[0]->registerData(P,buf);
    Dets[1]->registerData(P,buf);

    LogValue = evaluateLog(P,P.G,P.L);

    FirstAddressOfG = &myG[0][0];
    LastAddressOfG = FirstAddressOfG + P.getTotalNum()*DIM;

    buf.add(psiCurrent);
    buf.add(myL.first_address(), myL.last_address());
    buf.add(FirstAddressOfG,LastAddressOfG);

// debug, erase 
//    msd->registerData(P,buf);

    return LogValue;
  }

// FIX FIX FIX
  OrbitalBase::RealType MultiSlaterDeterminantFast::updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch)
  {
    UpdateTimer.start();
    if(fromscratch) {
      Dets[0]->updateBuffer(P,buf,true); 
      Dets[1]->updateBuffer(P,buf,true); 
    } else {
// FIX FIX FIX: right now, I need to allow to recalculate 
//    dets, grads and lapls without recomputing orbitals, 
//    for now i always recalculate
      Dets[0]->updateBuffer(P,buf,true); 
      Dets[1]->updateBuffer(P,buf,true); 
    }

    // can this change over time??? I don't know yet
    ValueVector_t& detValues_up = Dets[0]->detValues;
    ValueVector_t& detValues_dn = Dets[1]->detValues;
    GradMatrix_t& grads_up = Dets[0]->grads;
    GradMatrix_t& grads_dn = Dets[1]->grads;
    ValueMatrix_t& lapls_up = Dets[0]->lapls;
    ValueMatrix_t& lapls_dn = Dets[1]->lapls;
    int N1 = Dets[0]->FirstIndex;
    int N2 = Dets[1]->FirstIndex;
    int NP1 = Dets[0]->NumPtcls;
    int NP2 = Dets[1]->NumPtcls;

    psiCurrent=0.0;
    myG=0.0;
    myL=0.0;
    vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
    vector<RealType>::iterator it(C.begin()),last(C.end());
    while(it != last) {
      psiCurrent += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
      for(int k=0,n=N1; k<NP1; k++,n++) {
        myG(n) += (*it)*grads_up(*upC,k)*detValues_dn[*dnC];
        myL(n) += (*it)*lapls_up(*upC,k)*detValues_dn[*dnC];
      }
      for(int k=0,n=N2; k<NP2; k++,n++) {
        myG(n) += (*it)*grads_dn(*dnC,k)*detValues_up[*upC];
        myL(n) += (*it)*lapls_dn(*dnC,k)*detValues_up[*upC];
      }
      it++;upC++;dnC++;
    }
    ValueType psiinv = 1.0/psiCurrent;
    myG *= psiinv;
    myL *= psiinv;
    P.G += myG;
    for(int i=0; i<P.L.size(); i++)
      P.L(i) += myL[i] - dot(myG[i],myG[i]);

    buf.put(psiCurrent);
    buf.put(myL.first_address(), myL.last_address());
    buf.put(FirstAddressOfG,LastAddressOfG);

    UpdateTimer.stop();
    return LogValue = evaluateLogAndPhase(psiCurrent,PhaseValue);
  }

  void MultiSlaterDeterminantFast::copyFromBuffer(ParticleSet& P, BufferType& buf)
  {
    Dets[0]->copyFromBuffer(P,buf);
    Dets[1]->copyFromBuffer(P,buf);
    buf.get(psiCurrent);
    buf.get(myL.first_address(), myL.last_address());
    buf.get(FirstAddressOfG,LastAddressOfG);
  }


  void MultiSlaterDeterminantFast::checkInVariables(opt_variables_type& active)
  {
    if(Optimizable) 
    {
      if(myVars.size()) 
        active.insertFrom(myVars);
      else  
        Optimizable=false;
    }
  }

  void MultiSlaterDeterminantFast::checkOutVariables(const opt_variables_type& active)
  {
    if(Optimizable) myVars.getIndex(active);
  }

  /** resetParameters with optVariables
   *
   * USE_resetParameters
   */
  void MultiSlaterDeterminantFast::resetParameters(const opt_variables_type& active)
  {  
    if(Optimizable) 
    {
      for(int i=0; i<C.size(); i++) 
      {
        int loc=myVars.where(i);
        if(loc>=0) C[i]=myVars[i]=active[loc];
      }
      //for(int i=0; i<SDets.size(); i++) SDets[i]->resetParameters(active);
    }
  }
  void MultiSlaterDeterminantFast::reportStatus(ostream& os)
  {
  }

//   OrbitalBasePtr MultiSlaterDeterminantFast::makeClone(ParticleSet& tqp) const
//   {
//      APP_ABORT("IMPLEMENT OrbitalBase::makeClone");
//      return 0;
//   }

  void MultiSlaterDeterminantFast::evaluateDerivatives(ParticleSet& P, 
      const opt_variables_type& optvars,
      vector<RealType>& dlogpsi,
      vector<RealType>& dhpsioverpsi)
  {
    bool recalculate(false);
    for (int k=0; k<myVars.size(); ++k)
    {
      int kk=myVars.where(k);
      if (kk<0) continue;
      if (optvars.recompute(kk)) recalculate=true;
    }


// need to modify for CSF later on, right now assume Slater Det basis 
    if (recalculate)
    {

      if(laplSum_up.size() == 0)
         laplSum_up.resize(Dets[0]->detValues.size());
      if(laplSum_dn.size() == 0)
         laplSum_dn.resize(Dets[1]->detValues.size());

    Dets[0]->evaluateForWalkerMove(P);
    Dets[1]->evaluateForWalkerMove(P);

      ValueVector_t& detValues_up = Dets[0]->detValues;
      ValueVector_t& detValues_dn = Dets[1]->detValues;
      GradMatrix_t& grads_up = Dets[0]->grads;
      GradMatrix_t& grads_dn = Dets[1]->grads;
      ValueMatrix_t& lapls_up = Dets[0]->lapls;
      ValueMatrix_t& lapls_dn = Dets[1]->lapls;
      int N1 = Dets[0]->FirstIndex;
      int N2 = Dets[1]->FirstIndex;
      int NP1 = Dets[0]->NumPtcls;
      int NP2 = Dets[1]->NumPtcls;

      int n = P.getTotalNum();
      ValueType psiinv = 1.0/psiCurrent;
      ValueType lapl_sum=0.0;
      //ParticleSet::ParticleGradient_t g(n);
      ValueType gg=0.0, ggP=0.0;
      myG_temp=0.0;
      for(int i=0; i<C.size(); i++){
        int upC = C2node_up[i];
        int dnC = C2node_dn[i];
        ValueType tmp1 = C[i]*detValues_dn[dnC]*psiinv;
        ValueType tmp2 = C[i]*detValues_up[upC]*psiinv;
        //lapl_sum += tmp*(Sum(lapls_up[upC])+Sum(lapls_dn[dnC]));
        ValueType& dummy1 = laplSum_up[upC];
        ValueType& dummy2 = laplSum_dn[dnC];
        dummy1=dummy2=0.0; 
        //for(int k=0; k<nels_up; k++) dummy1 += lapls_up(upC,k); 
        ValueType* ptr = lapls_up[upC]; 
        for(int k=0; k<nels_up; k++,ptr++) dummy1 += *ptr; 
        //for(int k=0; k<nels_dn; k++) dummy2 += lapls_dn(dnC,k); 
        ptr = lapls_dn[dnC]; 
        for(int k=0; k<nels_dn; k++,ptr++) dummy2 += *ptr; 
        lapl_sum += tmp1*dummy1+tmp2*dummy2;  
        for(int k=0,j=N1; k<NP1; k++,j++)  
          myG_temp[j] += tmp1*grads_up(upC,k);
        for(int k=0,j=N2; k<NP2; k++,j++) 
          myG_temp[j] += tmp2*grads_dn(dnC,k);
      }
      //gg=Dot(myG_temp,myG_temp);
      //ggP=Dot(P.G,myG_temp);
      gg=ggP=0.0;
      for(int i=0; i<n; i++)
      {
        //gg += dot(myG_temp[i],myG_temp[i]); 
        //ggP += dot(P.G[i],myG_temp[i]); 
        gg += dot(myG_temp[i],myG_temp[i])-dot(P.G[i],myG_temp[i]);
      } 
 
      for(int i=0; i<C.size(); i++){
        int kk=myVars.where(i);
        if (kk<0) continue;
        //dlogpsi[kk] = cdet;
        int upC = C2node_up[i];
        int dnC = C2node_dn[i];
        ValueType cdet=detValues_up[upC]*detValues_dn[dnC]*psiinv;
        ValueType tmp1=detValues_dn[dnC]*psiinv;
        ValueType tmp2=detValues_up[upC]*psiinv;
        convert(cdet,dlogpsi[kk]);
        ValueType v1=0.0,v2=0.0;
        for(int k=0,j=N1; k<NP1; k++,j++) 
          v1 += (dot(P.G[j],grads_up(upC,k))-dot(myG_temp[j],grads_up(upC,k)) );
        for(int k=0,j=N2; k<NP2; k++,j++) 
          v2 += (dot(P.G[j],grads_dn(dnC,k))-dot(myG_temp[j],grads_dn(dnC,k)));
        ValueType dhpsi =  -0.5*(tmp1*laplSum_up[upC]+tmp2*laplSum_dn[dnC]
                                 -cdet*lapl_sum)  
                           -cdet*gg-(tmp1*v1+tmp2*v2);
//        ValueType dhpsi =  (-0.5*cdet)*
//                       ( laplSum_up[upC]+laplSum_dn[dnC]-lapl_sum
//                          +2.0*(gg-Dot(myG_temp,grads_up[upC])-Dot(myG_temp,grads_dn[dnC])
//                          +Dot(P.G,grads_up[upC])+Dot(P.G,grads_dn[dnC])-ggP));
        convert(dhpsi,dhpsioverpsi[kk]);
      }
    }
  }

  void MultiSlaterDeterminantFast::registerTimers()
  {
    RatioTimer.reset();
    RatioGradTimer.reset();
    RatioAllTimer.reset();
    UpdateTimer.reset();
    EvaluateTimer.reset();
    TimerManager.addTimer (&RatioTimer);
    TimerManager.addTimer (&RatioGradTimer);
    TimerManager.addTimer (&RatioAllTimer);
    TimerManager.addTimer (&UpdateTimer);
    TimerManager.addTimer (&EvaluateTimer);
  }


}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3416 $   $Date: 2008-12-07 11:34:49 -0600 (Sun, 07 Dec 2008) $
 * $Id: MultiSlaterDeterminantFast.cpp 3416 2008-12-07 17:34:49Z jnkim $
 ***************************************************************************/
