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

namespace qmcplusplus
{

MultiSlaterDeterminantFast::MultiSlaterDeterminantFast(ParticleSet& targetPtcl, MultiDiracDeterminantBase* up, MultiDiracDeterminantBase* dn):
  RatioTimer("MultiSlaterDeterminantFast::ratio"),
  RatioGradTimer("MultiSlaterDeterminantFast::ratioGrad"),
  RatioAllTimer("MultiSlaterDeterminantFast::ratio(all)"),
  Ratio1Timer("MultiSlaterDeterminantFast::detEval_ratio"),
  Ratio1GradTimer("MultiSlaterDeterminantFast::detEval_ratioGrad"),
  Ratio1AllTimer("MultiSlaterDeterminantFast::detEval_ratio(all)"),
  UpdateTimer("MultiSlaterDeterminantFast::updateBuffer"),
  EvaluateTimer("MultiSlaterDeterminantFast::evaluate"),
  AccRejTimer("MultiSlaterDeterminantFast::Accept_Reject")
{
  registerTimers();
  //Optimizable=true;
  Optimizable=true;
  OrbitalName="MultiSlaterDeterminantFast";
  usingCSF=false;
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
    for(int j=targetPtcl.first(i); j<targetPtcl.last(i); ++j)
      DetID[j]=i;
  usingBF=false;
  BFTrans=0;
}

OrbitalBasePtr MultiSlaterDeterminantFast::makeClone(ParticleSet& tqp) const
{
  MultiDiracDeterminantBase* up_clone = new MultiDiracDeterminantBase(*Dets[0]);
  MultiDiracDeterminantBase* dn_clone = new MultiDiracDeterminantBase(*Dets[1]);
  MultiSlaterDeterminantFast* clone = new MultiSlaterDeterminantFast(tqp,up_clone,dn_clone);
  if(usingBF)
  {
    BackflowTransformation *tr = BFTrans->makeClone(tqp);
    clone->setBF(tr);
  }
  clone->resetTargetParticleSet(tqp);
  clone->C2node_up=C2node_up;
  clone->C2node_dn=C2node_dn;
  clone->Optimizable=Optimizable;
  clone->C=C;
  clone->myVars=myVars;
  clone->usingCSF=usingCSF;
  clone->usingBF=usingBF;
  if (usingCSF)
  {
    clone->CSFcoeff=CSFcoeff;
    clone->CSFexpansion=CSFexpansion;
    clone->DetsPerCSF=DetsPerCSF;
  }
  return clone;
}

MultiSlaterDeterminantFast::~MultiSlaterDeterminantFast() { }

void MultiSlaterDeterminantFast::resetTargetParticleSet(ParticleSet& P)
{
  if(usingBF)
  {
    BFTrans->resetTargetParticleSet(P);
    for(int i=0; i<Dets.size(); i++)
      Dets[i]->resetTargetParticleSet(BFTrans->QP);
  }
  else
  {
    for(int i=0; i<Dets.size(); i++)
      Dets[i]->resetTargetParticleSet(P);
  }
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
//     log = msd->evaluate(P,G,L);
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
//     log = msd->evaluate(P,G,L);
  log0 = evaluate(P,G0,L0);
  PosType dr;
  dr[0] = 0.1;
  dr[1]=0.05;
  dr[2] = -0.01;
  PosType newpos(P.makeMove(iat,dr));
  app_log() <<"Testing ratio(P,dG,dL). \n";
  G=0;
  G0=0;
  L=0;
  L0=0;
//     log = msd->ratio(P,iat,G,L);
  log0 = ratio(P,iat,G0,L0);
  cout<<"Psi: " <<log <<"   " <<log0 <<"   " <<log/log0 <<endl;
  for(int i=0; i<n; i++)
  {
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
  while(it != last)
  {
    psiCurrent += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
    for(int k=0,n=N1; k<NP1; k++,n++)
    {
      myG(n) += (*it)*grads_up(*upC,k)*detValues_dn[*dnC];
      myL(n) += (*it)*lapls_up(*upC,k)*detValues_dn[*dnC];
    }
    for(int k=0,n=N2; k<NP2; k++,n++)
    {
      myG(n) += (*it)*grads_dn(*dnC,k)*detValues_up[*upC];
      myL(n) += (*it)*lapls_dn(*dnC,k)*detValues_up[*upC];
    }
    it++;
    upC++;
    dnC++;
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

OrbitalBase::RealType MultiSlaterDeterminantFast::evaluateLog(ParticleSet& P,
    ParticleSet::ParticleGradient_t& G,
    ParticleSet::ParticleLaplacian_t& L,
    PooledData<RealType>& buf,
    bool fillBuffer )
{
  if(fillBuffer)
  {
    Dets[0]->evaluateForWalkerMove(P);
    Dets[0]->copyToDerivativeBuffer(P,buf);
    Dets[1]->evaluateForWalkerMove(P);
    Dets[1]->copyToDerivativeBuffer(P,buf);
  }
  else
  {
    Dets[0]->copyFromDerivativeBuffer(P,buf);
    Dets[1]->copyFromDerivativeBuffer(P,buf);
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
  while(it != last)
  {
    psiCurrent += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
    for(int k=0,n=N1; k<NP1; k++,n++)
    {
      myG(n) += (*it)*grads_up(*upC,k)*detValues_dn[*dnC];
      myL(n) += (*it)*lapls_up(*upC,k)*detValues_dn[*dnC];
    }
    for(int k=0,n=N2; k<NP2; k++,n++)
    {
      myG(n) += (*it)*grads_dn(*dnC,k)*detValues_up[*upC];
      myL(n) += (*it)*lapls_dn(*dnC,k)*detValues_up[*upC];
    }
    it++;
    upC++;
    dnC++;
  }
  ValueType psiinv = 1.0/psiCurrent;
  myG *= psiinv;
  myL *= psiinv;
  G += myG;
  for(int i=0; i<L.size(); i++)
    L(i) += myL[i] - dot(myG[i],myG[i]);
  return evaluateLogAndPhase(psiCurrent,PhaseValue);
}


OrbitalBase::GradType MultiSlaterDeterminantFast::evalGrad(ParticleSet& P, int iat)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: evalGrad not implemented. \n");
  }
  GradType grad_iat;
  if(DetID[iat] == 0)
  {
    Dets[0]->evaluateGrads(P,iat);
    ValueVector_t& detValues_up = Dets[0]->detValues;
    ValueVector_t& detValues_dn = Dets[1]->detValues;
    GradMatrix_t& grads_up = Dets[0]->grads;
    int N1 = Dets[0]->FirstIndex;
    ValueType psi=0.0;
    vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
    vector<RealType>::iterator it(C.begin()),last(C.end());
    while(it != last)
    {
      psi += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
      grad_iat += (*it)*grads_up(*upC,iat-N1)*detValues_dn[*dnC];
      it++;
      upC++;
      dnC++;
    }
    grad_iat *= 1.0/psi;
    return grad_iat;
  }
  else
  {
    Dets[1]->evaluateGrads(P,iat);
    ValueType psi=0.0;
    ValueVector_t& detValues_up = Dets[0]->detValues;
    ValueVector_t& detValues_dn = Dets[1]->detValues;
    GradMatrix_t& grads_dn = Dets[1]->grads;
    int N2 = Dets[1]->FirstIndex;
    vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
    vector<RealType>::iterator it(C.begin()),last(C.end());
    while(it != last)
    {
      psi += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
      grad_iat += (*it)*grads_dn(*dnC,iat-N2)*detValues_up[*upC];
      it++;
      upC++;
      dnC++;
    }
    grad_iat *= 1.0/psi;
    return grad_iat;
  }
}

OrbitalBase::ValueType MultiSlaterDeterminantFast::ratioGrad(ParticleSet& P
    , int iat, GradType& grad_iat)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: ratioGrad not implemented. \n");
  }
  UpdateMode=ORB_PBYP_PARTIAL;
  if(DetID[iat] == 0)
  {
    RatioGradTimer.start();
    Ratio1GradTimer.start();
    Dets[0]->evaluateDetsAndGradsForPtclMove(P,iat);
    Ratio1GradTimer.stop();
    ValueVector_t& detValues_up = Dets[0]->new_detValues;
    ValueVector_t& detValues_dn = Dets[1]->detValues;
    GradMatrix_t& grads_up = Dets[0]->new_grads;
    int N1 = Dets[0]->FirstIndex;
    vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
    vector<RealType>::iterator it(C.begin()),last(C.end());
    ValueType psiNew=0.0;
    GradType dummy;
    it=C.begin();
    last=C.end();
    while(it != last)
    {
      psiNew += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
      dummy += (*it)*grads_up(*upC,iat-N1)*detValues_dn[*dnC];
      it++;
      upC++;
      dnC++;
    }
    grad_iat+=dummy/psiNew;
    curRatio = psiNew/psiCurrent;
    RatioGradTimer.stop();
    return curRatio;
  }
  else
  {
    RatioGradTimer.start();
    Ratio1GradTimer.start();
    Dets[1]->evaluateDetsAndGradsForPtclMove(P,iat);
    Ratio1GradTimer.stop();
    ValueVector_t& detValues_up = Dets[0]->detValues;
    ValueVector_t& detValues_dn = Dets[1]->new_detValues;
    GradMatrix_t& grads_dn = Dets[1]->new_grads;
    int N2 = Dets[1]->FirstIndex;
    vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
    vector<RealType>::iterator it(C.begin()),last(C.end());
    ValueType psiNew=0.0;
    GradType dummy;
    while(it != last)
    {
      psiNew += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
      dummy += (*it)*grads_dn(*dnC,iat-N2)*detValues_up[*upC];
      it++;
      upC++;
      dnC++;
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
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: ratio(P,dG,dL) not implemented. \n");
  }
  UpdateMode=ORB_PBYP_ALL;
  if(DetID[iat] == 0)
  {
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
    Ratio1AllTimer.start();
    Dets[0]->evaluateAllForPtclMove(P,iat);
    Ratio1AllTimer.stop();
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
    while(it != last)
    {
      psiNew += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
      for(int k=0,n=N1; k<NP1; k++,n++)
      {
        myG_temp(n) += (*it)*grads_up(*upC,k)*detValues_dn[*dnC];
        myL_temp(n) += (*it)*lapls_up(*upC,k)*detValues_dn[*dnC];
      }
      for(int k=0,n=N2; k<NP2; k++,n++)
      {
        myG_temp(n) += (*it)*grads_dn(*dnC,k)*detValues_up[*upC];
        myL_temp(n) += (*it)*lapls_dn(*dnC,k)*detValues_up[*upC];
      }
      it++;
      upC++;
      dnC++;
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
  else
  {
    RatioAllTimer.start();
    Ratio1AllTimer.start();
    Dets[1]->evaluateAllForPtclMove(P,iat);
    Ratio1AllTimer.stop();
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
    while(it != last)
    {
      psiNew += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
      for(int k=0,n=N1; k<NP1; k++,n++)
      {
        myG_temp(n) += (*it)*grads_up(*upC,k)*detValues_dn[*dnC];
        myL_temp(n) += (*it)*lapls_up(*upC,k)*detValues_dn[*dnC];
      }
      for(int k=0,n=N2; k<NP2; k++,n++)
      {
        myG_temp(n) += (*it)*grads_dn(*dnC,k)*detValues_up[*upC];
        myL_temp(n) += (*it)*lapls_dn(*dnC,k)*detValues_up[*upC];
      }
      it++;
      upC++;
      dnC++;
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
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: ratio not implemented. \n");
  }
  UpdateMode=ORB_PBYP_RATIO;
  if(DetID[iat] == 0)
  {
    RatioTimer.start();
    Ratio1Timer.start();
    Dets[0]->evaluateDetsForPtclMove(P,iat);
    Ratio1Timer.stop();
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
  }
  else
  {
    RatioTimer.start();
    Ratio1Timer.start();
    Dets[1]->evaluateDetsForPtclMove(P,iat);
    Ratio1Timer.stop();
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
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: acceptMove not implemented. \n");
  }
// update psiCurrent,myG_temp,myL_temp
  AccRejTimer.start();
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
  AccRejTimer.stop();
//    Dets[0]->evaluateForWalkerMove(P);
//    Dets[1]->evaluateForWalkerMove(P);
  // can this change over time??? I don't know yet
  /*
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

      ValueType psi=0.0;
      myG_temp=0.0;
      myL_temp=0.0;
      vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
      vector<RealType>::iterator it(C.begin()),last(C.end());
      while(it != last) {
        psi += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
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
      ValueType psiinv = 1.0/psi;
      myG_temp *= psiinv;
      myL_temp *= psiinv;
  */
}

void MultiSlaterDeterminantFast::restore(int iat)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: restore not implemented. \n");
  }
  AccRejTimer.start();
  Dets[DetID[iat]]->restore(iat);
  curRatio=1.0;
  AccRejTimer.stop();
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
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: restore not implemented. \n");
  }
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

// this routine does not initialize the data, just reserves the space
void MultiSlaterDeterminantFast::registerDataForDerivatives(ParticleSet& P, BufferType& buf, int storageType)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: registerDataForDerivatives not implemented. \n");
  }
  Dets[0]->registerDataForDerivatives(P,buf,storageType);
  Dets[1]->registerDataForDerivatives(P,buf,storageType);
}

// FIX FIX FIX
OrbitalBase::RealType MultiSlaterDeterminantFast::updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch)
{
  UpdateTimer.start();
  Dets[0]->updateBuffer(P,buf,fromscratch);
  Dets[1]->updateBuffer(P,buf,fromscratch);
  //Dets[0]->updateBuffer(P,buf,true);
  //Dets[1]->updateBuffer(P,buf,true);
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
  while(it != last)
  {
    psiCurrent += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
    for(int k=0,n=N1; k<NP1; k++,n++)
    {
      myG(n) += (*it)*grads_up(*upC,k)*detValues_dn[*dnC];
      myL(n) += (*it)*lapls_up(*upC,k)*detValues_dn[*dnC];
    }
    for(int k=0,n=N2; k<NP2; k++,n++)
    {
      myG(n) += (*it)*grads_dn(*dnC,k)*detValues_up[*upC];
      myL(n) += (*it)*lapls_dn(*dnC,k)*detValues_up[*upC];
    }
    it++;
    upC++;
    dnC++;
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
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: copyFromBuffer not implemented. \n");
  }
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
  if(Optimizable)
    myVars.getIndex(active);
}

/** resetParameters with optVariables
 *
 * USE_resetParameters
 */
void MultiSlaterDeterminantFast::resetParameters(const opt_variables_type& active)
{
  if(Optimizable)
  {
    if(usingCSF)
    {
      for(int i=0; i<CSFcoeff.size()-1; i++)
      {
        int loc=myVars.where(i);
        if(loc>=0)
          CSFcoeff[i+1]=myVars[i]=active[loc];
      }
      int cnt=0;
      for(int i=0; i<DetsPerCSF.size(); i++)
      {
        for(int k=0; k<DetsPerCSF[i]; k++)
        {
          C[cnt] = CSFcoeff[i]*CSFexpansion[cnt];
          cnt++;
        }
      }
      //for(int i=0; i<Dets.size(); i++) Dets[i]->resetParameters(active);
    }
    else
    {
      for(int i=0; i<C.size()-1; i++)
      {
        int loc=myVars.where(i);
        if(loc>=0)
          C[i+1]=myVars[i]=active[loc];
      }
      //for(int i=0; i<Dets.size(); i++) Dets[i]->resetParameters(active);
    }
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
    if (kk<0)
      continue;
    if (optvars.recompute(kk))
      recalculate=true;
  }
// need to modify for CSF later on, right now assume Slater Det basis
  if (recalculate)
  {
    if(usingCSF)
    {
      if(laplSum_up.size() == 0)
        laplSum_up.resize(Dets[0]->detValues.size());
      if(laplSum_dn.size() == 0)
        laplSum_dn.resize(Dets[1]->detValues.size());
      // assume that evaluateLog has been called in opt routine before
      //   Dets[0]->evaluateForWalkerMove(P);
      //   Dets[1]->evaluateForWalkerMove(P);
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
// myG,myL should already be calculated
      int n = P.getTotalNum();
      ValueType psiinv = 1.0/psiCurrent;
      ValueType lapl_sum=0.0;
      ValueType gg=0.0, ggP=0.0;
      myG_temp=0.0;
      int num=laplSum_up.size();
      ValueVector_t::iterator it(laplSum_up.begin());
      ValueVector_t::iterator last(laplSum_up.end());
      ValueType* ptr0 = lapls_up[0];
      while(it != last)
      {
        (*it)=0.0;
        for(int k=0; k<nels_up; k++,ptr0++)
          (*it) += *ptr0;
        it++;
      }
      it=laplSum_dn.begin();
      last=laplSum_dn.end();
      ptr0 = lapls_dn[0];
      while(it != last)
      {
        (*it)=0.0;
        for(int k=0; k<nels_dn; k++,ptr0++)
          (*it) += *ptr0;
        it++;
      }
      for(int i=0; i<C.size(); i++)
      {
        int upC = C2node_up[i];
        int dnC = C2node_dn[i];
        ValueType tmp1 = C[i]*detValues_dn[dnC]*psiinv;
        ValueType tmp2 = C[i]*detValues_up[upC]*psiinv;
        lapl_sum += tmp1*laplSum_up[upC]+tmp2*laplSum_dn[dnC];
        for(int k=0,j=N1; k<NP1; k++,j++)
          myG_temp[j] += tmp1*grads_up(upC,k);
        for(int k=0,j=N2; k<NP2; k++,j++)
          myG_temp[j] += tmp2*grads_dn(dnC,k);
      }
      gg=ggP=0.0;
      for(int i=0; i<n; i++)
      {
        gg += dot(myG_temp[i],myG_temp[i])-dot(P.G[i],myG_temp[i]);
      }
//       for(int i=0; i<C.size(); i++){
      num=CSFcoeff.size()-1;
      int cnt=0;
//        this one is not optable
      cnt+=DetsPerCSF[0];
      int ip(1);
      for(int i=0; i<num; i++,ip++)
      {
        int kk=myVars.where(i);
        if (kk<0)
        {
          cnt+=DetsPerCSF[ip];
          continue;
        }
        ValueType cdet=0.0,q0=0.0,v1=0.0,v2=0.0;
        for(int k=0; k<DetsPerCSF[ip]; k++)
        {
          int upC = C2node_up[cnt];
          int dnC = C2node_dn[cnt];
          ValueType tmp1=CSFexpansion[cnt]*detValues_dn[dnC]*psiinv;
          ValueType tmp2=CSFexpansion[cnt]*detValues_up[upC]*psiinv;
          cdet+=CSFexpansion[cnt]*detValues_up[upC]*detValues_dn[dnC]*psiinv;
          q0 += (tmp1*laplSum_up[upC] + tmp2*laplSum_dn[dnC]);
          for(int l=0,j=N1; l<NP1; l++,j++)
            v1 += tmp1*(dot(P.G[j],grads_up(upC,l))-dot(myG_temp[j],grads_up(upC,l)) );
          for(int l=0,j=N2; l<NP2; l++,j++)
            v2 += tmp2*(dot(P.G[j],grads_dn(dnC,l))-dot(myG_temp[j],grads_dn(dnC,l)));
          cnt++;
        }
        convert(cdet,dlogpsi[kk]);
        ValueType dhpsi =  -0.5*(q0-cdet*lapl_sum)
                           -cdet*gg-v1-v2;
        //ValueType dhpsi =  -0.5*(tmp1*laplSum_up[upC]+tmp2*laplSum_dn[dnC]
        //                         -cdet*lapl_sum)
        //                   -cdet*gg-(tmp1*v1+tmp2*v2);
        convert(dhpsi,dhpsioverpsi[kk]);
      }
    }
    else
      //usingCSF
    {
      if(laplSum_up.size() == 0)
        laplSum_up.resize(Dets[0]->detValues.size());
      if(laplSum_dn.size() == 0)
        laplSum_dn.resize(Dets[1]->detValues.size());
      // assume that evaluateLog has been called in opt routine before
      //   Dets[0]->evaluateForWalkerMove(P);
      //   Dets[1]->evaluateForWalkerMove(P);
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
      ValueType gg=0.0, ggP=0.0;
      myG_temp=0.0;
      int num=laplSum_up.size();
      ValueVector_t::iterator it(laplSum_up.begin());
      ValueVector_t::iterator last(laplSum_up.end());
      ValueType* ptr0 = lapls_up[0];
      while(it != last)
      {
        (*it)=0.0;
        for(int k=0; k<nels_up; k++,ptr0++)
          (*it) += *ptr0;
        it++;
      }
      it=laplSum_dn.begin();
      last=laplSum_dn.end();
      ptr0 = lapls_dn[0];
      while(it != last)
      {
        (*it)=0.0;
        for(int k=0; k<nels_dn; k++,ptr0++)
          (*it) += *ptr0;
        it++;
      }
      for(int i=0; i<C.size(); i++)
      {
        int upC = C2node_up[i];
        int dnC = C2node_dn[i];
        ValueType tmp1 = C[i]*detValues_dn[dnC]*psiinv;
        ValueType tmp2 = C[i]*detValues_up[upC]*psiinv;
        lapl_sum += tmp1*laplSum_up[upC]+tmp2*laplSum_dn[dnC];
        for(int k=0,j=N1; k<NP1; k++,j++)
          myG_temp[j] += tmp1*grads_up(upC,k);
        for(int k=0,j=N2; k<NP2; k++,j++)
          myG_temp[j] += tmp2*grads_dn(dnC,k);
      }
      gg=ggP=0.0;
      for(int i=0; i<n; i++)
      {
        gg += dot(myG_temp[i],myG_temp[i])-dot(P.G[i],myG_temp[i]);
      }
      for(int i=1; i<C.size(); i++)
      {
        int kk=myVars.where(i-1);
        if (kk<0)
          continue;
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
        convert(dhpsi,dhpsioverpsi[kk]);
      }
    } // usingCSF
  }
}

void MultiSlaterDeterminantFast::registerTimers()
{
  RatioTimer.reset();
  RatioGradTimer.reset();
  RatioAllTimer.reset();
  Ratio1Timer.reset();
  Ratio1GradTimer.reset();
  Ratio1AllTimer.reset();
  UpdateTimer.reset();
  EvaluateTimer.reset();
  AccRejTimer.reset();
  TimerManager.addTimer (&RatioTimer);
  TimerManager.addTimer (&RatioGradTimer);
  TimerManager.addTimer (&RatioAllTimer);
  TimerManager.addTimer (&Ratio1Timer);
  TimerManager.addTimer (&Ratio1GradTimer);
  TimerManager.addTimer (&Ratio1AllTimer);
  TimerManager.addTimer (&UpdateTimer);
  TimerManager.addTimer (&EvaluateTimer);
  TimerManager.addTimer (&AccRejTimer);
}


}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3416 $   $Date: 2008-12-07 11:34:49 -0600 (Sun, 07 Dec 2008) $
 * $Id: MultiSlaterDeterminantFast.cpp 3416 2008-12-07 17:34:49Z jnkim $
 ***************************************************************************/
