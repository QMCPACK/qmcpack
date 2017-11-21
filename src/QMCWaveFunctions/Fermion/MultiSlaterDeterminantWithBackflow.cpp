//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCWaveFunctions/Fermion/MultiSlaterDeterminantWithBackflow.h"
#include "ParticleBase/ParticleAttribOps.h"

namespace qmcplusplus
{

MultiSlaterDeterminantWithBackflow::MultiSlaterDeterminantWithBackflow(ParticleSet& targetPtcl, SPOSetProxyPtr upspo, SPOSetProxyPtr dnspo, BackflowTransformation *BF):MultiSlaterDeterminant(targetPtcl,upspo,dnspo),BFTrans(BF)
{
  Optimizable=false;
  OrbitalName="MultiSlaterDeterminantWithBackflow";
}

OrbitalBasePtr MultiSlaterDeterminantWithBackflow::makeClone(ParticleSet& tqp) const
{
  // mmorales: the proxy classes read from the particle set inside BFTrans
  BackflowTransformation *tr = BFTrans->makeClone(tqp);
  tr->resetTargetParticleSet(tqp);
  SPOSetProxyForMSD* spo_up_C = new SPOSetProxyForMSD(spo_up->refPhi->makeClone(),FirstIndex_up,LastIndex_up);
  SPOSetProxyForMSD* spo_dn_C = new SPOSetProxyForMSD(spo_dn->refPhi->makeClone(),FirstIndex_dn,LastIndex_dn);
  spo_up_C->occup= spo_up->occup;
  spo_dn_C->occup= spo_dn->occup;
  spo_up_C->refPhi->resetTargetParticleSet(tr->QP);
  spo_dn_C->refPhi->resetTargetParticleSet(tr->QP);
  MultiSlaterDeterminantWithBackflow* clone = new MultiSlaterDeterminantWithBackflow(tqp,spo_up_C,spo_dn_C,tr);
  clone->C2node_up=C2node_up;
  clone->C2node_dn=C2node_dn;
  clone->resize(dets_up.size(),dets_dn.size());
  if (usingCSF)
  {
    clone->CSFcoeff=CSFcoeff;
    clone->CSFexpansion=CSFexpansion;
    clone->DetsPerCSF=DetsPerCSF;
  }
  for(int i=0; i<dets_up.size(); i++)
  {
    DiracDeterminantWithBackflow* dclne = (DiracDeterminantWithBackflow*) dets_up[i]->makeCopy((SPOSetBasePtr) clone->spo_up);
    dclne->BFTrans=tr;
    dclne->resetTargetParticleSet(tr->QP);
    clone->dets_up.push_back(dclne);
  }
  for(int i=0; i<dets_dn.size(); i++)
  {
    DiracDeterminantWithBackflow* dclne = (DiracDeterminantWithBackflow*) dets_dn[i]->makeCopy((SPOSetBasePtr) clone->spo_dn);
    dclne->BFTrans=tr;
    dclne->resetTargetParticleSet(tr->QP);
    clone->dets_dn.push_back(dclne);
  }
  clone->Optimizable=Optimizable;
  clone->C=C;
  clone->myVars=myVars;
  clone->resetTargetParticleSet(tr->QP);
  return clone;
}


MultiSlaterDeterminantWithBackflow::~MultiSlaterDeterminantWithBackflow() { }

void MultiSlaterDeterminantWithBackflow::resize(int n1, int n2)
{
  NumUniqueDets_up=n1;
  NumUniqueDets_dn=n2;
  myG.resize(nels_up+nels_dn);
  myL.resize(nels_up+nels_dn);
  grads_up.resize(NumUniqueDets_up);
  grads_dn.resize(NumUniqueDets_dn);
  lapls_up.resize(NumUniqueDets_up);
  lapls_dn.resize(NumUniqueDets_dn);
  for(int i=0; i<NumUniqueDets_up; i++)
  {
    grads_up[i].resize(nels_up+nels_dn);
    lapls_up[i].resize(nels_up+nels_dn);
  }
  for(int i=0; i<NumUniqueDets_dn; i++)
  {
    grads_dn[i].resize(nels_up+nels_dn);
    lapls_dn[i].resize(nels_up+nels_dn);
  }
  detValues_up.resize(NumUniqueDets_up);
  detValues_dn.resize(NumUniqueDets_dn);
  int maxDet = NumUniqueDets_up;
  if(NumUniqueDets_dn > NumUniqueDets_up)
    maxDet=NumUniqueDets_dn;
  detsRatios.resize(maxDet);
  tempstorage_up.resize(NumUniqueDets_up);
  tempstorage_dn.resize(NumUniqueDets_dn);
  grad_temp.resize(maxDet);
  tempgrad.resize(maxDet);
  templapl.resize(maxDet);
  for(int i=0; i<maxDet; i++)
  {
    tempgrad[i].resize(nels_up+nels_dn);
    templapl[i].resize(nels_up+nels_dn);
  }
}

OrbitalBase::ValueType MultiSlaterDeterminantWithBackflow::evaluate(ParticleSet& P
    , ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
{
  EvaluateTimer.start();
  BFTrans->evaluate(P);
// mmorales: For now always assume 2 electron components, up/down
// trouble: if this is called during an optimization routine, you will need to
//          evaluate grad_grad_grad structures for BF derivatives
//  FIX FIX FIX
  spo_up->evaluateForWalkerMoveWithHessian(BFTrans->QP,FirstIndex_up,LastIndex_up);
  spo_dn->evaluateForWalkerMoveWithHessian(BFTrans->QP,FirstIndex_dn,LastIndex_dn);
  int numP = grads_up[0].size();
  for(int i=0; i<dets_up.size(); i++)
  {
    spo_up->prepareFor(i);
    grads_up[i]=0.0;
    lapls_up[i]=0.0;
    detValues_up[i]=dets_up[i]->evaluate(BFTrans->QP,grads_up[i],lapls_up[i]);
    // need \nabla^2 Det / Det
    for(int k=0; k<numP; k++)
      lapls_up[i][k] += dot(grads_up[i][k],grads_up[i][k]);
  }
  for(int i=0; i<dets_dn.size(); i++)
  {
    spo_dn->prepareFor(i);
    grads_dn[i]=0.0;
    lapls_dn[i]=0.0;
    detValues_dn[i]=dets_dn[i]->evaluate(BFTrans->QP,grads_dn[i],lapls_dn[i]);
    // need \nabla^2 Det / Det
    for(int k=0; k<numP; k++)
      lapls_dn[i][k] += dot(grads_dn[i][k],grads_dn[i][k]);
  }
  ValueType psi=0.0;
  myG=0.0;
  myL=0.0;
  for(int i=0; i<C.size(); i++)
  {
    int upC = C2node_up[i];
    int dnC = C2node_dn[i];
    ParticleSet::ParticleValue_t tmp = C[i]*detValues_up[upC]*detValues_dn[dnC];
    psi += tmp;
    myG += grads_up[upC]*tmp;
    myG += grads_dn[dnC]*tmp;
    myL += lapls_up[upC]*tmp;
    myL += lapls_dn[dnC]*tmp;
    for(int k=0; k<numP; k++)
      myL[k] += 2.0*static_cast<ParticleSet::ParticleValue_t>(tmp)*dot(grads_up[upC][k],grads_dn[dnC][k]);
  }
  ValueType psiinv = (RealType)1.0/psi;
  myG *= psiinv;
  myL *= psiinv;
  G += myG;
  for(int i=0; i<L.size(); i++)
    L[i] += myL[i] - dot(myG[i],myG[i]);
  EvaluateTimer.stop();
  return psi;
}

OrbitalBase::RealType MultiSlaterDeterminantWithBackflow::evaluateLog(ParticleSet& P
    , ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
{
  ValueType psi = evaluate(P,G,L);
  return LogValue = evaluateLogAndPhase(psi,PhaseValue);
}

OrbitalBase::GradType MultiSlaterDeterminantWithBackflow::evalGrad(ParticleSet& P, int iat)
{
  APP_ABORT("MultiSlaterDeterminantWithBackflow:: pbyp routines not implemented ");
  GradType grad_iat;
  if(DetID[iat] == 0)
  {
    for(int i=0; i<dets_up.size(); i++)
    {
      spo_up->prepareFor(i);
      grads_up[i][iat] = dets_up[i]->evalGrad(P,iat);
    }
    ValueType psi=0.0;
    for(int i=0; i<C.size(); i++)
    {
      int upC = C2node_up[i];
      int dnC = C2node_dn[i];
      ParticleSet::ParticleValue_t tmp = C[i]*detValues_up[upC]*detValues_dn[dnC];
      psi += tmp;
      grad_iat += grads_up[upC][iat]*tmp;
    }
    grad_iat *= (RealType)1.0/psi;
    return grad_iat;
  }
  else
  {
    ValueType psi=0.0;
    for(int i=0; i<dets_dn.size(); i++)
    {
      spo_dn->prepareFor(i);
      grads_dn[i][iat] = dets_dn[i]->evalGrad(P,iat);
    }
    for(int i=0; i<C.size(); i++)
    {
      int upC = C2node_up[i];
      int dnC = C2node_dn[i];
      ParticleSet::ParticleValue_t tmp = C[i]*detValues_up[upC]*detValues_dn[dnC];
      psi += tmp;
      grad_iat += grads_dn[dnC][iat]*tmp;
    }
    grad_iat *= (RealType)1.0/psi;
    return grad_iat;
  }
}

OrbitalBase::ValueType MultiSlaterDeterminantWithBackflow::ratioGrad(ParticleSet& P
    , int iat, GradType& grad_iat)
{
  APP_ABORT("MultiSlaterDeterminantWithBackflow:: pbyp routines not implemented ");
  UpdateMode=ORB_PBYP_PARTIAL;
  if(DetID[iat] == 0)
  {
    RatioGradTimer.start();
    evalOrbTimer.start();
    spo_up->evaluateAllForPtclMove(P,iat);
    evalOrbTimer.stop();
    Ratio1GradTimer.start();
    grad_temp=0.0;
    for(int i=0; i<dets_up.size(); i++)
    {
      spo_up->prepareFor(i);
      detsRatios[i]=dets_up[i]->ratioGrad(P,iat,grad_temp[i]);
    }
    Ratio1GradTimer.stop();
    ValueType psiOld=0.0,psiNew=0.0;
    GradType dummy;
    for(int i=0; i<C.size(); i++)
    {
      int upC = C2node_up[i];
      int dnC = C2node_dn[i];
      ValueType tmp2 = C[i]*detValues_up[upC]*detValues_dn[dnC];
      ValueType tmp = tmp2*detsRatios[upC];
      psiOld += tmp2;
      psiNew += tmp;
      dummy += tmp*grad_temp[upC];
    }
    grad_iat+=dummy/psiNew;
    curRatio = psiNew/psiOld;
    RatioGradTimer.stop();
    return curRatio;
  }
  else
  {
    RatioGradTimer.start();
    evalOrbTimer.start();
    spo_dn->evaluateAllForPtclMove(P,iat);
    evalOrbTimer.stop();
    Ratio1GradTimer.start();
    grad_temp=0.0;
    for(int i=0; i<dets_dn.size(); i++)
    {
      spo_dn->prepareFor(i);
      detsRatios[i]=dets_dn[i]->ratioGrad(P,iat,grad_temp[i]);
    }
    Ratio1GradTimer.stop();
    ValueType psiOld=0.0,psiNew=0.0;
    GradType dummy;
    for(int i=0; i<C.size(); i++)
    {
      int upC = C2node_up[i];
      int dnC = C2node_dn[i];
      ValueType tmp2 = C[i]*detValues_up[upC]*detValues_dn[dnC];
      ValueType tmp = tmp2*detsRatios[dnC];
      psiOld += tmp2;
      psiNew += tmp;
      dummy += tmp*grad_temp[dnC];
    }
    grad_iat+=dummy/psiNew;
    curRatio = psiNew/psiOld;
    RatioGradTimer.stop();
    return curRatio;
  }
}

// use ci_node for this routine only
OrbitalBase::ValueType MultiSlaterDeterminantWithBackflow::ratio(ParticleSet& P, int iat)
{
  APP_ABORT("MultiSlaterDeterminantWithBackflow:: pbyp routines not implemented ");
  UpdateMode=ORB_PBYP_RATIO;
  if(DetID[iat] == 0)
  {
    RatioTimer.start();
    spo_up->evaluateForPtclMove(P,iat);
    Ratio1Timer.start();
    for(int i=0; i<dets_up.size(); i++)
    {
      spo_up->prepareFor(i);
      detsRatios[i]=dets_up[i]->ratio(P,iat);
    }
    Ratio1Timer.stop();
    std::vector<size_t>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
    std::vector<RealType>::iterator it(C.begin()),last(C.end());
    ValueType psiOld=0.0,psiNew=0.0;
    while(it != last)
    {
      ValueType tmp = (*it)*detValues_up[*upC]*detValues_dn[*dnC];
      psiNew += tmp*detsRatios[*upC];
      psiOld += tmp;
      it++;
      upC++;
      dnC++;
    }
    curRatio = psiNew/psiOld;
    RatioTimer.stop();
    return curRatio;
  }
  else
  {
    RatioTimer.start();
    spo_dn->evaluateForPtclMove(P,iat);
    Ratio1Timer.start();
    for(int i=0; i<dets_dn.size(); i++)
    {
      spo_dn->prepareFor(i);
      detsRatios[i]=dets_dn[i]->ratio(P,iat);
    }
    Ratio1Timer.stop();
    std::vector<size_t>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
    std::vector<RealType>::iterator it(C.begin()),last(C.end());
    ValueType psiOld=0.0,psiNew=0.0;
    while(it != last)
    {
      ValueType tmp = (*it)*detValues_up[*upC]*detValues_dn[*dnC];
      psiNew += tmp*detsRatios[*dnC];
      psiOld += tmp;
      it++;
      upC++;
      dnC++;
    }
    curRatio = psiNew/psiOld;
    RatioTimer.stop();
    return curRatio;
  }
}

void MultiSlaterDeterminantWithBackflow::acceptMove(ParticleSet& P, int iat)
{
// this should depend on the type of update, ratio / ratioGrad
// for now is incorrect fot ratio(P,iat,dG,dL) updates
  APP_ABORT("MultiSlaterDeterminantWithBackflow:: pbyp routines not implemented ");
  //BFTrans->acceptMove(P,iat);
  AccRejTimer.start();
  if(DetID[iat] == 0)
  {
    for(int i=0; i<dets_up.size(); i++)
      dets_up[i]->acceptMove(P,iat);
    switch(UpdateMode)
    {
    case ORB_PBYP_RATIO:
      // ratio(P,iat)
      for(int i=0; i<detValues_up.size(); i++)
        detValues_up[i] *= detsRatios[i];
      PhaseValue += evaluatePhase(curRatio);
      LogValue +=std::log(std::abs(curRatio));
      curRatio=1.0;
      break;
    case ORB_PBYP_PARTIAL:
      // ratioGrad(P,iat,grad)
      for(int i=0; i<detValues_up.size(); i++)
      {
        detValues_up[i] *= detsRatios[i];
        grads_up[i][iat] = grad_temp[i];
      }
      PhaseValue += evaluatePhase(curRatio);
      LogValue +=std::log(std::abs(curRatio));
      curRatio=1.0;
      break;
    case ORB_PBYP_ALL:
      // ratio(P,iat,dG,dL)
      for(int i=0; i<detValues_up.size(); i++)
      {
        detValues_up[i] *= detsRatios[i];
        grads_up[i] = tempgrad[i];
        lapls_up[i] = templapl[i];
      }
      PhaseValue += evaluatePhase(curRatio);
      LogValue +=std::log(std::abs(curRatio));
      curRatio=1.0;
      break;
    default:
      for(int i=0; i<detValues_up.size(); i++)
        detValues_up[i] *= detsRatios[i];
      PhaseValue += evaluatePhase(curRatio);
      LogValue +=std::log(std::abs(curRatio));
      curRatio=1.0;
      break;
    }
  }
  else
  {
    for(int i=0; i<dets_dn.size(); i++)
      dets_dn[i]->acceptMove(P,iat);
    switch(UpdateMode)
    {
    case ORB_PBYP_RATIO:
      // ratio(P,iat)
      for(int i=0; i<detValues_dn.size(); i++)
        detValues_dn[i] *= detsRatios[i];
      PhaseValue += evaluatePhase(curRatio);
      LogValue +=std::log(std::abs(curRatio));
      curRatio=1.0;
      break;
    case ORB_PBYP_PARTIAL:
      // ratioGrad(P,iat,grad)
      for(int i=0; i<detValues_dn.size(); i++)
      {
        detValues_dn[i] *= detsRatios[i];
        grads_dn[i][iat] = grad_temp[i];
      }
      PhaseValue += evaluatePhase(curRatio);
      LogValue +=std::log(std::abs(curRatio));
      curRatio=1.0;
      break;
    case ORB_PBYP_ALL:
      // ratio(P,iat,dG,dL)
      for(int i=0; i<detValues_dn.size(); i++)
      {
        detValues_dn[i] *= detsRatios[i];
        grads_dn[i] = tempgrad[i];
        lapls_dn[i] = templapl[i];
      }
      PhaseValue += evaluatePhase(curRatio);
      LogValue +=std::log(std::abs(curRatio));
      curRatio=1.0;
      break;
    default:
      for(int i=0; i<detValues_dn.size(); i++)
        detValues_dn[i] *= detsRatios[i];
      PhaseValue += evaluatePhase(curRatio);
      LogValue +=std::log(std::abs(curRatio));
      curRatio=1.0;
      break;
    }
  }
  AccRejTimer.stop();
}

void MultiSlaterDeterminantWithBackflow::restore(int iat)
{
  APP_ABORT("MultiSlaterDeterminantWithBackflow:: pbyp routines not implemented ");
  AccRejTimer.start();
  if(DetID[iat] == 0)
  {
    for(int i=0; i<dets_up.size(); i++)
      dets_up[i]->restore(iat);
  }
  else
  {
    for(int i=0; i<dets_dn.size(); i++)
      dets_dn[i]->restore(iat);
  }
  curRatio=1.0;
  AccRejTimer.stop();
}

void MultiSlaterDeterminantWithBackflow::registerData(ParticleSet& P, WFBufferType& buf)
{
  BFTrans->evaluate(P);
// move resize of pbyp structures to here
  spo_up->evaluateForWalkerMoveWithHessian(BFTrans->QP,FirstIndex_up,LastIndex_up);
  spo_dn->evaluateForWalkerMoveWithHessian(BFTrans->QP,FirstIndex_dn,LastIndex_dn);
  myG = P.G;
  myL = P.L;
  for (int i=0; i<dets_up.size(); i++)
  {
    spo_up->prepareFor(i);
    dets_up[i]->registerData(BFTrans->QP,buf);
  }
  for (int i=0; i<dets_dn.size(); i++)
  {
    spo_dn->prepareFor(i);
    dets_dn[i]->registerData(BFTrans->QP,buf);
  }
  P.G = myG;
  P.L = myL;
  ValueType logpsi = evaluateLog(P,P.G,P.L);
  int TotalDim = PosType::Size*P.getTotalNum();
  //buf.add(detValues_up.begin(),detValues_up.end());
  //buf.add(detValues_dn.begin(),detValues_dn.end());
  buf.add(detValues_up.first_address(),detValues_up.last_address());
  buf.add(detValues_dn.first_address(),detValues_dn.last_address());
  for(int i=0; i<NumUniqueDets_up; i++)
  {
    buf.add(&(grads_up[i][0][0]), &(grads_up[i][0][0])+TotalDim);
    buf.add(lapls_up[i].first_address(),lapls_up[i].last_address());
  }
  for(int i=0; i<NumUniqueDets_dn; i++)
  {
    buf.add(&(grads_dn[i][0][0]), &(grads_dn[i][0][0])+TotalDim);
    buf.add(lapls_dn[i].first_address(),lapls_dn[i].last_address());
  }
}

// FIX FIX FIX
OrbitalBase::RealType MultiSlaterDeterminantWithBackflow::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch)
{
  UpdateTimer.start();
  if(fromscratch || UpdateMode == ORB_PBYP_RATIO)
  {
    BFTrans->evaluate(P);
    spo_up->evaluateForWalkerMoveWithHessian(BFTrans->QP,FirstIndex_up,LastIndex_up);
    spo_dn->evaluateForWalkerMoveWithHessian(BFTrans->QP,FirstIndex_dn,LastIndex_dn);
  }
  myG = P.G;
  myL = P.L;
  RealType logpsi(0.0);
  PhaseValue=0.0;
  for (int i=0; i<dets_up.size(); i++)
  {
    spo_up->prepareFor(i);
    logpsi = dets_up[i]->updateBuffer(BFTrans->QP,buf,fromscratch);
#if defined(QMC_COMPLEX)
    RealType ratioMag = std::exp(logpsi);
    detValues_up[i]= std::complex<OHMMS_PRECISION>(std::cos(dets_up[i]->PhaseValue)*ratioMag,std::sin(dets_up[i]->PhaseValue)*ratioMag);
#else
    detValues_up[i]=std::cos(dets_up[i]->PhaseValue)*std::exp(logpsi);
#endif
    grads_up[i]=dets_up[i]->myG;
    lapls_up[i]=dets_up[i]->myL;
    for(int k=FirstIndex_up; k<LastIndex_up; k++)
      lapls_up[i][k] += dot(grads_up[i][k],grads_up[i][k]);
  }
  for (int i=0; i<dets_dn.size(); i++)
  {
    spo_dn->prepareFor(i);
    logpsi = dets_dn[i]->updateBuffer(BFTrans->QP,buf,fromscratch);
#if defined(QMC_COMPLEX)
    RealType ratioMag = std::exp(logpsi);
    detValues_dn[i]= std::complex<OHMMS_PRECISION>(std::cos(dets_dn[i]->PhaseValue)*ratioMag,std::sin(dets_dn[i]->PhaseValue)*ratioMag);
#else
    detValues_dn[i]=std::cos(dets_dn[i]->PhaseValue)*std::exp(logpsi);
#endif
    grads_dn[i]=dets_dn[i]->myG;
    lapls_dn[i]=dets_dn[i]->myL;
    for(int k=FirstIndex_dn; k<LastIndex_dn; k++)
      lapls_dn[i][k] += dot(grads_dn[i][k],grads_dn[i][k]);
  }
  int TotalDim = PosType::Size*P.getTotalNum();
  buf.put(detValues_up.first_address(),detValues_up.last_address());
  buf.put(detValues_dn.first_address(),detValues_dn.last_address());
  for(int i=0; i<NumUniqueDets_up; i++)
  {
    buf.put(&(grads_up[i][0][0]), &(grads_up[i][0][0])+TotalDim);
    buf.put(lapls_up[i].first_address(),lapls_up[i].last_address());
  }
  for(int i=0; i<NumUniqueDets_dn; i++)
  {
    buf.put(&(grads_dn[i][0][0]), &(grads_dn[i][0][0])+TotalDim);
    buf.put(lapls_dn[i].first_address(),lapls_dn[i].last_address());
  }
  P.G=myG;
  P.L=myL;
  ValueType psi=0.0;
  myG=0.0;
  myL=0.0;
  for(int i=0; i<C.size(); i++)
  {
    int upC = C2node_up[i];
    int dnC = C2node_dn[i];
    ParticleSet::ParticleValue_t tmp = C[i]*detValues_up[upC]*detValues_dn[dnC];
    psi += tmp;
    myG += grads_up[upC]*tmp; // other spin sector should be zero
    myG += grads_dn[dnC]*tmp; // other spin sector should be zero
    myL += lapls_up[upC]*tmp;
    myL += lapls_dn[dnC]*tmp;
  }
  ValueType psiinv = (RealType)1.0/psi;
  myG *= psiinv;
  myL *= psiinv;
  P.G += myG;
  for(int i=0; i<P.L.size(); i++)
    P.L[i] += myL[i] - dot(myG[i],myG[i]);
  UpdateTimer.stop();
  return LogValue = evaluateLogAndPhase(psi,PhaseValue);;
}

void MultiSlaterDeterminantWithBackflow::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  BFTrans->evaluate(P);
  for (int i=0; i<dets_up.size(); i++)
    dets_up[i]->copyFromBuffer(BFTrans->QP,buf);
  for (int i=0; i<dets_dn.size(); i++)
    dets_dn[i]->copyFromBuffer(BFTrans->QP,buf);
  int TotalDim = PosType::Size*P.getTotalNum();
  buf.get(detValues_up.begin(),detValues_up.end());
  buf.get(detValues_dn.begin(),detValues_dn.end());
  for(int i=0; i<NumUniqueDets_up; i++)
  {
    buf.get(&(grads_up[i][0][0]), &(grads_up[i][0][0])+TotalDim);
    buf.get(lapls_up[i].first_address(),lapls_up[i].last_address());
  }
  for(int i=0; i<NumUniqueDets_dn; i++)
  {
    buf.get(&(grads_dn[i][0][0]), &(grads_dn[i][0][0])+TotalDim);
    buf.get(lapls_dn[i].first_address(),lapls_dn[i].last_address());
  }
}


/** resetParameters with optVariables
 *
 * USE_resetParameters
 */
void MultiSlaterDeterminantWithBackflow::resetParameters(const opt_variables_type& active)
{
  if(Optimizable)
  {
    BFTrans->resetParameters(active);
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
    }
    else
    {
      for(int i=0; i<C.size()-1; i++)
      {
        int loc=myVars.where(i);
        if(loc>=0)
          C[i+1]=myVars[i]=active[loc];
      }
    }
  }
}

void MultiSlaterDeterminantWithBackflow::reportStatus(std::ostream& os)
{
  if(Optimizable)
  {
    BFTrans->reportStatus(os);
    myVars.print(os);
  }
}

void MultiSlaterDeterminantWithBackflow::checkInVariables(opt_variables_type& active)
{
  if(Optimizable)
  {
    BFTrans->checkInVariables(active);
    if(myVars.size())
      active.insertFrom(myVars);
  }
}

void MultiSlaterDeterminantWithBackflow::checkOutVariables(const opt_variables_type& active)
{
  if(Optimizable)
  {
    BFTrans->checkOutVariables(active);
    if(myVars.size())
      myVars.getIndex(active);
  }
}

//   OrbitalBasePtr MultiSlaterDeterminantWithBackflow::makeClone(ParticleSet& tqp) const
//   {
//      APP_ABORT("IMPLEMENT OrbitalBase::makeClone");
//      return 0;
//   }

void MultiSlaterDeterminantWithBackflow::evaluateDerivatives(ParticleSet& P,
    const opt_variables_type& optvars,
    std::vector<RealType>& dlogpsi,
    std::vector<RealType>& dhpsioverpsi)
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
  bool optmBF=BFTrans->isOptimizable();
// need to modify for CSF later on, right now assume Slater Det basis
  if (recalculate || optmBF)
  {
    if(usingCSF)
    {
      int n = P.getTotalNum();
#if defined(QMC_COMPLEX)
      RealType ratioMag = std::exp(LogValue);
      ValueType psi = std::complex<OHMMS_PRECISION>(std::cos(PhaseValue)*ratioMag,std::sin(PhaseValue)*ratioMag);
#else
      ValueType psi = std::cos(PhaseValue)*std::exp(LogValue);
#endif
      ValueType psiinv = (RealType)1.0/psi;;
      ValueType lapl_sum=0.0;
      ParticleSet::ParticleGradient_t g(n),gmP(n);
      ValueType gg=0.0, ggP=0.0;
      g=0.0;
      gmP=0.0;
      for(int i=0; i<lapls_up.size(); i++)
        tempstorage_up[i] = Sum(lapls_up[i]);
      for(int i=0; i<lapls_dn.size(); i++)
        tempstorage_dn[i] = Sum(lapls_dn[i]);
      for(int i=0; i<C.size(); i++)
      {
        int upC = C2node_up[i];
        int dnC = C2node_dn[i];
        ValueType tmp = C[i]*detValues_up[upC]*detValues_dn[dnC]*psiinv;
        lapl_sum += tmp*(tempstorage_up[upC]+tempstorage_dn[dnC]
                         +static_cast<ValueType>(2.0*Dot(grads_up[upC],grads_dn[dnC])));
        g += grads_up[upC]*static_cast<ParticleSet::ParticleValue_t>(tmp);
        g += grads_dn[dnC]*static_cast<ParticleSet::ParticleValue_t>(tmp);
      }
      gmP=g-P.G;
      gg=Dot(gmP,g);
      int num=CSFcoeff.size()-1;
      int cnt=0;
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
          ValueType tmp=CSFexpansion[cnt]*detValues_up[upC]*detValues_dn[dnC]*psiinv;
          cdet+=tmp;
          q0 += tmp*(tempstorage_up[upC]+tempstorage_dn[dnC]+static_cast<ValueType>(2.0*Dot(grads_up[upC],grads_dn[dnC])));
          v1 += tmp*static_cast<ValueType>(Dot(gmP,grads_up[upC])+Dot(gmP,grads_dn[dnC]));
          cnt++;
        }
        convert(cdet,dlogpsi[kk]);
        ValueType dhpsi =  (RealType)-0.5*(q0-cdet*lapl_sum)
                           -cdet*gg+v1;
        convert(dhpsi,dhpsioverpsi[kk]);
      }
      if(optmBF)
      {
        // build QP,Amat,Bmat_full,Xmat,Cmat,Ymat
        BFTrans->evaluateDerivatives(P);
      }
    }
    else
    {
      int n = P.getTotalNum();
#if defined(QMC_COMPLEX)
      RealType ratioMag = std::exp(LogValue);
      ValueType psi = std::complex<OHMMS_PRECISION>(std::cos(PhaseValue)*ratioMag,std::sin(PhaseValue)*ratioMag);
#else
      ValueType psi = std::cos(PhaseValue)*std::exp(LogValue);
#endif
      ValueType psiinv = (RealType)1.0/psi;;
      ValueType lapl_sum=0.0;
      ParticleSet::ParticleGradient_t g(n),gmP(n);
      ValueType gg=0.0, ggP=0.0;
      g=0.0;
      gmP=0.0;
      for(int i=0; i<lapls_up.size(); i++)
        tempstorage_up[i] = Sum(lapls_up[i]);
      for(int i=0; i<lapls_dn.size(); i++)
        tempstorage_dn[i] = Sum(lapls_dn[i]);
      for(int i=0; i<C.size(); i++)
      {
        int upC = C2node_up[i];
        int dnC = C2node_dn[i];
        ValueType tmp = C[i]*detValues_up[upC]*detValues_dn[dnC]*psiinv;
        lapl_sum += tmp*(tempstorage_up[upC]+tempstorage_dn[dnC]+static_cast<ValueType>(2.0*Dot(grads_up[upC],grads_dn[dnC])));
        g += grads_up[upC]*static_cast<ParticleSet::ParticleValue_t>(tmp);
        g += grads_dn[dnC]*static_cast<ParticleSet::ParticleValue_t>(tmp);
      }
      gmP=g-P.G;
      ggP=Dot(gmP,g);
      // CI coefficients
      if(recalculate)
      {
        int ip(1);
        for(int i=0; i<C.size()-1; i++,ip++)
        {
          int kk=myVars.where(i);
          if (kk<0)
            continue;
          //dlogpsi[kk] = cdet;
          int upC = C2node_up[ip];
          int dnC = C2node_dn[ip];
          ValueType cdet=detValues_up[upC]*detValues_dn[dnC]*psiinv;
          convert(cdet,dlogpsi[kk]);
          ValueType dhpsi =  ((RealType)-0.5*cdet)*
                             ( tempstorage_up[upC]+tempstorage_dn[dnC]-lapl_sum
                               +static_cast<ValueType>(2.0*Dot(grads_up[upC],grads_dn[dnC]))
                               +(RealType)2.0*(ggP-static_cast<ValueType>(Dot(gmP,grads_up[upC])+Dot(gmP,grads_dn[dnC]))));
          convert(dhpsi,dhpsioverpsi[kk]);
        }
      }
      // BF parameters
      if(optmBF)
      {
        // build QP,Amat,Bmat_full,Xmat,Cmat,Ymat
        BFTrans->evaluateDerivatives(P);
        int numBFprm = BFTrans->optIndexMap.size();
        if(dpsia_up.rows() < dets_up.size() || dpsia_up.cols() < numBFprm)
          dpsia_up.resize(dets_up.size(),numBFprm);
        if(dpsia_dn.rows() < dets_dn.size() || dpsia_dn.cols() < numBFprm)
          dpsia_dn.resize(dets_dn.size(),numBFprm);
        if(dLa_up.rows() < dets_up.size() || dLa_up.cols() < numBFprm)
          dLa_up.resize(dets_up.size(),numBFprm);
        if(dLa_dn.rows() < dets_dn.size() || dLa_dn.cols() < numBFprm)
          dLa_dn.resize(dets_dn.size(),numBFprm);
        if(dGa_up.size(0) < dets_up.size() || dGa_up.size(1) < numBFprm || dGa_up.size(2) < n)
          dGa_up.resize(dets_up.size(),numBFprm,n);
        if(dGa_dn.size(0) < dets_dn.size() || dGa_dn.size(1) < numBFprm || dGa_dn.size(2) < n)
          dGa_dn.resize(dets_dn.size(),numBFprm,n);
        // avoid this in the future
        spo_up->evaluateForWalkerMoveWithThirdDeriv(BFTrans->QP,FirstIndex_up,LastIndex_up);
        spo_dn->evaluateForWalkerMoveWithThirdDeriv(BFTrans->QP,FirstIndex_dn,LastIndex_dn);
        for(int i=0; i<dets_up.size(); i++)
        {
          spo_up->prepareFor(i);
          dets_up[i]->evaluateDerivatives(BFTrans->QP,optvars,i,dpsia_up,dGa_up,dLa_up);
        }
        for(int i=0; i<dets_dn.size(); i++)
        {
          spo_dn->prepareFor(i);
          dets_dn[i]->evaluateDerivatives(BFTrans->QP,optvars,i,dpsia_dn,dGa_dn,dLa_dn);
        }
        for(int pa=0; pa<numBFprm; pa++)
        {
          int kk=BFTrans->optIndexMap[pa];
          if (kk<0)
            continue;
          ValueType dlog=0.0,dhpsi = 0.0;
          for(int i=0; i<C.size(); i++)
          {
            int upC = C2node_up[i];
            int dnC = C2node_dn[i];
            ValueType cdet=C[i]*detValues_up[upC]*detValues_dn[dnC]*psiinv;
            ParticleSet::ParticleValue_t dot1=0.0;
            ValueType dpsi1=dpsia_up(upC,pa);
            ValueType dpsi2=dpsia_dn(dnC,pa);
            ParticleSet::ParticleGradient_t& g1 = grads_up[upC];
            ParticleSet::ParticleGradient_t& g2 = grads_dn[dnC];
            for(int k=0; k<n; k++)
              dot1 += dot((g2[k]-gmP[k]),dGa_up(upC,pa,k))
                      + dot((g1[k]-gmP[k]),dGa_dn(dnC,pa,k))
                      - static_cast<ParticleSet::ParticleValue_t>(dpsi1)*(dot(gmP[k],g2[k]))
                      - static_cast<ParticleSet::ParticleValue_t>(dpsi2)*(dot(gmP[k],g1[k]));
            dlog += cdet*(dpsi1+dpsi2);
            dhpsi += cdet*( dLa_up(upC,pa) + dLa_dn(dnC,pa)
                            + dpsi2*tempstorage_up[upC] + dpsi1*tempstorage_dn[dnC]
                            + static_cast<ValueType>(2.0*dot1));
          } // i
          dhpsi = (RealType)-0.5*(dhpsi + dlog*((RealType)2.0*ggP-lapl_sum));
          convert(dlog,dlogpsi[kk]);
          convert(dhpsi,dhpsioverpsi[kk]);
        }  // pa
      }
    }
  }
}

}
