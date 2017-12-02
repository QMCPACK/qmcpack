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
    
    
#include "QMCWaveFunctions/Fermion/MultiSlaterDeterminant.h"
#include "ParticleBase/ParticleAttribOps.h"

namespace qmcplusplus
{

MultiSlaterDeterminant::MultiSlaterDeterminant(ParticleSet& targetPtcl, SPOSetProxyPtr upspo, SPOSetProxyPtr dnspo):spo_up(upspo),spo_dn(dnspo),
  RatioTimer("MultiSlaterDeterminant::ratio"),
  Ratio1Timer("MultiSlaterDeterminant::detEval_ratio"),
  RatioGradTimer("MultiSlaterDeterminant::ratioGrad"),
  Ratio1GradTimer("MultiSlaterDeterminant::detEval_ratioGrad"),
  RatioAllTimer("MultiSlaterDeterminant::ratio(all)"),
  Ratio1AllTimer("MultiSlaterDeterminant::detEval_ratio(all)"),
  UpdateTimer("MultiSlaterDeterminant::updateBuffer"),
  AccRejTimer("MultiSlaterDeterminant::Accept_Reject"),
  EvaluateTimer("MultiSlaterDeterminant::evaluate"),
  evalOrbTimer("MultiSlaterDeterminant::evalOrbGrad")
{
  registerTimers();
  //Optimizable=true;
  Optimizable=false;
  OrbitalName="MultiSlaterDeterminant";
  usingCSF=false;
  FirstIndex_up = targetPtcl.first(0);
  LastIndex_up = targetPtcl.last(0);
  FirstIndex_dn = targetPtcl.first(1);
  LastIndex_dn = targetPtcl.last(1);
  nels_up = LastIndex_up-FirstIndex_up;
  nels_dn = LastIndex_dn-FirstIndex_dn;
  DetID.resize(targetPtcl.getTotalNum());
  for(int i=0; i<targetPtcl.groups(); ++i)
    for(int j=targetPtcl.first(i); j<targetPtcl.last(i); ++j)
      DetID[j]=i;
}

OrbitalBasePtr MultiSlaterDeterminant::makeClone(ParticleSet& tqp) const
{
  SPOSetProxyForMSD* spo_up_C = new SPOSetProxyForMSD(spo_up->refPhi->makeClone(),FirstIndex_up,LastIndex_up);
  SPOSetProxyForMSD* spo_dn_C = new SPOSetProxyForMSD(spo_dn->refPhi->makeClone(),FirstIndex_dn,LastIndex_dn);
  spo_up_C->occup= spo_up->occup;
  spo_dn_C->occup= spo_dn->occup;
  spo_up_C->refPhi->resetTargetParticleSet(tqp);
  spo_dn_C->refPhi->resetTargetParticleSet(tqp);
  MultiSlaterDeterminant* clone = new MultiSlaterDeterminant(tqp,spo_up_C,spo_dn_C);
  clone->C2node_up=C2node_up;
  clone->C2node_dn=C2node_dn;
  clone->resize(dets_up.size(),dets_dn.size());
  if (usingCSF)
  {
    clone->CSFcoeff=CSFcoeff;
    clone->CSFexpansion=CSFexpansion;
    clone->DetsPerCSF=DetsPerCSF;
  }
//     SPOSetProxyForMSD* spo = clone->spo_up;
//     spo->occup.resize(uniqueConfg_up.size(),clone->nels_up);
  for(int i=0; i<dets_up.size(); i++)
  {
//       int nq=0;
// //       configuration& ci = uniqueConfg_up[i];
//       for(int k=0; k<uniqueConfg_up[i].occup.size(); k++) {
//         if(uniqueConfg_up[i].occup[k]) {
//           spo->occup(i,nq++) = k;
//         }
//       }
    DiracDeterminantBase* adet = new DiracDeterminantBase((SPOSetBasePtr) clone->spo_up,0);
    adet->set(clone->FirstIndex_up,clone->nels_up);
    adet->resetTargetParticleSet(tqp);
    clone->dets_up.push_back(adet);
  }
//     spo = clone->spo_dn;
//     spo->occup.resize(uniqueConfg_dn.size(),clone->nels_dn);
  for(int i=0; i<dets_dn.size(); i++)
  {
//       int nq=0;
// //       configuration& ci = uniqueConfg_dn[i];
//       for(int k=0; k<uniqueConfg_dn[i].occup.size(); k++) {
//         if(uniqueConfg_dn[i].occup[k]) {
//           spo->occup(i,nq++) = k;
//         }
//       }
    DiracDeterminantBase* adet = new DiracDeterminantBase((SPOSetBasePtr) clone->spo_dn,0);
    adet->set(clone->FirstIndex_dn,clone->nels_dn);
    adet->resetTargetParticleSet(tqp);
    clone->dets_dn.push_back(adet);
  }
  clone->Optimizable=Optimizable;
  clone->C=C;
  clone->myVars=myVars;
  return clone;
}


MultiSlaterDeterminant::~MultiSlaterDeterminant() { }
void MultiSlaterDeterminant::resetTargetParticleSet(ParticleSet& P)
{
  for(int i=0; i<dets_up.size(); i++)
    dets_up[i]->resetTargetParticleSet(P);
  for(int i=0; i<dets_dn.size(); i++)
    dets_dn[i]->resetTargetParticleSet(P);
}

void MultiSlaterDeterminant::resize(int n1, int n2)
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

OrbitalBase::ValueType MultiSlaterDeterminant::evaluate(ParticleSet& P
    , ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
{
  EvaluateTimer.start();
// mmorales: For now always assume 2 electron components, up/down
  spo_up->evaluateForWalkerMove(P,FirstIndex_up,LastIndex_up);
  spo_dn->evaluateForWalkerMove(P,FirstIndex_dn,LastIndex_dn);
  for(int i=0; i<dets_up.size(); i++)
  {
    spo_up->prepareFor(i);
    grads_up[i]=0.0;
    lapls_up[i]=0.0;
    detValues_up[i]=dets_up[i]->evaluate(P,grads_up[i],lapls_up[i]);
    // need \nabla^2 Det / Det
    for(int k=FirstIndex_up; k<LastIndex_up; k++)
      lapls_up[i][k] += dot(grads_up[i][k],grads_up[i][k]);
  }
  for(int i=0; i<dets_dn.size(); i++)
  {
    spo_dn->prepareFor(i);
    grads_dn[i]=0.0;
    lapls_dn[i]=0.0;
    detValues_dn[i]=dets_dn[i]->evaluate(P,grads_dn[i],lapls_dn[i]);
    // need \nabla^2 Det / Det
    for(int k=FirstIndex_dn; k<LastIndex_dn; k++)
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
    //for(int n=FirstIndex_up; n<LastIndex_up; n++) {
    myG += grads_up[upC]*tmp; // other spin sector should be zero
    myL += lapls_up[upC]*tmp;
    //}
    //for(int n=FirstIndex_dn; n<LastIndex_dn; n++) {
    myG += grads_dn[dnC]*tmp; // other spin sector should be zero
    myL += lapls_dn[dnC]*tmp;
    //}
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

OrbitalBase::RealType MultiSlaterDeterminant::evaluateLog(ParticleSet& P
    , ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
{
  ValueType psi = evaluate(P,G,L);
  return LogValue = evaluateLogAndPhase(psi,PhaseValue);
}

OrbitalBase::GradType MultiSlaterDeterminant::evalGrad(ParticleSet& P, int iat)
{
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
      ValueType tmp = C[i]*detValues_up[upC]*detValues_dn[dnC];
      psi += tmp;
      GradType grad_temp = grads_up[upC][iat];
      grad_iat += grad_temp*tmp;
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
      ValueType tmp = C[i]*detValues_up[upC]*detValues_dn[dnC];
      psi += tmp;
      GradType grad_temp = grads_dn[dnC][iat];
      grad_iat += grad_temp*tmp;
    }
    grad_iat *= (RealType)1.0/psi;
    return grad_iat;
  }
}

OrbitalBase::ValueType MultiSlaterDeterminant::ratioGrad(ParticleSet& P
    , int iat, GradType& grad_iat)
{
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
OrbitalBase::ValueType MultiSlaterDeterminant::ratio(ParticleSet& P, int iat)
{
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

void MultiSlaterDeterminant::acceptMove(ParticleSet& P, int iat)
{
// this should depend on the type of update, ratio / ratioGrad
// for now is incorrect fot ratio(P,iat,dG,dL) updates
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

void MultiSlaterDeterminant::restore(int iat)
{
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

void MultiSlaterDeterminant::registerData(ParticleSet& P, WFBufferType& buf)
{
// move resize of pbyp structures to here
  spo_up->evaluateForWalkerMove(P,FirstIndex_up,LastIndex_up);
  spo_dn->evaluateForWalkerMove(P,FirstIndex_dn,LastIndex_dn);
  myG = P.G;
  myL = P.L;
  RealType logpsi(0.0);
  PhaseValue=0.0;
  for (int i=0; i<dets_up.size(); i++)
  {
    spo_up->prepareFor(i);
    dets_up[i]->registerData(P,buf);
  }
  for (int i=0; i<dets_dn.size(); i++)
  {
    spo_dn->prepareFor(i);
    dets_dn[i]->registerData(P,buf);
  }
  P.G = myG;
  P.L = myL;
  logpsi = evaluateLog(P,P.G,P.L);
  int TotalDim = PosType::Size*P.getTotalNum();
  //buf.add(detValues_up.begin(),detValues_up.end());
  //buf.add(detValues_dn.begin(),detValues_dn.end());
  buf.add(detValues_up.first_address(),detValues_up.last_address());
  buf.add(detValues_dn.first_address(),detValues_dn.last_address());
  for(int i=0; i<NumUniqueDets_up; i++)
  {
    buf.add(&(grads_up[i][0][0]), &(grads_up[i][0][0])+TotalDim);
//      buf.add(&(lapls_up[i][0]), &(lapls_up[i][P.getTotalNum()]));
    buf.add(lapls_up[i].first_address(),lapls_up[i].last_address());
  }
  for(int i=0; i<NumUniqueDets_dn; i++)
  {
    buf.add(&(grads_dn[i][0][0]), &(grads_dn[i][0][0])+TotalDim);
//      buf.add(&(lapls_dn[i][0]), &(lapls_dn[i][P.getTotalNum()]));
    buf.add(lapls_dn[i].first_address(),lapls_dn[i].last_address());
  }
}

// FIX FIX FIX
OrbitalBase::RealType MultiSlaterDeterminant::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch)
{
  UpdateTimer.start();
  if(fromscratch || UpdateMode == ORB_PBYP_RATIO)
  {
    spo_up->evaluateForWalkerMove(P,FirstIndex_up,LastIndex_up);
    spo_dn->evaluateForWalkerMove(P,FirstIndex_dn,LastIndex_dn);
  }
  myG = P.G;
  myL = P.L;
  RealType logpsi(0.0);
  PhaseValue=0.0;
  for (int i=0; i<dets_up.size(); i++)
  {
    spo_up->prepareFor(i);
    logpsi = dets_up[i]->updateBuffer(P,buf,fromscratch);
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
    logpsi = dets_dn[i]->updateBuffer(P,buf,fromscratch);
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
  //buf.put(detValues_up.begin(),detValues_up.end());
  //buf.put(detValues_dn.begin(),detValues_dn.end());
  buf.put(detValues_up.first_address(),detValues_up.last_address());
  buf.put(detValues_dn.first_address(),detValues_dn.last_address());
  for(int i=0; i<NumUniqueDets_up; i++)
  {
    buf.put(&(grads_up[i][0][0]), &(grads_up[i][0][0])+TotalDim);
//      buf.put(&(lapls_up[i][0]), &(lapls_up[i][P.getTotalNum()]));
    buf.put(lapls_up[i].first_address(),lapls_up[i].last_address());
  }
  for(int i=0; i<NumUniqueDets_dn; i++)
  {
    buf.put(&(grads_dn[i][0][0]), &(grads_dn[i][0][0])+TotalDim);
//      buf.put(&(lapls_dn[i][0]), &(lapls_dn[i][P.getTotalNum()]));
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
    //for(int n=FirstIndex_up; n<LastIndex_up; n++) {
    myG += grads_up[upC]*tmp; // other spin sector should be zero
    myL += lapls_up[upC]*tmp;
    //}
    //for(int n=FirstIndex_dn; n<LastIndex_dn; n++) {
    myG += grads_dn[dnC]*tmp; // other spin sector should be zero
    myL += lapls_dn[dnC]*tmp;
    //}
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

void MultiSlaterDeterminant::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  for (int i=0; i<dets_up.size(); i++)
    dets_up[i]->copyFromBuffer(P,buf);
  for (int i=0; i<dets_dn.size(); i++)
    dets_dn[i]->copyFromBuffer(P,buf);
  int TotalDim = PosType::Size*P.getTotalNum();
  buf.get(detValues_up.begin(),detValues_up.end());
  buf.get(detValues_dn.begin(),detValues_dn.end());
  for(int i=0; i<NumUniqueDets_up; i++)
  {
    buf.get(&(grads_up[i][0][0]), &(grads_up[i][0][0])+TotalDim);
    //buf.get(&(lapls_up[i][0]), &(lapls_up[i][P.getTotalNum()]));
    buf.get(lapls_up[i].first_address(),lapls_up[i].last_address());
  }
  for(int i=0; i<NumUniqueDets_dn; i++)
  {
    buf.get(&(grads_dn[i][0][0]), &(grads_dn[i][0][0])+TotalDim);
    //buf.get(&(lapls_dn[i][0]), &(lapls_dn[i][P.getTotalNum()]));
    buf.get(lapls_dn[i].first_address(),lapls_dn[i].last_address());
  }
}


void MultiSlaterDeterminant::checkInVariables(opt_variables_type& active)
{
  if(Optimizable)
  {
    if(myVars.size())
      active.insertFrom(myVars);
    else
      Optimizable=false;
  }
}

void MultiSlaterDeterminant::checkOutVariables(const opt_variables_type& active)
{
  if(Optimizable)
    myVars.getIndex(active);
}

/** resetParameters with optVariables
 *
 * USE_resetParameters
 */
void MultiSlaterDeterminant::resetParameters(const opt_variables_type& active)
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
    //for(int i=0; i<SDets.size(); i++) SDets[i]->resetParameters(active);
  }
}
void MultiSlaterDeterminant::reportStatus(std::ostream& os)
{
}

//   OrbitalBasePtr MultiSlaterDeterminant::makeClone(ParticleSet& tqp) const
//   {
//      APP_ABORT("IMPLEMENT OrbitalBase::makeClone");
//      return 0;
//   }

void MultiSlaterDeterminant::evaluateDerivatives(ParticleSet& P,
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
// need to modify for CSF later on, right now assume Slater Det basis
  if (recalculate)
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
        lapl_sum += tmp*(tempstorage_up[upC]+tempstorage_dn[dnC]);
        ParticleSet::ParticleValue_t tmp_DP(tmp);
        g += grads_up[upC]*tmp_DP;
        g += grads_dn[dnC]*tmp_DP;
      }
      gmP=g-P.G;
      gg=Dot(gmP,g);
      //gg=Dot(g,g)-Dot(P.G,g);
//        ggP=Dot(P.G,g);
      int num=CSFcoeff.size()-1;
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
          ValueType tmp=CSFexpansion[cnt]*detValues_up[upC]*detValues_dn[dnC]*psiinv;
          cdet+=tmp;
          q0 += tmp*(tempstorage_up[upC]+tempstorage_dn[dnC]);
          ParticleSet::ParticleValue_t tmp_DP(tmp);
          v1 += tmp_DP*(Dot(gmP,grads_up[upC])+Dot(gmP,grads_dn[dnC]));
//           v1 += tmp*(Dot(P.G,grads_up[upC])-Dot(g,grads_up[upC]));
//           v2 += tmp*(Dot(P.G,grads_dn[dnC])-Dot(g,grads_dn[dnC]));
          cnt++;
        }
        convert(cdet,dlogpsi[kk]);
        ValueType dhpsi =  (RealType)(-0.5)*(q0-cdet*lapl_sum)
                           -cdet*gg+v1;
//                            -cdet*gg-v1-v2;
        //ValueType dhpsi =  -0.5*(tmp1*laplSum_up[upC]+tmp2*laplSum_dn[dnC]
        //                         -cdet*lapl_sum)
        //                   -cdet*gg-(tmp1*v1+tmp2*v2);
        convert(dhpsi,dhpsioverpsi[kk]);
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
        lapl_sum += tmp*(tempstorage_up[upC]+tempstorage_dn[dnC]);
        ParticleSet::ParticleValue_t tmp_DP(tmp);
        g += grads_up[upC]*tmp_DP;
        g += grads_dn[dnC]*tmp_DP;
      }
      gmP=g-P.G;
      gg=Dot(gmP,g);
      //gg=Dot(g,g);
      //ggP=Dot(P.G,g);
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
        ValueType dhpsi =  ((RealType)(-0.5)*cdet)*
                           ( tempstorage_up[upC]+tempstorage_dn[dnC]-lapl_sum
                             +(RealType)2.0*(gg-static_cast<ValueType>(Dot(gmP,grads_up[upC])+Dot(gmP,grads_dn[dnC]))));
        //+2.0*(gg-Dot(g,grads_up[upC])-Dot(g,grads_dn[dnC])
        //+Dot(P.G,grads_up[upC])+Dot(P.G,grads_dn[dnC])-ggP));
        convert(dhpsi,dhpsioverpsi[kk]);
      }
    }
  }
}

void MultiSlaterDeterminant::registerTimers()
{
  RatioTimer.reset();
  Ratio1Timer.reset();
  RatioGradTimer.reset();
  Ratio1GradTimer.reset();
  RatioAllTimer.reset();
  Ratio1AllTimer.reset();
  UpdateTimer.reset();
  AccRejTimer.reset();
  EvaluateTimer.reset();
  evalOrbTimer.reset();
  TimerManager.addTimer (&RatioTimer);
  TimerManager.addTimer (&RatioGradTimer);
  TimerManager.addTimer (&RatioAllTimer);
  TimerManager.addTimer (&Ratio1Timer);
  TimerManager.addTimer (&Ratio1GradTimer);
  TimerManager.addTimer (&Ratio1AllTimer);
  TimerManager.addTimer (&UpdateTimer);
  TimerManager.addTimer (&AccRejTimer);
  TimerManager.addTimer (&EvaluateTimer);
  TimerManager.addTimer (&evalOrbTimer);
}

}
