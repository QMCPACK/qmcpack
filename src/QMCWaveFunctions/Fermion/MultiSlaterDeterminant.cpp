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


#include "MultiSlaterDeterminant.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "type_traits/ConvertToReal.h"

namespace qmcplusplus
{
MultiSlaterDeterminant::MultiSlaterDeterminant(ParticleSet& targetPtcl,
                                               std::vector<std::unique_ptr<SPOSetProxyForMSD>> spos,
                                               const std::string& class_name)
    : WaveFunctionComponent(class_name),
      RatioTimer(*timer_manager.createTimer(ClassName + "ratio")),
      RatioGradTimer(*timer_manager.createTimer(ClassName + "ratioGrad")),
      RatioAllTimer(*timer_manager.createTimer(ClassName + "ratio(all)")),
      UpdateTimer(*timer_manager.createTimer(ClassName + "updateBuffer")),
      EvaluateTimer(*timer_manager.createTimer(ClassName + "evaluate")),
      Ratio1Timer(*timer_manager.createTimer(ClassName + "detEval_ratio")),
      Ratio1GradTimer(*timer_manager.createTimer(ClassName + "detEval_ratioGrad")),
      Ratio1AllTimer(*timer_manager.createTimer(ClassName + "detEval_ratio(all)")),
      AccRejTimer(*timer_manager.createTimer(ClassName + "Accept_Reject")),
      evalOrbTimer(*timer_manager.createTimer(ClassName + "evalOrbGrad")),
      spo_up(std::move(spos[0])),
      spo_dn(std::move(spos[1]))
{
  assert(spos.size() == targetPtcl.groups());
  registerTimers();

  Optimizable   = false;
  is_fermionic  = true;
  usingCSF      = false;
  FirstIndex_up = targetPtcl.first(0);
  LastIndex_up  = targetPtcl.last(0);
  FirstIndex_dn = targetPtcl.first(1);
  LastIndex_dn  = targetPtcl.last(1);
  nels_up       = LastIndex_up - FirstIndex_up;
  nels_dn       = LastIndex_dn - FirstIndex_dn;
  DetID.resize(targetPtcl.getTotalNum());
  for (int i = 0; i < targetPtcl.groups(); ++i)
    for (int j = targetPtcl.first(i); j < targetPtcl.last(i); ++j)
      DetID[j] = i;
}

std::unique_ptr<WaveFunctionComponent> MultiSlaterDeterminant::makeClone(ParticleSet& tqp) const
{
  typedef DiracDeterminant<> Det_t;
  auto spo_up_C    = std::make_unique<SPOSetProxyForMSD>(std::unique_ptr<SPOSet>(spo_up->refPhi->makeClone()),
                                                      FirstIndex_up, LastIndex_up);
  auto spo_dn_C    = std::make_unique<SPOSetProxyForMSD>(std::unique_ptr<SPOSet>(spo_dn->refPhi->makeClone()),
                                                      FirstIndex_dn, LastIndex_dn);
  spo_up_C->occup  = spo_up->occup;
  spo_dn_C->occup  = spo_dn->occup;
  std::vector<std::unique_ptr<SPOSetProxyForMSD>> spos;
  spos.push_back(std::move(spo_up_C));
  spos.push_back(std::move(spo_dn_C));
  auto clone       = std::make_unique<MultiSlaterDeterminant>(tqp, std::move(spos));
  clone->C2node_up = C2node_up;
  clone->C2node_dn = C2node_dn;
  clone->resize(dets_up.size(), dets_dn.size());
  if (usingCSF)
  {
    clone->CSFcoeff     = CSFcoeff;
    clone->CSFexpansion = CSFexpansion;
    clone->DetsPerCSF   = DetsPerCSF;
  }

  for (int i = 0; i < dets_up.size(); i++)
    clone->dets_up.push_back(std::make_unique<Det_t>(std::static_pointer_cast<SPOSet>(clone->spo_up), clone->FirstIndex_up, clone->FirstIndex_up + clone->nels_up));
  for (int i = 0; i < dets_dn.size(); i++)
    clone->dets_dn.emplace_back(std::make_unique<Det_t>(std::static_pointer_cast<SPOSet>(clone->spo_dn), clone->FirstIndex_dn, clone->FirstIndex_dn + clone->nels_dn));

  clone->Optimizable = Optimizable;
  clone->C           = C;
  clone->myVars      = myVars;
  return clone;
}


MultiSlaterDeterminant::~MultiSlaterDeterminant() {}

void MultiSlaterDeterminant::resize(int n1, int n2)
{
  NumUniqueDets_up = n1;
  NumUniqueDets_dn = n2;
  myG.resize(nels_up + nels_dn);
  myL.resize(nels_up + nels_dn);
  grads_up.resize(NumUniqueDets_up);
  grads_dn.resize(NumUniqueDets_dn);
  lapls_up.resize(NumUniqueDets_up);
  lapls_dn.resize(NumUniqueDets_dn);
  for (int i = 0; i < NumUniqueDets_up; i++)
  {
    grads_up[i].resize(nels_up + nels_dn);
    lapls_up[i].resize(nels_up + nels_dn);
  }
  for (int i = 0; i < NumUniqueDets_dn; i++)
  {
    grads_dn[i].resize(nels_up + nels_dn);
    lapls_dn[i].resize(nels_up + nels_dn);
  }
  detValues_up.resize(NumUniqueDets_up);
  detValues_dn.resize(NumUniqueDets_dn);
  int maxDet = NumUniqueDets_up;
  if (NumUniqueDets_dn > NumUniqueDets_up)
    maxDet = NumUniqueDets_dn;
  detsRatios.resize(maxDet);
  tempstorage_up.resize(NumUniqueDets_up);
  tempstorage_dn.resize(NumUniqueDets_dn);
  grad_temp.resize(maxDet);
  tempgrad.resize(maxDet);
  templapl.resize(maxDet);
  for (int i = 0; i < maxDet; i++)
  {
    tempgrad[i].resize(nels_up + nels_dn);
    templapl[i].resize(nels_up + nels_dn);
  }
}

WaveFunctionComponent::ValueType MultiSlaterDeterminant::evaluate(const ParticleSet& P,
                                                                  ParticleSet::ParticleGradient_t& G,
                                                                  ParticleSet::ParticleLaplacian_t& L)
{
  EvaluateTimer.start();
  // mmorales: For now always assume 2 electron components, up/down
  spo_up->evaluateForWalkerMove(P, FirstIndex_up, LastIndex_up);
  spo_dn->evaluateForWalkerMove(P, FirstIndex_dn, LastIndex_dn);
  for (int i = 0; i < dets_up.size(); i++)
  {
    spo_up->prepareFor(i);
    grads_up[i] = 0.0;
    lapls_up[i] = 0.0;
    dets_up[i]->evaluateLog(P, grads_up[i], lapls_up[i]);
    detValues_up[i] = dets_up[i]->getValue();
    // need \nabla^2 Det / Det
    for (int k = FirstIndex_up; k < LastIndex_up; k++)
      lapls_up[i][k] += dot(grads_up[i][k], grads_up[i][k]);
  }
  for (int i = 0; i < dets_dn.size(); i++)
  {
    spo_dn->prepareFor(i);
    grads_dn[i] = 0.0;
    lapls_dn[i] = 0.0;
    dets_dn[i]->evaluateLog(P, grads_dn[i], lapls_dn[i]);
    detValues_dn[i] = dets_dn[i]->getValue();
    // need \nabla^2 Det / Det
    for (int k = FirstIndex_dn; k < LastIndex_dn; k++)
      lapls_dn[i][k] += dot(grads_dn[i][k], grads_dn[i][k]);
  }
  ValueType psi = 0.0;
  myG           = 0.0;
  myL           = 0.0;
  for (int i = 0; i < C.size(); i++)
  {
    int upC                                = C2node_up[i];
    int dnC                                = C2node_dn[i];
    ParticleSet::SingleParticleValue_t tmp = C[i] * detValues_up[upC] * detValues_dn[dnC];
    psi += tmp;
    //for(int n=FirstIndex_up; n<LastIndex_up; n++) {
    myG += grads_up[upC] * tmp; // other spin sector should be zero
    myL += lapls_up[upC] * tmp;
    //}
    //for(int n=FirstIndex_dn; n<LastIndex_dn; n++) {
    myG += grads_dn[dnC] * tmp; // other spin sector should be zero
    myL += lapls_dn[dnC] * tmp;
    //}
  }
  ValueType psiinv = (RealType)1.0 / psi;
  myG *= psiinv;
  myL *= psiinv;
  G += myG;
  for (int i = 0; i < L.size(); i++)
    L[i] += myL[i] - dot(myG[i], myG[i]);
  EvaluateTimer.stop();
  return psi;
}

WaveFunctionComponent::LogValueType MultiSlaterDeterminant::evaluateLog(const ParticleSet& P,
                                                                        ParticleSet::ParticleGradient_t& G,
                                                                        ParticleSet::ParticleLaplacian_t& L)
{
  return log_value_ = convertValueToLog(evaluate(P, G, L));
}

WaveFunctionComponent::GradType MultiSlaterDeterminant::evalGrad(ParticleSet& P, int iat)
{
  GradType grad_iat;
  if (DetID[iat] == 0)
  {
    for (int i = 0; i < dets_up.size(); i++)
    {
      spo_up->prepareFor(i);
      grads_up[i][iat] = dets_up[i]->evalGrad(P, iat);
    }
    ValueType psi = 0.0;
    for (int i = 0; i < C.size(); i++)
    {
      int upC       = C2node_up[i];
      int dnC       = C2node_dn[i];
      ValueType tmp = C[i] * detValues_up[upC] * detValues_dn[dnC];
      psi += tmp;
      GradType grad_temp = grads_up[upC][iat];
      grad_iat += grad_temp * tmp;
    }
    grad_iat *= (RealType)1.0 / psi;
    return grad_iat;
  }
  else
  {
    ValueType psi = 0.0;
    for (int i = 0; i < dets_dn.size(); i++)
    {
      spo_dn->prepareFor(i);
      grads_dn[i][iat] = dets_dn[i]->evalGrad(P, iat);
    }
    for (int i = 0; i < C.size(); i++)
    {
      int upC       = C2node_up[i];
      int dnC       = C2node_dn[i];
      ValueType tmp = C[i] * detValues_up[upC] * detValues_dn[dnC];
      psi += tmp;
      GradType grad_temp = grads_dn[dnC][iat];
      grad_iat += grad_temp * tmp;
    }
    grad_iat *= (RealType)1.0 / psi;
    return grad_iat;
  }
}

WaveFunctionComponent::PsiValueType MultiSlaterDeterminant::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  UpdateMode = ORB_PBYP_PARTIAL;
  if (DetID[iat] == 0)
  {
    RatioGradTimer.start();
    evalOrbTimer.start();
    spo_up->evaluateAllForPtclMove(P, iat);
    evalOrbTimer.stop();
    Ratio1GradTimer.start();
    grad_temp = 0.0;
    for (int i = 0; i < dets_up.size(); i++)
    {
      spo_up->prepareFor(i);
      detsRatios[i] = dets_up[i]->ratioGrad(P, iat, grad_temp[i]);
    }
    Ratio1GradTimer.stop();
    ValueType psiOld = 0.0, psiNew = 0.0;
    GradType dummy;
    for (int i = 0; i < C.size(); i++)
    {
      int upC        = C2node_up[i];
      int dnC        = C2node_dn[i];
      ValueType tmp2 = C[i] * detValues_up[upC] * detValues_dn[dnC];
      ValueType tmp  = tmp2 * detsRatios[upC];
      psiOld += tmp2;
      psiNew += tmp;
      dummy += tmp * grad_temp[upC];
    }
    grad_iat += dummy / psiNew;
    curRatio = psiNew / psiOld;
    RatioGradTimer.stop();
    return curRatio;
  }
  else
  {
    RatioGradTimer.start();
    evalOrbTimer.start();
    spo_dn->evaluateAllForPtclMove(P, iat);
    evalOrbTimer.stop();
    Ratio1GradTimer.start();
    grad_temp = 0.0;
    for (int i = 0; i < dets_dn.size(); i++)
    {
      spo_dn->prepareFor(i);
      detsRatios[i] = dets_dn[i]->ratioGrad(P, iat, grad_temp[i]);
    }
    Ratio1GradTimer.stop();
    ValueType psiOld = 0.0, psiNew = 0.0;
    GradType dummy;
    for (int i = 0; i < C.size(); i++)
    {
      int upC        = C2node_up[i];
      int dnC        = C2node_dn[i];
      ValueType tmp2 = C[i] * detValues_up[upC] * detValues_dn[dnC];
      ValueType tmp  = tmp2 * detsRatios[dnC];
      psiOld += tmp2;
      psiNew += tmp;
      dummy += tmp * grad_temp[dnC];
    }
    grad_iat += dummy / psiNew;
    curRatio = psiNew / psiOld;
    RatioGradTimer.stop();
    return curRatio;
  }
}


// use ci_node for this routine only
WaveFunctionComponent::PsiValueType MultiSlaterDeterminant::ratio(ParticleSet& P, int iat)
{
  UpdateMode = ORB_PBYP_RATIO;
  if (DetID[iat] == 0)
  {
    RatioTimer.start();
    spo_up->evaluateForPtclMove(P, iat);
    Ratio1Timer.start();
    for (int i = 0; i < dets_up.size(); i++)
    {
      spo_up->prepareFor(i);
      detsRatios[i] = dets_up[i]->ratio(P, iat);
    }
    Ratio1Timer.stop();
    std::vector<size_t>::iterator upC(C2node_up.begin()), dnC(C2node_dn.begin());
    std::vector<ValueType>::iterator it(C.begin()), last(C.end());
    ValueType psiOld = 0.0, psiNew = 0.0;
    while (it != last)
    {
      ValueType tmp = (*it) * detValues_up[*upC] * detValues_dn[*dnC];
      psiNew += tmp * detsRatios[*upC];
      psiOld += tmp;
      it++;
      upC++;
      dnC++;
    }
    curRatio = psiNew / psiOld;
    RatioTimer.stop();
    return curRatio;
  }
  else
  {
    RatioTimer.start();
    spo_dn->evaluateForPtclMove(P, iat);
    Ratio1Timer.start();
    for (int i = 0; i < dets_dn.size(); i++)
    {
      spo_dn->prepareFor(i);
      detsRatios[i] = dets_dn[i]->ratio(P, iat);
    }
    Ratio1Timer.stop();
    std::vector<size_t>::iterator upC(C2node_up.begin()), dnC(C2node_dn.begin());
    std::vector<ValueType>::iterator it(C.begin()), last(C.end());
    ValueType psiOld = 0.0, psiNew = 0.0;
    while (it != last)
    {
      ValueType tmp = (*it) * detValues_up[*upC] * detValues_dn[*dnC];
      psiNew += tmp * detsRatios[*dnC];
      psiOld += tmp;
      it++;
      upC++;
      dnC++;
    }
    curRatio = psiNew / psiOld;
    RatioTimer.stop();
    return curRatio;
  }
}

void MultiSlaterDeterminant::acceptMove(ParticleSet& P, int iat, bool safe_to_delay)
{
  // this should depend on the type of update, ratio / ratioGrad
  // for now is incorrect fot ratio(P,iat,dG,dL) updates
  AccRejTimer.start();
  if (DetID[iat] == 0)
  {
    for (int i = 0; i < dets_up.size(); i++)
      dets_up[i]->acceptMove(P, iat);
    switch (UpdateMode)
    {
    case ORB_PBYP_RATIO:
      // ratio(P,iat)
      for (int i = 0; i < detValues_up.size(); i++)
        detValues_up[i] *= detsRatios[i];
      log_value_ += convertValueToLog(curRatio);
      curRatio = 1.0;
      break;
    case ORB_PBYP_PARTIAL:
      // ratioGrad(P,iat,grad)
      for (int i = 0; i < detValues_up.size(); i++)
      {
        detValues_up[i] *= detsRatios[i];
        grads_up[i][iat] = grad_temp[i];
      }
      log_value_ += convertValueToLog(curRatio);
      curRatio = 1.0;
      break;
    case ORB_PBYP_ALL:
      // ratio(P,iat,dG,dL)
      for (int i = 0; i < detValues_up.size(); i++)
      {
        detValues_up[i] *= detsRatios[i];
        grads_up[i] = tempgrad[i];
        lapls_up[i] = templapl[i];
      }
      log_value_ += convertValueToLog(curRatio);
      curRatio = 1.0;
      break;
    default:
      for (int i = 0; i < detValues_up.size(); i++)
        detValues_up[i] *= detsRatios[i];
      log_value_ += convertValueToLog(curRatio);
      curRatio = 1.0;
      break;
    }
  }
  else
  {
    for (int i = 0; i < dets_dn.size(); i++)
      dets_dn[i]->acceptMove(P, iat);
    switch (UpdateMode)
    {
    case ORB_PBYP_RATIO:
      // ratio(P,iat)
      for (int i = 0; i < detValues_dn.size(); i++)
        detValues_dn[i] *= detsRatios[i];
      log_value_ += convertValueToLog(curRatio);
      curRatio = 1.0;
      break;
    case ORB_PBYP_PARTIAL:
      // ratioGrad(P,iat,grad)
      for (int i = 0; i < detValues_dn.size(); i++)
      {
        detValues_dn[i] *= detsRatios[i];
        grads_dn[i][iat] = grad_temp[i];
      }
      log_value_ += convertValueToLog(curRatio);
      curRatio = 1.0;
      break;
    case ORB_PBYP_ALL:
      // ratio(P,iat,dG,dL)
      for (int i = 0; i < detValues_dn.size(); i++)
      {
        detValues_dn[i] *= detsRatios[i];
        grads_dn[i] = tempgrad[i];
        lapls_dn[i] = templapl[i];
      }
      log_value_ += convertValueToLog(curRatio);
      curRatio = 1.0;
      break;
    default:
      for (int i = 0; i < detValues_dn.size(); i++)
        detValues_dn[i] *= detsRatios[i];
      log_value_ += convertValueToLog(curRatio);
      curRatio = 1.0;
      break;
    }
  }
  AccRejTimer.stop();
}

void MultiSlaterDeterminant::restore(int iat)
{
  AccRejTimer.start();
  if (DetID[iat] == 0)
  {
    for (int i = 0; i < dets_up.size(); i++)
      dets_up[i]->restore(iat);
  }
  else
  {
    for (int i = 0; i < dets_dn.size(); i++)
      dets_dn[i]->restore(iat);
  }
  curRatio = 1.0;
  AccRejTimer.stop();
}

void MultiSlaterDeterminant::registerData(ParticleSet& P, WFBufferType& buf)
{
  for (int i = 0; i < dets_up.size(); i++)
    dets_up[i]->registerData(P, buf);
  for (int i = 0; i < dets_dn.size(); i++)
    dets_dn[i]->registerData(P, buf);

  buf.add(detValues_up.data(), detValues_up.end());
  buf.add(detValues_dn.data(), detValues_dn.end());

  for (int i = 0; i < NumUniqueDets_up; i++)
  {
    buf.add(grads_up[i].data(), grads_up[i].end());
    buf.add(lapls_up[i].data(), lapls_up[i].end());
  }
  for (int i = 0; i < NumUniqueDets_dn; i++)
  {
    buf.add(grads_dn[i].data(), grads_dn[i].end());
    buf.add(lapls_dn[i].data(), lapls_dn[i].end());
  }
}

// FIX FIX FIX
WaveFunctionComponent::LogValueType MultiSlaterDeterminant::updateBuffer(ParticleSet& P,
                                                                         WFBufferType& buf,
                                                                         bool fromscratch)
{
  UpdateTimer.start();
  if (fromscratch || UpdateMode == ORB_PBYP_RATIO)
  {
    spo_up->evaluateForWalkerMove(P, FirstIndex_up, LastIndex_up);
    spo_dn->evaluateForWalkerMove(P, FirstIndex_dn, LastIndex_dn);
  }
  myG = P.G;
  myL = P.L;
  LogValueType logpsi(0.0);
  for (int i = 0; i < dets_up.size(); i++)
  {
    P.G = 0.0;
    P.L = 0.0;
    spo_up->prepareFor(i);
    logpsi          = dets_up[i]->updateBuffer(P, buf, fromscratch);
    detValues_up[i] = dets_up[i]->getValue();
    grads_up[i]     = P.G;
    lapls_up[i]     = P.L;
    for (int k = FirstIndex_up; k < LastIndex_up; k++)
      lapls_up[i][k] += dot(grads_up[i][k], grads_up[i][k]);
  }
  for (int i = 0; i < dets_dn.size(); i++)
  {
    P.G = 0.0;
    P.L = 0.0;
    spo_dn->prepareFor(i);
    logpsi          = dets_dn[i]->updateBuffer(P, buf, fromscratch);
    detValues_dn[i] = dets_dn[i]->getValue();
    grads_dn[i]     = P.G;
    lapls_dn[i]     = P.L;
    for (int k = FirstIndex_dn; k < LastIndex_dn; k++)
      lapls_dn[i][k] += dot(grads_dn[i][k], grads_dn[i][k]);
  }
  int TotalDim = PosType::Size * P.getTotalNum();
  //buf.put(detValues_up.begin(),detValues_up.end());
  //buf.put(detValues_dn.begin(),detValues_dn.end());
  buf.put(detValues_up.first_address(), detValues_up.last_address());
  buf.put(detValues_dn.first_address(), detValues_dn.last_address());
  for (int i = 0; i < NumUniqueDets_up; i++)
  {
    buf.put(&(grads_up[i][0][0]), &(grads_up[i][0][0]) + TotalDim);
    //      buf.put(&(lapls_up[i][0]), &(lapls_up[i][P.getTotalNum()]));
    buf.put(lapls_up[i].first_address(), lapls_up[i].last_address());
  }
  for (int i = 0; i < NumUniqueDets_dn; i++)
  {
    buf.put(&(grads_dn[i][0][0]), &(grads_dn[i][0][0]) + TotalDim);
    //      buf.put(&(lapls_dn[i][0]), &(lapls_dn[i][P.getTotalNum()]));
    buf.put(lapls_dn[i].first_address(), lapls_dn[i].last_address());
  }
  P.G           = myG;
  P.L           = myL;
  ValueType psi = 0.0;
  myG           = 0.0;
  myL           = 0.0;
  for (int i = 0; i < C.size(); i++)
  {
    int upC                                = C2node_up[i];
    int dnC                                = C2node_dn[i];
    ParticleSet::SingleParticleValue_t tmp = C[i] * detValues_up[upC] * detValues_dn[dnC];
    psi += tmp;
    //for(int n=FirstIndex_up; n<LastIndex_up; n++) {
    myG += grads_up[upC] * tmp; // other spin sector should be zero
    myL += lapls_up[upC] * tmp;
    //}
    //for(int n=FirstIndex_dn; n<LastIndex_dn; n++) {
    myG += grads_dn[dnC] * tmp; // other spin sector should be zero
    myL += lapls_dn[dnC] * tmp;
    //}
  }
  ValueType psiinv = (RealType)1.0 / psi;
  myG *= psiinv;
  myL *= psiinv;
  P.G += myG;
  for (int i = 0; i < P.L.size(); i++)
    P.L[i] += myL[i] - dot(myG[i], myG[i]);
  UpdateTimer.stop();
  return log_value_ = convertValueToLog(psi);
  ;
}

void MultiSlaterDeterminant::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  for (int i = 0; i < dets_up.size(); i++)
    dets_up[i]->copyFromBuffer(P, buf);
  for (int i = 0; i < dets_dn.size(); i++)
    dets_dn[i]->copyFromBuffer(P, buf);
  int TotalDim = PosType::Size * P.getTotalNum();
  buf.get(detValues_up.begin(), detValues_up.end());
  buf.get(detValues_dn.begin(), detValues_dn.end());
  for (int i = 0; i < NumUniqueDets_up; i++)
  {
    buf.get(&(grads_up[i][0][0]), &(grads_up[i][0][0]) + TotalDim);
    //buf.get(&(lapls_up[i][0]), &(lapls_up[i][P.getTotalNum()]));
    buf.get(lapls_up[i].first_address(), lapls_up[i].last_address());
  }
  for (int i = 0; i < NumUniqueDets_dn; i++)
  {
    buf.get(&(grads_dn[i][0][0]), &(grads_dn[i][0][0]) + TotalDim);
    //buf.get(&(lapls_dn[i][0]), &(lapls_dn[i][P.getTotalNum()]));
    buf.get(lapls_dn[i].first_address(), lapls_dn[i].last_address());
  }
}


void MultiSlaterDeterminant::checkInVariables(opt_variables_type& active)
{
  if (Optimizable)
  {
    if (myVars.size())
      active.insertFrom(myVars);
    else
      Optimizable = false;
  }
}

void MultiSlaterDeterminant::checkOutVariables(const opt_variables_type& active)
{
  if (Optimizable)
    myVars.getIndex(active);
}

/** resetParameters with optVariables
 *
 * USE_resetParameters
 */
void MultiSlaterDeterminant::resetParameters(const opt_variables_type& active)
{
  if (Optimizable)
  {
    if (usingCSF)
    {
      for (int i = 0; i < CSFcoeff.size() - 1; i++)
      {
        int loc = myVars.where(i);
        if (loc >= 0)
          CSFcoeff[i + 1] = myVars[i] = active[loc];
      }
      int cnt = 0;
      for (int i = 0; i < DetsPerCSF.size(); i++)
      {
        for (int k = 0; k < DetsPerCSF[i]; k++)
        {
          C[cnt] = CSFcoeff[i] * CSFexpansion[cnt];
          cnt++;
        }
      }
    }
    else
    {
      for (int i = 0; i < C.size() - 1; i++)
      {
        int loc = myVars.where(i);
        if (loc >= 0)
          C[i + 1] = myVars[i] = active[loc];
      }
    }
    //for(int i=0; i<SDets.size(); i++) SDets[i]->resetParameters(active);
  }
}
void MultiSlaterDeterminant::reportStatus(std::ostream& os) {}

//   WaveFunctionComponentPtr MultiSlaterDeterminant::makeClone(ParticleSet& tqp) const
//   {
//      APP_ABORT("IMPLEMENT WaveFunctionComponent::makeClone");
//      return 0;
//   }

void MultiSlaterDeterminant::evaluateDerivatives(ParticleSet& P,
                                                 const opt_variables_type& optvars,
                                                 std::vector<ValueType>& dlogpsi,
                                                 std::vector<ValueType>& dhpsioverpsi)
{
  bool recalculate(false);
  for (int k = 0; k < myVars.size(); ++k)
  {
    int kk = myVars.where(k);
    if (kk < 0)
      continue;
    if (optvars.recompute(kk))
      recalculate = true;
  }
  // need to modify for CSF later on, right now assume Slater Det basis
  if (recalculate)
  {
    if (usingCSF)
    {
      int n            = P.getTotalNum();
      ValueType psiinv = ValueType(1) / LogToValue<ValueType>::convert(log_value_);

      ValueType lapl_sum = 0.0;
      ParticleSet::ParticleGradient_t g(n), gmP(n);
      ValueType gg = 0.0;
      g            = 0.0;
      gmP          = 0.0;
      for (int i = 0; i < lapls_up.size(); i++)
        tempstorage_up[i] = Sum(lapls_up[i]);
      for (int i = 0; i < lapls_dn.size(); i++)
        tempstorage_dn[i] = Sum(lapls_dn[i]);
      for (int i = 0; i < C.size(); i++)
      {
        int upC       = C2node_up[i];
        int dnC       = C2node_dn[i];
        ValueType tmp = C[i] * detValues_up[upC] * detValues_dn[dnC] * psiinv;
        lapl_sum += tmp * (tempstorage_up[upC] + tempstorage_dn[dnC]);
        ParticleSet::SingleParticleValue_t tmp_DP(tmp);
        g += grads_up[upC] * tmp_DP;
        g += grads_dn[dnC] * tmp_DP;
      }
      gmP = g - P.G;
      gg  = Dot(gmP, g);
      //gg=Dot(g,g)-Dot(P.G,g);
      //        ggP=Dot(P.G,g);
      int num = CSFcoeff.size() - 1;
      int cnt = 0;
      //        this one is not optable
      cnt += DetsPerCSF[0];
      int ip(1);
      for (int i = 0; i < num; i++, ip++)
      {
        int kk = myVars.where(i);
        if (kk < 0)
        {
          cnt += DetsPerCSF[ip];
          continue;
        }
        ValueType cdet = 0.0, q0 = 0.0, v1 = 0.0;
        for (int k = 0; k < DetsPerCSF[ip]; k++)
        {
          int upC       = C2node_up[cnt];
          int dnC       = C2node_dn[cnt];
          ValueType tmp = CSFexpansion[cnt] * detValues_up[upC] * detValues_dn[dnC] * psiinv;
          cdet += tmp;
          q0 += tmp * (tempstorage_up[upC] + tempstorage_dn[dnC]);
          ParticleSet::SingleParticleValue_t tmp_DP(tmp);
          v1 += tmp_DP * (Dot(gmP, grads_up[upC]) + Dot(gmP, grads_dn[dnC]));
          //           v1 += tmp*(Dot(P.G,grads_up[upC])-Dot(g,grads_up[upC]));
          //           v2 += tmp*(Dot(P.G,grads_dn[dnC])-Dot(g,grads_dn[dnC]));
          cnt++;
        }
        dlogpsi[kk] = cdet;
        ValueType dhpsi = (RealType)(-0.5) * (q0 - cdet * lapl_sum) - cdet * gg + v1;
        //                            -cdet*gg-v1-v2;
        //ValueType dhpsi =  -0.5*(tmp1*laplSum_up[upC]+tmp2*laplSum_dn[dnC]
        //                         -cdet*lapl_sum)
        //                   -cdet*gg-(tmp1*v1+tmp2*v2);
        dhpsioverpsi[kk] = dhpsi;
      }
    }
    else
    {
      int n            = P.getTotalNum();
      ValueType psiinv = ValueType(1) / LogToValue<ValueType>::convert(log_value_);

      ValueType lapl_sum = 0.0;
      ParticleSet::ParticleGradient_t g(n), gmP(n);
      ValueType gg = 0.0;
      g            = 0.0;
      gmP          = 0.0;
      for (int i = 0; i < lapls_up.size(); i++)
        tempstorage_up[i] = Sum(lapls_up[i]);
      for (int i = 0; i < lapls_dn.size(); i++)
        tempstorage_dn[i] = Sum(lapls_dn[i]);
      for (int i = 0; i < C.size(); i++)
      {
        int upC       = C2node_up[i];
        int dnC       = C2node_dn[i];
        ValueType tmp = C[i] * detValues_up[upC] * detValues_dn[dnC] * psiinv;
        lapl_sum += tmp * (tempstorage_up[upC] + tempstorage_dn[dnC]);
        ParticleSet::SingleParticleValue_t tmp_DP(tmp);
        g += grads_up[upC] * tmp_DP;
        g += grads_dn[dnC] * tmp_DP;
      }
      gmP = g - P.G;
      gg  = Dot(gmP, g);
      //gg=Dot(g,g);
      //ggP=Dot(P.G,g);
      int ip(1);
      for (int i = 0; i < C.size() - 1; i++, ip++)
      {
        int kk = myVars.where(i);
        if (kk < 0)
          continue;
        //dlogpsi[kk] = cdet;
        int upC        = C2node_up[ip];
        int dnC        = C2node_dn[ip];
        ValueType cdet = detValues_up[upC] * detValues_dn[dnC] * psiinv;
        dlogpsi[kk] = cdet;
        ValueType dhpsi = ((RealType)(-0.5) * cdet) *
            (tempstorage_up[upC] + tempstorage_dn[dnC] - lapl_sum +
             (RealType)2.0 * (gg - static_cast<ValueType>(Dot(gmP, grads_up[upC]) + Dot(gmP, grads_dn[dnC]))));
        //+2.0*(gg-Dot(g,grads_up[upC])-Dot(g,grads_dn[dnC])
        //+Dot(P.G,grads_up[upC])+Dot(P.G,grads_dn[dnC])-ggP));
        dhpsioverpsi[kk] = dhpsi;
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
}

} // namespace qmcplusplus
