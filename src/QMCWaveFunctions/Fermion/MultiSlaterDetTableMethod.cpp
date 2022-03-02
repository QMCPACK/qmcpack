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


#include "MultiSlaterDetTableMethod.h"
#include "QMCWaveFunctions/Fermion/MultiDiracDeterminant.h"
#include "ParticleBase/ParticleAttribOps.h"

namespace qmcplusplus
{
MultiSlaterDetTableMethod::MultiSlaterDetTableMethod(ParticleSet& targetPtcl,
                                                     std::vector<std::unique_ptr<MultiDiracDeterminant>>&& dets,
                                                     bool use_pre_computing)
    : WaveFunctionComponent("MultiSlaterDetTableMethod"),
      RatioTimer(*timer_manager.createTimer(ClassName + "::ratio")),
      MWRatioTimer(*timer_manager.createTimer(ClassName + "::mwratio")),
      OffloadRatioTimer(*timer_manager.createTimer(ClassName + "::offloadRatio")),
      OffloadGradTimer(*timer_manager.createTimer(ClassName + "::offloadGrad")),
      EvalGradTimer(*timer_manager.createTimer(ClassName + "::evalGrad")),
      MWEvalGradTimer(*timer_manager.createTimer(ClassName + "::mwevalGrad")),
      RatioGradTimer(*timer_manager.createTimer(ClassName + "::ratioGrad")),
      MWRatioGradTimer(*timer_manager.createTimer(ClassName + "::mwratioGrad")),
      PrepareGroupTimer(*timer_manager.createTimer(ClassName + "::prepareGroup")),
      UpdateTimer(*timer_manager.createTimer(ClassName + "::updateBuffer")),
      AccRejTimer(*timer_manager.createTimer(ClassName + "::Accept_Reject")),
      EvaluateTimer(*timer_manager.createTimer(ClassName + "::evaluate")),
      CI_Optimizable(false),
      use_pre_computing_(use_pre_computing)
{
  registerTimers();
  //Optimizable=true;
  Optimizable  = false;
  is_fermionic = true;
  Dets         = std::move(dets);
  C_otherDs.resize(Dets.size());
  int NP = targetPtcl.getTotalNum();
  myG.resize(NP);
  myL.resize(NP);
  myG_temp.resize(NP);
  myL_temp.resize(NP);

  Last.resize(targetPtcl.groups());
  for (int i = 0; i < Last.size(); ++i)
    Last[i] = targetPtcl.last(i) - 1;
}

void MultiSlaterDetTableMethod::initialize(std::unique_ptr<std::vector<std::vector<size_t>>> C2node_in,
                                           std::unique_ptr<std::vector<ValueType>> C_in,
                                           std::unique_ptr<opt_variables_type> myVars_in,
                                           std::unique_ptr<CSFData> csf_data_in,
                                           bool optimizable,
                                           bool CI_optimizable)
{
  C2node         = std::move(C2node_in);
  C              = std::move(C_in);
  myVars         = std::move(myVars_in);
  csf_data_      = std::move(csf_data_in);
  Optimizable    = optimizable;
  CI_Optimizable = CI_optimizable;
}

MultiSlaterDetTableMethod::~MultiSlaterDetTableMethod() = default;

std::unique_ptr<WaveFunctionComponent> MultiSlaterDetTableMethod::makeClone(ParticleSet& tqp) const
{
  std::vector<std::unique_ptr<MultiDiracDeterminant>> dets_clone;
  for (auto& det : Dets)
    dets_clone.emplace_back(std::make_unique<MultiDiracDeterminant>(*det));

  auto clone = std::make_unique<MultiSlaterDetTableMethod>(tqp, std::move(dets_clone), use_pre_computing_);

  clone->CI_Optimizable = CI_Optimizable;
  clone->C2node         = C2node;
  clone->C              = C;
  clone->myVars         = myVars;

  clone->Optimizable = Optimizable;
  clone->csf_data_   = csf_data_;

  return clone;
}

/** Compute VGL of this MultiSlaterDetTableMethod
 *
 * THis is introduced to remove redundant code in 
 * - evaluate(P,G,L)
 * - evaluateLog(P,G,L,buf,fillbuffer)
 * Miguel's note: can this change over time??? I don't know yet
 */
WaveFunctionComponent::LogValueType MultiSlaterDetTableMethod::evaluate_vgl_impl(const ParticleSet& P,
                                                                                 ParticleSet::ParticleGradient& g_tmp,
                                                                                 ParticleSet::ParticleLaplacian& l_tmp)
{
  const ValueType czero(0);
  psi_ratio_to_ref_det_ = czero;
  g_tmp                 = czero;
  l_tmp                 = czero;

  for (size_t ig = 0; ig < Dets.size(); ig++)
    precomputeC_otherDs(P, ig);

  for (size_t i = 0; i < Dets[0]->getNumDets(); ++i)
    psi_ratio_to_ref_det_ += C_otherDs[0][i] * Dets[0]->getRatiosToRefDet()[i];

  for (size_t id = 0; id < Dets.size(); id++)
    for (size_t i = 0; i < Dets[id]->getNumDets(); ++i)
      for (int k = 0, n = Dets[id]->getFirstIndex(); k < Dets[id]->getNumPtcls(); k++, n++)
      {
        g_tmp[n] += C_otherDs[id][i] * Dets[id]->getGrads()(i, k);
        l_tmp[n] += C_otherDs[id][i] * Dets[id]->getLapls()(i, k);
      }

  ValueType psiinv = static_cast<ValueType>(PsiValueType(1.0) / psi_ratio_to_ref_det_);
  g_tmp *= psiinv;
  l_tmp *= psiinv;
  LogValueType log_psi = convertValueToLog(psi_ratio_to_ref_det_);
  for (size_t id = 0; id < Dets.size(); id++)
    log_psi += Dets[id]->getLogValueRefDet();
  return log_psi;
}

WaveFunctionComponent::LogValueType MultiSlaterDetTableMethod::evaluateLog(const ParticleSet& P,
                                                                           ParticleSet::ParticleGradient& G,
                                                                           ParticleSet::ParticleLaplacian& L)
{
  ScopedTimer local_timer(EvaluateTimer);
  for (size_t id = 0; id < Dets.size(); id++)
  {
    if (P.isSpinor())
      Dets[id]->evaluateForWalkerMoveWithSpin(P);
    else
      Dets[id]->evaluateForWalkerMove(P);
  }

  log_value_ = evaluate_vgl_impl(P, myG, myL);

  G += myG;
  for (size_t i = 0; i < L.size(); i++)
    L[i] += myL[i] - dot(myG[i], myG[i]);

  return log_value_;
}

WaveFunctionComponent::PsiValueType MultiSlaterDetTableMethod::evalGrad_impl(ParticleSet& P,
                                                                             int iat,
                                                                             bool newpos,
                                                                             GradType& g_at)
{
  const int det_id = getDetID(iat);

  if (newpos)
    Dets[det_id]->evaluateDetsAndGradsForPtclMove(P, iat);
  else
    Dets[det_id]->evaluateGrads(P, iat);

  const auto& grads             = (newpos) ? Dets[det_id]->getNewGrads() : Dets[det_id]->getGrads();
  const ValueVector& detValues0 = (newpos) ? Dets[det_id]->getNewRatiosToRefDet() : Dets[det_id]->getRatiosToRefDet();
  const size_t noffset          = Dets[det_id]->getFirstIndex();

  PsiValueType psi(0);
  // enforce full precision reduction due to numerical sensitivity
  QTFull::GradType g_sum;
  for (size_t i = 0; i < Dets[det_id]->getNumDets(); i++)
  {
    psi += detValues0[i] * C_otherDs[det_id][i];
    g_sum += C_otherDs[det_id][i] * grads(i, iat - noffset);
  }

  g_at = g_sum / psi;
  return psi;
}

WaveFunctionComponent::PsiValueType MultiSlaterDetTableMethod::evalGradWithSpin_impl(ParticleSet& P,
                                                                                     int iat,
                                                                                     bool newpos,
                                                                                     GradType& g_at,
                                                                                     ComplexType& sg_at)
{
  const int det_id = getDetID(iat);

  if (newpos)
    Dets[det_id]->evaluateDetsAndGradsForPtclMoveWithSpin(P, iat);
  else
    Dets[det_id]->evaluateGradsWithSpin(P, iat);

  const auto& grads             = (newpos) ? Dets[det_id]->getNewGrads() : Dets[det_id]->getGrads();
  const ValueVector& detValues0 = (newpos) ? Dets[det_id]->getNewRatiosToRefDet() : Dets[det_id]->getRatiosToRefDet();
  const ValueMatrix& spingrads  = (newpos) ? Dets[det_id]->getNewSpinGrads() : Dets[det_id]->getSpinGrads();
  const size_t noffset          = Dets[det_id]->getFirstIndex();

  PsiValueType psi(0);
  for (size_t i = 0; i < Dets[det_id]->getNumDets(); i++)
  {
    psi += detValues0[i] * C_otherDs[det_id][i];
    g_at += C_otherDs[det_id][i] * grads(i, iat - noffset);
    sg_at += C_otherDs[det_id][i] * spingrads(i, iat - noffset);
  }
  g_at *= PsiValueType(1.0) / psi;
  sg_at *= PsiValueType(1.0) / psi;
  return psi;
}

void MultiSlaterDetTableMethod::mw_evalGrad_impl(const RefVectorWithLeader<WaveFunctionComponent>& WFC_list,
                                                 const RefVectorWithLeader<ParticleSet>& P_list,
                                                 int iat,
                                                 bool newpos,
                                                 std::vector<GradType>& grad_now,
                                                 std::vector<PsiValueType>& psi_list)
{
  auto& det_leader         = WFC_list.getCastedLeader<MultiSlaterDetTableMethod>();
  auto& det_value_ptr_list = det_leader.det_value_ptr_list;
  auto& C_otherDs_ptr_list = det_leader.C_otherDs_ptr_list;
  const int det_id         = det_leader.getDetID(iat);
  const int nw             = WFC_list.size();
  const int ndets          = det_leader.Dets[det_id]->getNumDets();


  RefVectorWithLeader<MultiDiracDeterminant> det_list(*det_leader.Dets[det_id]);
  det_list.reserve(WFC_list.size());
  ScopedTimer local_timer(det_leader.MWEvalGradTimer);
  for (int iw = 0; iw < WFC_list.size(); iw++)
  {
    auto& det = WFC_list.getCastedElement<MultiSlaterDetTableMethod>(iw);
    det_list.push_back(*det.Dets[det_id]);
  }

  if (newpos)
    det_leader.Dets[det_id]->mw_evaluateDetsAndGradsForPtclMove(det_list, P_list, iat);
  else
    det_leader.Dets[det_id]->mw_evaluateGrads(det_list, P_list, iat);

  det_value_ptr_list.resize(nw);
  C_otherDs_ptr_list.resize(nw);
  //Data layout change for grads function
  Matrix<ValueType, OMPallocator<ValueType, PinnedAlignedAllocator<ValueType>>> Grads_copy;
  Grads_copy.resize(3 * nw, ndets);

  for (size_t iw = 0; iw < nw; iw++)
  {
    auto& det = WFC_list.getCastedElement<MultiSlaterDetTableMethod>(iw);

    const size_t noffset = det.Dets[det_id]->getFirstIndex();
    const auto& grads    = (newpos) ? det.Dets[det_id]->getNewGrads() : det.Dets[det_id]->getGrads();

    for (size_t i = 0; i < ndets; i++)
    {
      Grads_copy[3 * iw + 0][i] = grads(i, iat - noffset)[0];
      Grads_copy[3 * iw + 1][i] = grads(i, iat - noffset)[1];
      Grads_copy[3 * iw + 2][i] = grads(i, iat - noffset)[2];
    }

    const ValueType* restrict detValues0 =
        (newpos) ? det.Dets[det_id]->getNewRatiosToRefDet().data() : det.Dets[det_id]->getRatiosToRefDet().data();
    // allocate device memory and transfer content to device
    PRAGMA_OFFLOAD("omp target enter data map(to : detValues0[:ndets])")
    det_value_ptr_list[iw] = getOffloadDevicePtr(detValues0);
    //Put C_otherDs in host
    C_otherDs_ptr_list[iw] = det.C_otherDs[det_id].device_data();
  }

  std::vector<ValueType> grad_now_list(nw * 3, 0);
  auto* grad_now_list_ptr      = grad_now_list.data();
  auto* Grads_copy_ptr         = Grads_copy.data();
  auto* psi_list_ptr           = psi_list.data();
  auto* C_otherDs_ptr_list_ptr = C_otherDs_ptr_list.data();
  auto* det_value_ptr_list_ptr = det_value_ptr_list.data();
  {
    ScopedTimer local_timer(det_leader.OffloadGradTimer);
    PRAGMA_OFFLOAD("omp target teams distribute map(from: psi_list_ptr[:nw]) \
                    map(from: grad_now_list_ptr[:3 * nw]) \
                    map(always, to: det_value_ptr_list_ptr[:nw], C_otherDs_ptr_list_ptr[:nw]) \
                    map(always, to: Grads_copy_ptr[:Grads_copy.size()])")
    for (size_t iw = 0; iw < nw; iw++)
    {
      // enforce full precision reduction due to numerical sensitivity
      PsiValueType psi_local(0);
      PsiValueType grad_local_x(0);
      PsiValueType grad_local_y(0);
      PsiValueType grad_local_z(0);
      PRAGMA_OFFLOAD("omp parallel for reduction(+:psi_local, grad_local_x, grad_local_y, grad_local_z)")
      for (size_t i = 0; i < ndets; i++)
      {
        psi_local += det_value_ptr_list_ptr[iw][i] * C_otherDs_ptr_list_ptr[iw][i];
        grad_local_x += C_otherDs_ptr_list_ptr[iw][i] * Grads_copy_ptr[(3 * iw + 0) * ndets + i];
        grad_local_y += C_otherDs_ptr_list_ptr[iw][i] * Grads_copy_ptr[(3 * iw + 1) * ndets + i];
        grad_local_z += C_otherDs_ptr_list_ptr[iw][i] * Grads_copy_ptr[(3 * iw + 2) * ndets + i];
      }
      psi_list_ptr[iw]              = psi_local;
      grad_now_list_ptr[iw * 3 + 0] = grad_local_x;
      grad_now_list_ptr[iw * 3 + 1] = grad_local_y;
      grad_now_list_ptr[iw * 3 + 2] = grad_local_z;
    }
  }

  for (size_t iw = 0; iw < nw; iw++)
  {
    auto psi_inv    = static_cast<ValueType>(PsiValueType(1.0) / psi_list[iw]);
    grad_now[iw][0] = grad_now_list[iw * 3 + 0] * psi_inv;
    grad_now[iw][1] = grad_now_list[iw * 3 + 1] * psi_inv;
    grad_now[iw][2] = grad_now_list[iw * 3 + 2] * psi_inv;

    //Free Memory
    auto& det = WFC_list.getCastedElement<MultiSlaterDetTableMethod>(iw);
    const ValueType* restrict detValues0 =
        (newpos) ? det.Dets[det_id]->getNewRatiosToRefDet().data() : det.Dets[det_id]->getRatiosToRefDet().data();
    PRAGMA_OFFLOAD("omp target exit data map(delete : detValues0[:ndets])") //free memory on device
  }
}

WaveFunctionComponent::PsiValueType MultiSlaterDetTableMethod::evalGrad_impl_no_precompute(ParticleSet& P,
                                                                                           int iat,
                                                                                           bool newpos,
                                                                                           GradType& g_at)
{
  const int det_id = getDetID(iat);

  if (newpos)
    Dets[det_id]->evaluateDetsAndGradsForPtclMove(P, iat);
  else
    Dets[det_id]->evaluateGrads(P, iat);

  const auto& grads              = (newpos) ? Dets[det_id]->getNewGrads() : Dets[det_id]->getGrads();
  const auto& detValues0         = (newpos) ? Dets[det_id]->getNewRatiosToRefDet() : Dets[det_id]->getRatiosToRefDet();
  const size_t* restrict det0    = (*C2node)[det_id].data();
  const ValueType* restrict cptr = C->data();
  const size_t nc                = C->size();
  const size_t noffset           = Dets[det_id]->getFirstIndex();
  PsiValueType psi(0);
  for (size_t i = 0; i < nc; ++i)
  {
    const size_t d0 = det0[i];
    ValueType t     = cptr[i];
    for (size_t id = 0; id < Dets.size(); id++)
      if (id != det_id)
        t *= Dets[id]->getRatiosToRefDet()[(*C2node)[id][i]];
    psi += t * detValues0[d0];
    g_at += t * grads(d0, iat - noffset);
  }
  g_at *= PsiValueType(1.0) / psi;
  return psi;
}

WaveFunctionComponent::PsiValueType MultiSlaterDetTableMethod::evalGradWithSpin_impl_no_precompute(ParticleSet& P,
                                                                                                   int iat,
                                                                                                   bool newpos,
                                                                                                   GradType& g_at,
                                                                                                   ComplexType& sg_at)
{
  const int det_id = getDetID(iat);

  if (newpos)
    Dets[det_id]->evaluateDetsAndGradsForPtclMoveWithSpin(P, iat);
  else
    Dets[det_id]->evaluateGradsWithSpin(P, iat);

  const auto& grads              = (newpos) ? Dets[det_id]->getNewGrads() : Dets[det_id]->getGrads();
  const auto& detValues0         = (newpos) ? Dets[det_id]->getNewRatiosToRefDet() : Dets[det_id]->getRatiosToRefDet();
  const auto& spingrads          = (newpos) ? Dets[det_id]->getNewSpinGrads() : Dets[det_id]->getSpinGrads();
  const size_t* restrict det0    = (*C2node)[det_id].data();
  const ValueType* restrict cptr = C->data();
  const size_t nc                = C->size();
  const size_t noffset           = Dets[det_id]->getFirstIndex();
  PsiValueType psi(0);
  for (size_t i = 0; i < nc; ++i)
  {
    const size_t d0 = det0[i];
    ValueType t     = cptr[i];
    for (size_t id = 0; id < Dets.size(); id++)
      if (id != det_id)
        t *= Dets[id]->getRatiosToRefDet()[(*C2node)[id][i]];
    psi += t * detValues0[d0];
    g_at += t * grads(d0, iat - noffset);
    sg_at += t * spingrads(d0, iat - noffset);
  }
  g_at *= PsiValueType(1.0) / psi;
  sg_at *= PsiValueType(1.0) / psi;
  return psi;
}

WaveFunctionComponent::GradType MultiSlaterDetTableMethod::evalGrad(ParticleSet& P, int iat)
{
  ScopedTimer local_timer(EvalGradTimer);

  GradType grad_iat;
  if (use_pre_computing_)
    evalGrad_impl(P, iat, false, grad_iat);
  else
    evalGrad_impl_no_precompute(P, iat, false, grad_iat);

  return grad_iat;
}

WaveFunctionComponent::GradType MultiSlaterDetTableMethod::evalGradWithSpin(ParticleSet& P,
                                                                            int iat,
                                                                            ComplexType& spingrad)
{
  ScopedTimer local_timer(EvalGradTimer);

  GradType grad_iat;
  ComplexType spingrad_iat;
  if (use_pre_computing_)
    evalGradWithSpin_impl(P, iat, false, grad_iat, spingrad_iat);
  else
    evalGradWithSpin_impl_no_precompute(P, iat, false, grad_iat, spingrad_iat);

  return grad_iat;
}

void MultiSlaterDetTableMethod::mw_evalGrad(const RefVectorWithLeader<WaveFunctionComponent>& WFC_list,
                                            const RefVectorWithLeader<ParticleSet>& P_list,
                                            int iat,
                                            std::vector<GradType>& grad_now) const
{
  if (!use_pre_computing_)
  {
    WaveFunctionComponent::mw_evalGrad(WFC_list, P_list, iat, grad_now);
    return;
  }

  const int nw = WFC_list.size();

  std::vector<PsiValueType> psi_list(nw, 0);
  mw_evalGrad_impl(WFC_list, P_list, iat, false, grad_now, psi_list);
}


WaveFunctionComponent::PsiValueType MultiSlaterDetTableMethod::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  ScopedTimer local_timer(RatioGradTimer);
  UpdateMode = ORB_PBYP_PARTIAL;

  GradType dummy;
  if (use_pre_computing_)
    new_psi_ratio_to_new_ref_det_ = evalGrad_impl(P, iat, true, dummy);
  else
    new_psi_ratio_to_new_ref_det_ = evalGrad_impl_no_precompute(P, iat, true, dummy);

  const int det_id = getDetID(iat);
  curRatio         = Dets[det_id]->getRefDetRatio() * new_psi_ratio_to_new_ref_det_ / psi_ratio_to_ref_det_;
  grad_iat += dummy;
  return curRatio;
}

WaveFunctionComponent::PsiValueType MultiSlaterDetTableMethod::ratioGradWithSpin(ParticleSet& P,
                                                                                 int iat,
                                                                                 GradType& grad_iat,
                                                                                 ComplexType& spingrad_iat)
{
  ScopedTimer local_timer(RatioGradTimer);
  UpdateMode = ORB_PBYP_PARTIAL;

  GradType dummy;
  ComplexType spindummy;
  if (use_pre_computing_)
    new_psi_ratio_to_new_ref_det_ = evalGradWithSpin_impl(P, iat, true, dummy, spindummy);
  else
    new_psi_ratio_to_new_ref_det_ = evalGradWithSpin_impl_no_precompute(P, iat, true, dummy, spindummy);

  const int det_id = getDetID(iat);
  curRatio         = Dets[det_id]->getRefDetRatio() * new_psi_ratio_to_new_ref_det_ / psi_ratio_to_ref_det_;
  grad_iat += dummy;
  spingrad_iat += spindummy;
  return curRatio;
}

void MultiSlaterDetTableMethod::mw_ratioGrad(const RefVectorWithLeader<WaveFunctionComponent>& WFC_list,
                                             const RefVectorWithLeader<ParticleSet>& P_list,
                                             int iat,
                                             std::vector<WaveFunctionComponent::PsiValueType>& ratios,
                                             std::vector<GradType>& grad_new) const
{
  if (!use_pre_computing_)
  {
    WaveFunctionComponent::mw_ratioGrad(WFC_list, P_list, iat, ratios, grad_new);
    return;
  }

  auto& det_leader         = WFC_list.getCastedLeader<MultiSlaterDetTableMethod>();
  auto& det_value_ptr_list = det_leader.det_value_ptr_list;
  auto& C_otherDs_ptr_list = det_leader.C_otherDs_ptr_list;
  const int nw             = WFC_list.size();

  ScopedTimer local_timer(det_leader.MWRatioGradTimer);
  std::vector<PsiValueType> psi_list(nw, 0);
  std::vector<GradType> dummy;
  dummy.resize(nw);

  mw_evalGrad_impl(WFC_list, P_list, iat, true, dummy, psi_list);

  const int det_id = getDetID(iat);
  for (size_t iw = 0; iw < nw; iw++)
  {
    auto& det                         = WFC_list.getCastedElement<MultiSlaterDetTableMethod>(iw);
    det.new_psi_ratio_to_new_ref_det_ = psi_list[iw];
    grad_new[iw] += dummy[iw];
    ratios[iw] = det.curRatio = det.Dets[det_id]->getRefDetRatio() * psi_list[iw] / det.psi_ratio_to_ref_det_;
  }
}

WaveFunctionComponent::PsiValueType MultiSlaterDetTableMethod::ratio_impl(ParticleSet& P, int iat)
{
  const int det_id = getDetID(iat);

  Dets[det_id]->evaluateDetsForPtclMove(P, iat);

  const ValueVector& detValues0 = Dets[det_id]->getNewRatiosToRefDet();

  PsiValueType psi = 0;
  // This function computes
  // psi=Det_Coeff[i]*Det_Value[unique_det_up]*Det_Value[unique_det_dn]*Det_Value[unique_det_AnyOtherType]
  // Since only one electron group is moved at the time, identified by det_id, We precompute:
  // C_otherDs[det_id][i]=Det_Coeff[i]*Det_Value[unique_det_dn]*Det_Value[unique_det_AnyOtherType]
  for (size_t i = 0; i < Dets[det_id]->getNumDets(); i++)
    psi += detValues0[i] * C_otherDs[det_id][i];

  return psi;
}


WaveFunctionComponent::PsiValueType MultiSlaterDetTableMethod::ratio_impl_no_precompute(ParticleSet& P, int iat)
{
  const int det_id = getDetID(iat);
  Dets[det_id]->evaluateDetsForPtclMove(P, iat);

  const ValueVector& detValues0  = Dets[det_id]->getNewRatiosToRefDet(); //always new
  const size_t* restrict det0    = (*C2node)[det_id].data();
  const ValueType* restrict cptr = C->data();
  const size_t nc                = C->size();

  PsiValueType psi = 0;
  for (size_t i = 0; i < nc; ++i)
  {
    ValueType t = cptr[i];
    for (size_t id = 0; id < Dets.size(); id++)
      if (id != det_id)
        t *= Dets[id]->getRatiosToRefDet()[(*C2node)[id][i]];
    t *= detValues0[det0[i]];
    psi += t;
  }
  return psi;
}

// use ci_node for this routine only
WaveFunctionComponent::PsiValueType MultiSlaterDetTableMethod::ratio(ParticleSet& P, int iat)
{
  ScopedTimer local_timer(RatioTimer);
  UpdateMode = ORB_PBYP_RATIO;

  if (use_pre_computing_)
    new_psi_ratio_to_new_ref_det_ = ratio_impl(P, iat);
  else
    new_psi_ratio_to_new_ref_det_ = ratio_impl_no_precompute(P, iat);

  const int det_id = getDetID(iat);
  curRatio         = Dets[det_id]->getRefDetRatio() * new_psi_ratio_to_new_ref_det_ / psi_ratio_to_ref_det_;
  return curRatio;
}

void MultiSlaterDetTableMethod::mw_calcRatio(const RefVectorWithLeader<WaveFunctionComponent>& WFC_list,
                                             const RefVectorWithLeader<ParticleSet>& P_list,
                                             int iat,
                                             std::vector<PsiValueType>& ratios) const
{
  if (!use_pre_computing_)
  {
    WaveFunctionComponent::mw_calcRatio(WFC_list, P_list, iat, ratios);
    return;
  }

  const int det_id = getDetID(iat);

  const int nw             = WFC_list.size();
  auto& det_leader         = WFC_list.getCastedLeader<MultiSlaterDetTableMethod>();
  auto& det_value_ptr_list = det_leader.det_value_ptr_list;
  auto& C_otherDs_ptr_list = det_leader.C_otherDs_ptr_list;
  const int ndets          = det_leader.Dets[det_id]->getNumDets();

  ScopedTimer local_timer(det_leader.MWRatioTimer);

  RefVectorWithLeader<MultiDiracDeterminant> det_list(*det_leader.Dets[det_id]);
  det_list.reserve(WFC_list.size());
  for (int iw = 0; iw < WFC_list.size(); iw++)
  {
    auto& det = WFC_list.getCastedElement<MultiSlaterDetTableMethod>(iw);
    det_list.push_back(*det.Dets[det_id]);
  }

  det_leader.Dets[det_id]->mw_evaluateDetsForPtclMove(det_list, P_list, iat);

  det_value_ptr_list.resize(nw);
  C_otherDs_ptr_list.resize(nw);

  for (size_t iw = 0; iw < nw; iw++)
  {
    auto& det      = WFC_list.getCastedElement<MultiSlaterDetTableMethod>(iw);
    det.UpdateMode = ORB_PBYP_RATIO;

    const ValueType* restrict detValues0 = det.Dets[det_id]->getNewRatiosToRefDet().data(); //always new
    // allocate device memory and transfer content to device
    PRAGMA_OFFLOAD("omp target enter data map(to : detValues0[:ndets])")
    det_value_ptr_list[iw] = getOffloadDevicePtr(detValues0);
    C_otherDs_ptr_list[iw] = det.C_otherDs[det_id].device_data();
  }

  std::vector<PsiValueType> psi_list(nw, 0);
  auto* psi_list_ptr           = psi_list.data();
  auto* C_otherDs_ptr_list_ptr = C_otherDs_ptr_list.data();
  auto* det_value_ptr_list_ptr = det_value_ptr_list.data();
  OffloadRatioTimer.start();
  PRAGMA_OFFLOAD("omp target teams distribute map(from: psi_list_ptr[:nw]) \
          map(always, to: det_value_ptr_list_ptr[:nw], C_otherDs_ptr_list_ptr[:nw])")
  for (size_t iw = 0; iw < nw; iw++)
  {
    PsiValueType psi_local(0);
    PRAGMA_OFFLOAD("omp parallel for reduction(+ : psi_local)")
    for (size_t i = 0; i < ndets; i++)
      psi_local += det_value_ptr_list_ptr[iw][i] * C_otherDs_ptr_list_ptr[iw][i];
    psi_list_ptr[iw] = psi_local;
  }
  OffloadRatioTimer.stop();
  for (size_t iw = 0; iw < nw; iw++)
  {
    auto& det                         = WFC_list.getCastedElement<MultiSlaterDetTableMethod>(iw);
    det.new_psi_ratio_to_new_ref_det_ = psi_list[iw];
    ratios[iw] = det.curRatio = det.Dets[det_id]->getRefDetRatio() * psi_list[iw] / det.psi_ratio_to_ref_det_;

    const ValueType* restrict detValues0 = det.Dets[det_id]->getNewRatiosToRefDet().data(); //always new
    PRAGMA_OFFLOAD("omp target exit data map(delete : detValues0[:ndets])")                 //free memory on device.
  }
}

void MultiSlaterDetTableMethod::evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios)
{
  ScopedTimer local_timer(RatioTimer);

  const int det_id = getDetID(VP.refPtcl);

  for (size_t iat = 0; iat < VP.getTotalNum(); ++iat)
  {
    Dets[det_id]->evaluateDetsForPtclMove(VP, iat, VP.refPtcl);
    const ValueVector& detValues0 = Dets[det_id]->getNewRatiosToRefDet();

    PsiValueType psiNew(0);
    if (use_pre_computing_)
      for (size_t i = 0; i < Dets[det_id]->getNumDets(); i++)
        psiNew += detValues0[i] * C_otherDs[det_id][i];
    else
    {
      const size_t* restrict det0    = (*C2node)[det_id].data();
      const ValueType* restrict cptr = C->data();
      const size_t nc                = C->size();

      for (size_t i = 0; i < nc; ++i)
      {
        ValueType t = cptr[i];
        for (size_t id = 0; id < Dets.size(); id++)
          if (id != det_id)
            t *= Dets[id]->getRatiosToRefDet()[(*C2node)[id][i]];
        t *= detValues0[det0[i]];
        psiNew += t;
      }
    }
    ratios[iat] = Dets[det_id]->getRefDetRatio() * psiNew / psi_ratio_to_ref_det_;
  }
}

void MultiSlaterDetTableMethod::acceptMove(ParticleSet& P, int iat, bool safe_to_delay)
{
  // this should depend on the type of update, ratio / ratioGrad
  // for now is incorrect fot ratio(P,iat,dG,dL) updates
  ScopedTimer local_timer(AccRejTimer);
  // update psi_ratio_to_ref_det_,myG_temp,myL_temp
  psi_ratio_to_ref_det_ = new_psi_ratio_to_new_ref_det_;
  log_value_ += convertValueToLog(curRatio);
  curRatio = 1.0;

  Dets[getDetID(iat)]->acceptMove(P, iat, safe_to_delay);
}

void MultiSlaterDetTableMethod::restore(int iat)
{
  ScopedTimer local_timer(AccRejTimer);
  Dets[getDetID(iat)]->restore(iat);
  curRatio = 1.0;
}

void MultiSlaterDetTableMethod::registerData(ParticleSet& P, WFBufferType& buf)
{
  for (size_t id = 0; id < Dets.size(); id++)
    Dets[id]->registerData(P, buf);

  buf.add(log_value_);
  buf.add(psi_ratio_to_ref_det_);
}

WaveFunctionComponent::LogValueType MultiSlaterDetTableMethod::updateBuffer(ParticleSet& P,
                                                                            WFBufferType& buf,
                                                                            bool fromscratch)
{
  ScopedTimer local_timer(UpdateTimer);

  for (size_t id = 0; id < Dets.size(); id++)
    Dets[id]->updateBuffer(P, buf, fromscratch);

  log_value_ = evaluate_vgl_impl(P, myG, myL);

  P.G += myG;
  for (int i = 0; i < P.L.size(); i++)
    P.L[i] += myL[i] - dot(myG[i], myG[i]);

  buf.put(log_value_);
  buf.put(psi_ratio_to_ref_det_);

  return log_value_;
}

void MultiSlaterDetTableMethod::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  for (size_t id = 0; id < Dets.size(); id++)
    Dets[id]->copyFromBuffer(P, buf);

  buf.get(log_value_);
  buf.get(psi_ratio_to_ref_det_);
}


void MultiSlaterDetTableMethod::checkInVariables(opt_variables_type& active)
{
  if (CI_Optimizable)
  {
    if (myVars->size())
      active.insertFrom(*myVars);
    else
      Optimizable = false;
  }
  bool all_Optimizable = true;
  for (size_t id = 0; id < Dets.size() && all_Optimizable; id++)
    all_Optimizable = Dets[id]->Optimizable;

  if (all_Optimizable)
    for (size_t id = 0; id < Dets.size(); id++)
      Dets[id]->checkInVariables(active);
}

void MultiSlaterDetTableMethod::checkOutVariables(const opt_variables_type& active)
{
  if (CI_Optimizable)
    myVars->getIndex(active);

  bool all_Optimizable = true;
  for (size_t id = 0; id < Dets.size() && all_Optimizable; id++)
    all_Optimizable = Dets[id]->Optimizable;

  if (all_Optimizable)
    for (size_t id = 0; id < Dets.size(); id++)
      Dets[id]->checkOutVariables(active);
}

/** resetParameters with optVariables
 *
 * USE_resetParameters
 */
void MultiSlaterDetTableMethod::resetParameters(const opt_variables_type& active)
{
  if (CI_Optimizable)
  {
    if (csf_data_)
    {
      ValueType* restrict CSFcoeff_p = csf_data_->coeffs.data();
      for (int i = 0; i < csf_data_->coeffs.size() - 1; i++)
      {
        int loc = myVars->where(i);
        if (loc >= 0)
        {
          CSFcoeff_p[i + 1] = (*myVars)[i] = active[loc];
        }
      }
      int cnt                                 = 0;
      ValueType* restrict C_p                 = C->data();
      const RealType* restrict CSFexpansion_p = csf_data_->expansion.data();
      for (int i = 0; i < csf_data_->dets_per_csf.size(); i++)
      {
        for (int k = 0; k < csf_data_->dets_per_csf[i]; k++)
        {
          C_p[cnt] = CSFcoeff_p[i] * CSFexpansion_p[cnt];
          cnt++;
        }
      }
      //for(int i=0; i<Dets.size(); i++) Dets[i]->resetParameters(active);
    }
    else
    {
      ValueType* restrict C_p = C->data();
      for (int i = 0; i < C->size() - 1; i++)
      {
        int loc = myVars->where(i);
        if (loc >= 0)
        {
          C_p[i + 1] = (*myVars)[i] = active[loc];
        }
      }
      //for(int i=0; i<Dets.size(); i++) Dets[i]->resetParameters(active);
    }
  }
  bool all_Optimizable = true;
  for (size_t id = 0; id < Dets.size() && all_Optimizable; id++)
    all_Optimizable = Dets[id]->Optimizable;

  if (all_Optimizable)
    for (size_t id = 0; id < Dets.size(); id++)
      Dets[id]->resetParameters(active);
}
void MultiSlaterDetTableMethod::reportStatus(std::ostream& os) {}


void MultiSlaterDetTableMethod::evaluateDerivatives(ParticleSet& P,
                                                    const opt_variables_type& optvars,
                                                    std::vector<ValueType>& dlogpsi,
                                                    std::vector<ValueType>& dhpsioverpsi)
{
  evaluateDerivativesWF(P, optvars, dlogpsi);
  if (CI_Optimizable)
  {
    bool recalculate(false);
    for (int k = 0; k < myVars->size(); ++k)
    {
      int kk = myVars->where(k);
      if (kk < 0)
        continue;
      if (optvars.recompute(kk))
        recalculate = true;
    }
    // need to modify for CSF later on, right now assume Slater Det basis
    if (recalculate)
    {
      ValueType psiinv = static_cast<ValueType>(PsiValueType(1.0) / psi_ratio_to_ref_det_);
      laplSum.resize(Dets.size());
      for (size_t id = 0; id < Dets.size(); id++)
      {
        laplSum[id].resize(Dets[id]->getNumDets());
        // assume that evaluateLog has been called in opt routine before
        //   Dets[id]->evaluateForWalkerMove(P);
        // myG,myL should already be calculated
        for (size_t i = 0; i < laplSum[id].size(); i++)
        {
          laplSum[id][i] = 0.0;
          for (size_t k = 0; k < Dets[id]->getNumPtcls(); k++)
            laplSum[id][i] += Dets[id]->getLapls()[i][k];
        }
      }

      ValueType lapl_sum = 0.0;
      myG_temp           = 0.0;
      for (size_t id = 0; id < Dets.size(); id++)
        for (size_t i = 0; i < Dets[id]->getNumDets(); i++)
        {
          // assume C_otherDs prepared by evaluateLog already
          ValueType tmp = C_otherDs[id][i] * psiinv;
          lapl_sum += tmp * laplSum[id][i];
          for (size_t k = 0, j = Dets[id]->getFirstIndex(); k < Dets[id]->getNumPtcls(); k++, j++)
            myG_temp[j] += tmp * Dets[id]->getGrads()(i, k);
        }

      ValueType gg = 0.0;
      for (size_t i = 0; i < P.getTotalNum(); i++)
        gg += dot(myG_temp[i], myG_temp[i]) - dot(P.G[i], myG_temp[i]);

      if (csf_data_)
      {
        const int num = csf_data_->coeffs.size() - 1;
        int cnt       = 0;
        //        this one is not optable
        cnt += csf_data_->dets_per_csf[0];
        int ip(1);
        for (int i = 0; i < num; i++, ip++)
        {
          int kk = myVars->where(i);
          if (kk < 0)
          {
            cnt += csf_data_->dets_per_csf[ip];
            continue;
          }
          ValueType q0 = 0.0;
          std::vector<ValueType> v(Dets.size());
          const RealType* restrict CSFexpansion_p = csf_data_->expansion.data();
          for (int k = 0; k < csf_data_->dets_per_csf[ip]; k++)
          {
            for (size_t id = 0; id < Dets.size(); id++)
            {
              const auto& grads_spin = Dets[id]->getGrads();
              size_t spinC           = (*C2node)[id][cnt];
              ValueType tmp          = CSFexpansion_p[cnt] * psiinv;
              for (size_t other_id = 0; other_id < Dets.size(); other_id++)
              {
                if (id == other_id)
                  continue;
                const ValueVector& detValues_otherspin = Dets[other_id]->getRatiosToRefDet();
                size_t otherspinC                      = (*C2node)[other_id][cnt];
                tmp *= detValues_otherspin[otherspinC];
              }
              q0 += tmp * laplSum[id][spinC];
              for (size_t l = 0, j = Dets[id]->getFirstIndex(); l < Dets[id]->getNumPtcls(); l++, j++)
                v[id] += tmp *
                    static_cast<ValueType>(dot(P.G[j], grads_spin(spinC, l)) - dot(myG_temp[j], grads_spin(spinC, l)));
            }
            cnt++;
          }
          ValueType dhpsi = (RealType)-0.5 * (q0 - dlogpsi[kk] * lapl_sum) - dlogpsi[kk] * gg;
          for (size_t id = 0; id < Dets.size(); id++)
            dhpsi -= v[id];
          dhpsioverpsi[kk] = dhpsi;
        }
      }
      else
      { //usingDETS
        for (size_t i = 1; i < C->size(); i++)
        {
          int kk = myVars->where(i - 1);
          if (kk < 0)
            continue;

          ValueType q0 = 0.0;
          std::vector<ValueType> v(Dets.size());
          for (size_t id = 0; id < Dets.size(); id++)
          {
            const auto& grads_spin = Dets[id]->getGrads();
            size_t spinC           = (*C2node)[id][i];
            ValueType tmp          = psiinv;
            for (size_t other_id = 0; other_id < Dets.size(); other_id++)
            {
              if (id == other_id)
                continue;
              size_t otherspinC = (*C2node)[other_id][i];
              tmp *= Dets[other_id]->getRatiosToRefDet()[otherspinC];
            }
            q0 += tmp * laplSum[id][spinC];
            for (size_t l = 0, j = Dets[id]->getFirstIndex(); l < Dets[id]->getNumPtcls(); l++, j++)
              v[id] += tmp *
                  static_cast<ValueType>(dot(P.G[j], grads_spin(spinC, l)) - dot(myG_temp[j], grads_spin(spinC, l)));
          }
          ValueType dhpsi = (RealType)-0.5 * (q0 - dlogpsi[kk] * lapl_sum) - dlogpsi[kk] * gg;
          for (size_t id = 0; id < Dets.size(); id++)
            dhpsi -= v[id];
          dhpsioverpsi[kk] = dhpsi;
        }
      }
    }
  }

  evaluateMultiDiracDeterminantDerivatives(P, optvars, dlogpsi, dhpsioverpsi);
}

void MultiSlaterDetTableMethod::evaluateMultiDiracDeterminantDerivatives(ParticleSet& P,
                                                                         const opt_variables_type& optvars,
                                                                         std::vector<ValueType>& dlogpsi,
                                                                         std::vector<ValueType>& dhpsioverpsi)
{
  //Currently, the MultiDiracDeterminant::evaluateDerivatives works with a legacy design, essentially requiring only up and down determinants.
  //e.g. for spinor cases, we only have one determinant so this interface doesn't work.
  //Here we throw an error only if the optimization is turned on for MultiDiracDeterminants until the code is updated
  for (auto const& det : Dets)
    if (!det->Optimizable)
      return;

  if (Dets.size() != 2)
  {
    throw std::runtime_error(
        "MultiSlaterDetTableMethod::evaluateDerivatives only compatible with two quantum particle types.");
  }
  else
  {
    Dets[0]->evaluateDerivatives(P, optvars, dlogpsi, dhpsioverpsi, *Dets[1],
                                 static_cast<ValueType>(psi_ratio_to_ref_det_), *C, (*C2node)[0], (*C2node)[1]);
    Dets[1]->evaluateDerivatives(P, optvars, dlogpsi, dhpsioverpsi, *Dets[0],
                                 static_cast<ValueType>(psi_ratio_to_ref_det_), *C, (*C2node)[1], (*C2node)[0]);
  }

  //note: the future redesign of MultiDiracDeterminant::evaluateDerivatives should look something like this
  //for (size_t id = 0; id < Dets.size(); id++)
  //  Dets[id]->evaluateDerivatives(P, optvars, dlogpsi, dhpsioverpsi, *Dets, static_cast<ValueType>(psi_ratio_to_ref_det_), *C, *C2node, id);
}

void MultiSlaterDetTableMethod::evaluateDerivativesWF(ParticleSet& P,
                                                      const opt_variables_type& optvars,
                                                      std::vector<ValueType>& dlogpsi)
{
  if (CI_Optimizable)
  {
    bool recalculate(false);
    for (int k = 0; k < myVars->size(); ++k)
    {
      int kk = myVars->where(k);
      if (kk < 0)
        continue;
      if (optvars.recompute(kk))
        recalculate = true;
    }
    // need to modify for CSF later on, right now assume Slater Det basis
    if (recalculate)
    {
      if (csf_data_)
      {
        ValueType psiinv = static_cast<ValueType>(PsiValueType(1.0) / psi_ratio_to_ref_det_);

        const int num = csf_data_->coeffs.size() - 1;
        int cnt       = 0;
        //        this one is not optable
        cnt += csf_data_->dets_per_csf[0];
        int ip(1);
        for (int i = 0; i < num; i++, ip++)
        {
          int kk = myVars->where(i);
          if (kk < 0)
          {
            cnt += csf_data_->dets_per_csf[ip];
            continue;
          }
          ValueType cdet                          = 0.0;
          const RealType* restrict CSFexpansion_p = csf_data_->expansion.data();
          for (int k = 0; k < csf_data_->dets_per_csf[ip]; k++)
          {
            ValueType t = CSFexpansion_p[cnt] * psiinv;
            // assume that evaluateLog has been called in opt routine before
            for (size_t id = 0; id < Dets.size(); id++)
              t *= Dets[id]->getRatiosToRefDet()[(*C2node)[id][cnt]];
            cdet += t;
            cnt++;
          }
          dlogpsi[kk] = cdet;
        }
      }
      else
      //usingDETS
      {
        ValueType psiinv = static_cast<ValueType>(PsiValueType(1.0) / psi_ratio_to_ref_det_);
        for (size_t i = 1; i < C->size(); i++)
        {
          int kk = myVars->where(i - 1);
          if (kk < 0)
            continue;
          ValueType cdet = psiinv;
          // assume that evaluateLog has been called in opt routine before
          for (size_t id = 0; id < Dets.size(); id++)
            cdet *= Dets[id]->getRatiosToRefDet()[(*C2node)[id][i]];
          dlogpsi[kk] = cdet;
        }
      }
    }
  }

  evaluateMultiDiracDeterminantDerivativesWF(P, optvars, dlogpsi);
}

void MultiSlaterDetTableMethod::evaluateMultiDiracDeterminantDerivativesWF(ParticleSet& P,
                                                                           const opt_variables_type& optvars,
                                                                           std::vector<ValueType>& dlogpsi)
{
  //Currently, the MultiDiracDeterminant::evaluateDerivativesWF works with a legacy design, essentially requiring only up and down determinants.
  //e.g. for spinor cases, we only have one determinant so this interface doesn't work.
  //Here we throw an error only if the optimization is turned on for MultiDiracDeterminants until the code is updated
  for (auto const& det : Dets)
    if (!det->Optimizable)
      return;

  if (Dets.size() != 2)
  {
    throw std::runtime_error(
        "MultiSlaterDetTableMethod::evaluateDerivativesWF only compatible with two quantum particle types.");
  }
  else
  {
    // FIXME this needs to be fixed by SPF to separate evaluateDerivatives and evaluateDerivativesWF for orbital rotation matrix
    Dets[0]->evaluateDerivativesWF(P, optvars, dlogpsi, *Dets[1], psi_ratio_to_ref_det_, *C, (*C2node)[0],
                                   (*C2node)[1]);
    Dets[1]->evaluateDerivativesWF(P, optvars, dlogpsi, *Dets[0], psi_ratio_to_ref_det_, *C, (*C2node)[1],
                                   (*C2node)[0]);
  }
  //note: the future redesign of MultiDiracDeterminant::evaluateDerivativesWF should look something like this
  // for (size_t id = 0; id < Dets.size(); id++)
  //   Dets[id]->evaluateDerivativesWF(P, optvars, dlogpsi, *Dets, psi_ratio_to_ref_det_, *C, *C2node, id);
}

void MultiSlaterDetTableMethod::buildOptVariables()
{
  for (size_t id = 0; id < Dets.size(); id++)
    Dets[id]->buildOptVariables((*C2node)[id]);
}

void MultiSlaterDetTableMethod::registerTimers()
{
  RatioTimer.reset();
  EvalGradTimer.reset();
  MWEvalGradTimer.reset();
  RatioGradTimer.reset();
  PrepareGroupTimer.reset();
  UpdateTimer.reset();
  EvaluateTimer.reset();
  AccRejTimer.reset();
}

void MultiSlaterDetTableMethod::prepareGroup(ParticleSet& P, int ig)
{
  if (!use_pre_computing_)
    return;
  precomputeC_otherDs(P, ig);
}

void MultiSlaterDetTableMethod::precomputeC_otherDs(const ParticleSet& P, int ig)
{
  // This function computes
  // C_otherDs[det_id][i]=Det_Coeff[i]*Det_Value[unique_det_dn]*Det_Value[unique_det_AnyOtherType]
  // Since only one electron group is moved at the time, identified by det_id, We precompute C_otherDs[det_id][i]:
  // psi=Det_Coeff[i]*Det_Value[unique_det_up]*Det_Value[unique_det_dn]*Det_Value[unique_det_AnyOtherType]
  // becomes:
  // psi=Det_Value[unique_det_up]*C_otherDs[det_id][i]
  // ig is the id of the group electron being moved. In this function, we compute the other groups
  // of electrons.
  // We loop over the number of type of determinants (up, diwn, positrons, etc), but we only multiply for ll types BUT ig
  // C_otherDs(0, :) stores C x D_dn x D_pos
  // C_otherDs(1, :) stores C x D_up x D_pos
  // C_otherDs(2, :) stores C x D_up x D_dn

  ScopedTimer local_timer(PrepareGroupTimer);
  C_otherDs[ig].resize(Dets[ig]->getNumDets());
  std::fill(C_otherDs[ig].begin(), C_otherDs[ig].end(), ValueType(0));
  for (size_t i = 0; i < C->size(); i++)
  {
    // enforce full precision reduction on C_otherDs due to numerical sensitivity
    PsiValueType product = (*C)[i];
    for (size_t id = 0; id < Dets.size(); id++)
      if (id != ig)
        product *= Dets[id]->getRatiosToRefDet()[(*C2node)[id][i]];
    C_otherDs[ig][(*C2node)[ig][i]] += product;
  }
  //put C_otherDs in host
  auto* C_otherDs_ptr = C_otherDs[ig].data();
  PRAGMA_OFFLOAD("omp target update to(C_otherDs_ptr[:Dets[ig]->getNumDets()])") //transfer content to device
}


} // namespace qmcplusplus
