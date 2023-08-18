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
#include "Platforms/OMPTarget/ompReductionComplex.hpp"

namespace qmcplusplus
{
struct MultiSlaterDetTableMethod::MultiSlaterDetTableMethodMultiWalkerResource : public Resource
{
  MultiSlaterDetTableMethodMultiWalkerResource() : Resource("MultiSlaterDetTableMethod") {}
  MultiSlaterDetTableMethodMultiWalkerResource(const MultiSlaterDetTableMethodMultiWalkerResource&)
      : MultiSlaterDetTableMethodMultiWalkerResource()
  {}

  std::unique_ptr<Resource> makeClone() const override
  {
    return std::make_unique<MultiSlaterDetTableMethodMultiWalkerResource>(*this);
  }

  /// grads of each unique determinants for multiple walkers
  Matrix<ValueType, OffloadAllocator<ValueType>> mw_grads;
  /// a collection of device pointers of multiple walkers fused for fast H2D transfer.
  OffloadVector<const ValueType*> C_otherDs_ptr_list;
  OffloadVector<const ValueType*> det_value_ptr_list;
};

MultiSlaterDetTableMethod::MultiSlaterDetTableMethod(ParticleSet& targetPtcl,
                                                     std::vector<std::unique_ptr<MultiDiracDeterminant>>&& dets,
                                                     bool use_pre_computing)
    : OptimizableObject("CI"),
      RatioTimer(createGlobalTimer(getClassName() + "::ratio")),
      offload_timer(createGlobalTimer(getClassName() + "::offload")),
      EvalGradTimer(createGlobalTimer(getClassName() + "::evalGrad")),
      RatioGradTimer(createGlobalTimer(getClassName() + "::ratioGrad")),
      PrepareGroupTimer(createGlobalTimer(getClassName() + "::prepareGroup")),
      UpdateTimer(createGlobalTimer(getClassName() + "::updateBuffer")),
      AccRejTimer(createGlobalTimer(getClassName() + "::Accept_Reject")),
      EvaluateTimer(createGlobalTimer(getClassName() + "::evaluate")),
      CI_Optimizable(false),
      use_pre_computing_(use_pre_computing)
{
  Dets = std::move(dets);
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

  clone->csf_data_ = csf_data_;

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
/*
void MultiSlaterDetTableMethod::mw_evaluateLog(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                      const RefVectorWithLeader<ParticleSet>& p_list,
                      const RefVector<ParticleSet::ParticleGradient>& G_list,
                      const RefVector<ParticleSet::ParticleLaplacian>& L_list)
{
	APP_ABORT("Iam In MW_LOG");
}
*/


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

  const auto& grads = (newpos) ? Dets[det_id]->getNewGrads() : Dets[det_id]->getGrads();
  const OffloadVector<ValueType>& detValues0 =
      (newpos) ? Dets[det_id]->getNewRatiosToRefDet() : Dets[det_id]->getRatiosToRefDet();
  const size_t noffset = Dets[det_id]->getFirstIndex();

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

  const auto& grads = (newpos) ? Dets[det_id]->getNewGrads() : Dets[det_id]->getGrads();
  const OffloadVector<ValueType>& detValues0 =
      (newpos) ? Dets[det_id]->getNewRatiosToRefDet() : Dets[det_id]->getRatiosToRefDet();
  const Matrix<ValueType>& spingrads = (newpos) ? Dets[det_id]->getNewSpinGrads() : Dets[det_id]->getSpinGrads();
  const size_t noffset               = Dets[det_id]->getFirstIndex();

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
  auto& det_leader = WFC_list.getCastedLeader<MultiSlaterDetTableMethod>();
  const int det_id = det_leader.getDetID(iat);
  const int nw     = WFC_list.size();
  const int ndets  = det_leader.Dets[det_id]->getNumDets();

  RefVectorWithLeader<MultiDiracDeterminant> det_list(*det_leader.Dets[det_id]);
  det_list.reserve(WFC_list.size());
  for (int iw = 0; iw < WFC_list.size(); iw++)
  {
    auto& det = WFC_list.getCastedElement<MultiSlaterDetTableMethod>(iw);
    det_list.push_back(*det.Dets[det_id]);
  }

  auto& mw_res   = det_leader.mw_res_handle_.getResource();
  auto& mw_grads = mw_res.mw_grads;
  mw_grads.resize(3 * nw, ndets);
  if (newpos)
    det_leader.Dets[det_id]->mw_evaluateDetsAndGradsForPtclMove(det_list, P_list, iat, mw_grads);
  else
    det_leader.Dets[det_id]->mw_evaluateGrads(det_list, P_list, iat, mw_grads);

  auto& det_value_ptr_list = mw_res.det_value_ptr_list;
  det_value_ptr_list.resize(nw);
  for (size_t iw = 0; iw < nw; iw++)
  {
    auto& det              = WFC_list.getCastedElement<MultiSlaterDetTableMethod>(iw);
    det_value_ptr_list[iw] = (newpos) ? det.Dets[det_id]->getNewRatiosToRefDet().device_data()
                                      : det.Dets[det_id]->getRatiosToRefDet().device_data();
  }

  std::vector<ValueType> grad_now_list(nw * 3, 0);
  auto* grad_now_list_ptr      = grad_now_list.data();
  auto* mw_grads_ptr           = mw_grads.data();
  auto* psi_list_ptr           = psi_list.data();
  auto* C_otherDs_ptr_list_ptr = mw_res.C_otherDs_ptr_list.data();
  auto* det_value_ptr_list_ptr = det_value_ptr_list.data();
  {
    ScopedTimer local_timer(det_leader.offload_timer);
    PRAGMA_OFFLOAD("omp target teams distribute map(from: psi_list_ptr[:nw]) \
                    map(from: grad_now_list_ptr[:3 * nw]) \
                    map(always, to: det_value_ptr_list_ptr[:nw]) \
                    map(to: mw_grads_ptr[:mw_grads.size()])")
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
        grad_local_x += C_otherDs_ptr_list_ptr[iw][i] * mw_grads_ptr[(3 * iw + 0) * ndets + i];
        grad_local_y += C_otherDs_ptr_list_ptr[iw][i] * mw_grads_ptr[(3 * iw + 1) * ndets + i];
        grad_local_z += C_otherDs_ptr_list_ptr[iw][i] * mw_grads_ptr[(3 * iw + 2) * ndets + i];
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

  spingrad += spingrad_iat;

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

  auto& det_leader = WFC_list.getCastedLeader<MultiSlaterDetTableMethod>();
  ScopedTimer local_timer(det_leader.EvalGradTimer);

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

  auto& det_leader = WFC_list.getCastedLeader<MultiSlaterDetTableMethod>();
  const int nw     = WFC_list.size();

  ScopedTimer local_timer(det_leader.RatioGradTimer);
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

WaveFunctionComponent::PsiValueType MultiSlaterDetTableMethod::computeRatio_NewMultiDet_to_NewRefDet(int det_id) const
{
  const auto& detValues0 = Dets[det_id]->getNewRatiosToRefDet();

  PsiValueType psi = 0;
  if (use_pre_computing_)
  {
    // This function computes
    // psi=Det_Coeff[i]*Det_Value[unique_det_up]*Det_Value[unique_det_dn]*Det_Value[unique_det_AnyOtherType]
    // Since only one electron group is moved at the time, identified by det_id, We precompute:
    // C_otherDs[det_id][i]=Det_Coeff[i]*Det_Value[unique_det_dn]*Det_Value[unique_det_AnyOtherType]
    for (size_t i = 0; i < Dets[det_id]->getNumDets(); i++)
      psi += detValues0[i] * C_otherDs[det_id][i];
  }
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
      psi += t;
    }
  }
  return psi;
}

// use ci_node for this routine only
WaveFunctionComponent::PsiValueType MultiSlaterDetTableMethod::ratio(ParticleSet& P, int iat)
{
  ScopedTimer local_timer(RatioTimer);
  UpdateMode = ORB_PBYP_RATIO;

  const int det_id = getDetID(iat);
  Dets[det_id]->evaluateDetsForPtclMove(P, iat);

  new_psi_ratio_to_new_ref_det_ = computeRatio_NewMultiDet_to_NewRefDet(det_id);
  curRatio = Dets[det_id]->getRefDetRatio() * new_psi_ratio_to_new_ref_det_ / psi_ratio_to_ref_det_;
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


  auto& det_leader = WFC_list.getCastedLeader<MultiSlaterDetTableMethod>();
  ScopedTimer local_timer(det_leader.RatioTimer);

  const int det_id = getDetID(iat);
  const int nw     = WFC_list.size();
  const int ndets  = det_leader.Dets[det_id]->getNumDets();

  RefVectorWithLeader<MultiDiracDeterminant> det_list(*det_leader.Dets[det_id]);
  det_list.reserve(WFC_list.size());
  for (int iw = 0; iw < WFC_list.size(); iw++)
  {
    auto& det = WFC_list.getCastedElement<MultiSlaterDetTableMethod>(iw);
    det_list.push_back(*det.Dets[det_id]);
  }

  det_leader.Dets[det_id]->mw_evaluateDetsForPtclMove(det_list, P_list, iat);

  auto& mw_res             = det_leader.mw_res_handle_.getResource();
  auto& det_value_ptr_list = mw_res.det_value_ptr_list;
  det_value_ptr_list.resize(nw);
  for (size_t iw = 0; iw < nw; iw++)
  {
    auto& det      = WFC_list.getCastedElement<MultiSlaterDetTableMethod>(iw);
    det.UpdateMode = ORB_PBYP_RATIO;

    det_value_ptr_list[iw] = det.Dets[det_id]->getNewRatiosToRefDet().device_data();
  }

  std::vector<PsiValueType> psi_list(nw, 0);
  auto* psi_list_ptr           = psi_list.data();
  auto* C_otherDs_ptr_list_ptr = mw_res.C_otherDs_ptr_list.data();
  auto* det_value_ptr_list_ptr = det_value_ptr_list.data();
  {
    ScopedTimer local_timer(det_leader.offload_timer);
    PRAGMA_OFFLOAD("omp target teams distribute map(always,from: psi_list_ptr[:nw]) \
          map(always, to: det_value_ptr_list_ptr[:nw])")
    for (size_t iw = 0; iw < nw; iw++)
    {
      PsiValueType psi_local(0);
      PRAGMA_OFFLOAD("omp parallel for reduction(+ : psi_local)")
      for (size_t i = 0; i < ndets; i++)
      {
        psi_local += det_value_ptr_list_ptr[iw][i] * C_otherDs_ptr_list_ptr[iw][i];
      }
      psi_list_ptr[iw] = psi_local;
    }
  }

  for (size_t iw = 0; iw < nw; iw++)
  {
    auto& det                         = WFC_list.getCastedElement<MultiSlaterDetTableMethod>(iw);
    det.new_psi_ratio_to_new_ref_det_ = psi_list[iw];
    ratios[iw] = det.curRatio = det.Dets[det_id]->getRefDetRatio() * psi_list[iw] / det.psi_ratio_to_ref_det_;
  }
}

void MultiSlaterDetTableMethod::evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios)
{
  ScopedTimer local_timer(RatioTimer);

  const int det_id = getDetID(VP.refPtcl);

  for (size_t iat = 0; iat < VP.getTotalNum(); ++iat)
  {
    Dets[det_id]->evaluateDetsForPtclMove(VP, iat, VP.refPtcl);
    const OffloadVector<ValueType>& detValues0 = Dets[det_id]->getNewRatiosToRefDet();

    PsiValueType psiNew = computeRatio_NewMultiDet_to_NewRefDet(det_id);
    ratios[iat]         = Dets[det_id]->getRefDetRatio() * psiNew / psi_ratio_to_ref_det_;
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

void MultiSlaterDetTableMethod::mw_accept_rejectMove(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                     const RefVectorWithLeader<ParticleSet>& p_list,
                                                     int iat,
                                                     const std::vector<bool>& isAccepted,
                                                     bool safe_to_delay) const
{
  ScopedTimer local_timer(AccRejTimer);
  for (size_t iw = 0; iw < isAccepted.size(); iw++)
  {
    auto& det = wfc_list.getCastedElement<MultiSlaterDetTableMethod>(iw);
    if (isAccepted[iw])
    {
      det.psi_ratio_to_ref_det_ = det.new_psi_ratio_to_new_ref_det_;
      det.log_value_ += convertValueToLog(det.curRatio);
    }
    det.curRatio = 1.0;
  }
  const auto det_id = getDetID(iat);
  const auto det_list(extract_DetRef_list(wfc_list, det_id));
  Dets[det_id]->mw_accept_rejectMove(det_list, p_list, iat, isAccepted);
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

void MultiSlaterDetTableMethod::extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs)
{
  opt_obj_refs.push_back(*this);
  for (int i = 0; i < Dets.size(); i++)
    Dets[i]->extractOptimizableObjectRefs(opt_obj_refs);
}

void MultiSlaterDetTableMethod::checkInVariablesExclusive(opt_variables_type& active)
{
  if (CI_Optimizable && myVars->size())
  {
    myVars->setIndexDefault();
    active.insertFrom(*myVars);
  }
}

void MultiSlaterDetTableMethod::checkOutVariables(const opt_variables_type& active)
{
  if (CI_Optimizable)
    myVars->getIndex(active);

  for (size_t id = 0; id < Dets.size(); id++)
    if (Dets[id]->isOptimizable())
      Dets[id]->checkOutVariables(active);
}

void MultiSlaterDetTableMethod::resetParametersExclusive(const opt_variables_type& active)
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
}

void MultiSlaterDetTableMethod::evaluateDerivatives(ParticleSet& P,
                                                    const opt_variables_type& optvars,
                                                    Vector<ValueType>& dlogpsi,
                                                    Vector<ValueType>& dhpsioverpsi)
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
                const OffloadVector<ValueType>& detValues_otherspin = Dets[other_id]->getRatiosToRefDet();
                size_t otherspinC                                   = (*C2node)[other_id][cnt];
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
                                                                         Vector<ValueType>& dlogpsi,
                                                                         Vector<ValueType>& dhpsioverpsi)
{
  //Currently, the MultiDiracDeterminant::evaluateDerivatives works with a legacy design, essentially requiring only up and down determinants.
  //e.g. for spinor cases, we only have one determinant so this interface doesn't work.
  //Here we throw an error only if the optimization is turned on for MultiDiracDeterminants until the code is updated
  for (auto const& det : Dets)
    if (!det->isOptimizable())
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
                                                      Vector<ValueType>& dlogpsi)
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

    if (recalculate)
    {
      Vector<ValueType> dlogpsi_local;
      evaluateDerivativesMSD(psi_ratio_to_ref_det_, dlogpsi_local);

      const size_t nparams = csf_data_ ? csf_data_->coeffs.size() - 1 : C->size() - 1;
      assert(dlogpsi_local.size() == nparams);
      for (int i = 0; i < nparams; i++)
      {
        int kk = myVars->where(i);
        if (kk < 0)
          continue;
        dlogpsi[kk] = dlogpsi_local[i];
      }
    }
  }

  evaluateMultiDiracDeterminantDerivativesWF(P, optvars, dlogpsi);
}

void MultiSlaterDetTableMethod::evaluateDerivativesMSD(const PsiValueType& multi_det_to_ref,
                                                       Vector<ValueType>& dlogpsi,
                                                       int det_id) const
{
  const bool newpos = det_id < 0 ? false : true;
  // when not using a new position, the result doesn't get affected by det_id, thus choose 0.
  if (det_id < 0)
    det_id = 0;

  ValueType psiinv       = static_cast<ValueType>(PsiValueType(1.0) / multi_det_to_ref);
  const auto& detValues0 = newpos ? Dets[det_id]->getNewRatiosToRefDet() : Dets[det_id]->getRatiosToRefDet();

  if (csf_data_) // CSF
  {
    dlogpsi.resize(csf_data_->coeffs.size() - 1);
    // this one is not optimizable
    int cnt = csf_data_->dets_per_csf[0];
    for (int i = 1; i < csf_data_->coeffs.size(); i++)
    {
      ValueType cdet = 0.0;
      for (int k = 0; k < csf_data_->dets_per_csf[i]; k++)
      {
        ValueType t = csf_data_->expansion[cnt] * psiinv * detValues0[(*C2node)[det_id][cnt]];
        // assume that evaluateLog has been called in opt routine before
        for (size_t id = 0; id < Dets.size(); id++)
          if (id != det_id)
            t *= Dets[id]->getRatiosToRefDet()[(*C2node)[id][cnt]];
        cdet += t;
        cnt++;
      }
      dlogpsi[i - 1] = cdet;
    }
  }
  else // CI
  {
    dlogpsi.resize(C->size() - 1);
    for (size_t i = 1; i < C->size(); i++)
    {
      ValueType cdet = psiinv * detValues0[(*C2node)[det_id][i]];
      // assume that evaluateLog has been called in opt routine before
      for (size_t id = 0; id < Dets.size(); id++)
        if (id != det_id)
          cdet *= Dets[id]->getRatiosToRefDet()[(*C2node)[id][i]];
      dlogpsi[i - 1] = cdet;
    }
  }
}

void MultiSlaterDetTableMethod::evaluateDerivRatios(const VirtualParticleSet& VP,
                                                    const opt_variables_type& optvars,
                                                    std::vector<ValueType>& ratios,
                                                    Matrix<ValueType>& dratios)
{
  const int det_id = getDetID(VP.refPtcl);

  bool recalculate(false);
  if (CI_Optimizable)
    for (int k = 0; k < myVars->size(); ++k)
    {
      int kk = myVars->where(k);
      if (kk < 0)
        continue;
      if (optvars.recompute(kk))
        recalculate = true;
    }

  // calculate derivatives based on the reference electron position
  Vector<ValueType> dlogpsi_ref, dlogpsi_vp;
  if (recalculate)
    evaluateDerivativesMSD(psi_ratio_to_ref_det_, dlogpsi_ref);

  for (size_t iat = 0; iat < VP.getTotalNum(); ++iat)
  {
    Dets[det_id]->evaluateDetsForPtclMove(VP, iat, VP.refPtcl);
    const OffloadVector<ValueType>& detValues0 = Dets[det_id]->getNewRatiosToRefDet();

    // calculate VP ratios
    PsiValueType psiNew = computeRatio_NewMultiDet_to_NewRefDet(det_id);
    ratios[iat]         = Dets[det_id]->getRefDetRatio() * psiNew / psi_ratio_to_ref_det_;

    // calculate VP ratios derivatives
    if (recalculate)
    {
      evaluateDerivativesMSD(psiNew, dlogpsi_vp, det_id);

      const size_t nparams = csf_data_ ? csf_data_->coeffs.size() - 1 : C->size() - 1;
      assert(dlogpsi_vp.size() == nparams);

      for (int i = 0; i < nparams; i++)
      {
        int kk = myVars->where(i);
        if (kk < 0)
          continue;
        dratios[iat][kk] = dlogpsi_vp[i] - dlogpsi_ref[i];
      }
    }
  }
}

void MultiSlaterDetTableMethod::evaluateMultiDiracDeterminantDerivativesWF(ParticleSet& P,
                                                                           const opt_variables_type& optvars,
                                                                           Vector<ValueType>& dlogpsi)
{
  //Currently, the MultiDiracDeterminant::evaluateDerivativesWF works with a legacy design, essentially requiring only up and down determinants.
  //e.g. for spinor cases, we only have one determinant so this interface doesn't work.
  //Here we throw an error only if the optimization is turned on for MultiDiracDeterminants until the code is updated
  for (auto const& det : Dets)
    if (!det->isOptimizable())
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

void MultiSlaterDetTableMethod::prepareGroup(ParticleSet& P, int ig)
{
  if (!use_pre_computing_)
    return;
  precomputeC_otherDs(P, ig);
}

void MultiSlaterDetTableMethod::mw_prepareGroup(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                const RefVectorWithLeader<ParticleSet>& p_list,
                                                int ig) const
{
  if (!use_pre_computing_)
    return;

  auto& det_leader = wfc_list.getCastedLeader<MultiSlaterDetTableMethod>();
  const size_t nw  = wfc_list.size();
  assert(this == &det_leader);

  auto& C_otherDs_ptr_list = det_leader.mw_res_handle_.getResource().C_otherDs_ptr_list;
  C_otherDs_ptr_list.resize(nw);
  for (int iw = 0; iw < nw; iw++)
  {
    auto& det = wfc_list.getCastedElement<MultiSlaterDetTableMethod>(iw);
    det.prepareGroup(p_list[iw], ig);
    C_otherDs_ptr_list[iw] = det.C_otherDs[ig].device_data();
  }
  C_otherDs_ptr_list.updateTo();
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
  //put C_otherDs in device
  C_otherDs[ig].updateTo();
}

void MultiSlaterDetTableMethod::createResource(ResourceCollection& collection) const
{
  collection.addResource(std::make_unique<MultiSlaterDetTableMethodMultiWalkerResource>());
  for (auto& det : Dets)
    det->createResource(collection);
}

void MultiSlaterDetTableMethod::acquireResource(ResourceCollection& collection,
                                                const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const
{
  auto& wfc_leader          = wfc_list.getCastedLeader<MultiSlaterDetTableMethod>();
  wfc_leader.mw_res_handle_ = collection.lendResource<MultiSlaterDetTableMethodMultiWalkerResource>();
  for (int idet = 0; idet < Dets.size(); idet++)
  {
    const auto det_list(extract_DetRef_list(wfc_list, idet));
    Dets[idet]->acquireResource(collection, det_list);
  }
}

void MultiSlaterDetTableMethod::releaseResource(ResourceCollection& collection,
                                                const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const
{
  auto& wfc_leader = wfc_list.getCastedLeader<MultiSlaterDetTableMethod>();
  collection.takebackResource(wfc_leader.mw_res_handle_);
  for (int idet = 0; idet < Dets.size(); idet++)
  {
    const auto det_list(extract_DetRef_list(wfc_list, idet));
    Dets[idet]->releaseResource(collection, det_list);
  }
}
} // namespace qmcplusplus
