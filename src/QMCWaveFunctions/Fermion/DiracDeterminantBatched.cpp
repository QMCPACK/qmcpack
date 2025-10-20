//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "DiracDeterminantBatched.h"
#include <cassert>
#include "Numerics/DeterminantOperators.h"
#include "CPU/BLAS.hpp"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Numerics/MatrixOperators.h"
#include "QMCWaveFunctions/TWFFastDerivWrapper.h"
#ifndef QMC_COMPLEX
#include "QMCWaveFunctions/RotatedSPOs.h"
#endif
#include "CPU/SIMD/inner_product.hpp"
#include "DiracMatrixInverterOMPTarget.hpp"
#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
#include "DiracMatrixInverterCUDA.hpp"
#endif

namespace qmcplusplus
{

template<PlatformKind UEPL, typename FPVT, typename VT>
struct DetInverterSelector
{
  using Inverter = DiracMatrixInverterOMPTarget<FPVT, VT>;
};

#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
template<typename FPVT, typename VT>
struct DetInverterSelector<PlatformKind::CUDA, FPVT, VT>
{
  using Inverter = DiracMatrixInverterCUDA<FPVT, VT>;
};
#endif

template<PlatformKind PL, typename VT, typename FPVT>
struct DiracDeterminantBatched<PL, VT, FPVT>::DiracDeterminantBatchedMultiWalkerResource : public Resource
{
  DiracDeterminantBatchedMultiWalkerResource() : Resource("DiracDeterminantBatched") {}
  DiracDeterminantBatchedMultiWalkerResource(const DiracDeterminantBatchedMultiWalkerResource&)
      : DiracDeterminantBatchedMultiWalkerResource()
  {}

  std::unique_ptr<Resource> makeClone() const override
  {
    return std::make_unique<DiracDeterminantBatchedMultiWalkerResource>(*this);
  }
  DualVector<LogValue> log_values;
  /// value, grads, laplacian of single-particle orbital for particle-by-particle update and multi walker [5][nw][norb]
  OffloadMWVGLArray phi_vgl_v;
  /// multi walker of ratio
  std::vector<Value> ratios_local;
  /// multi walker of grads
  std::vector<Grad> grad_new_local;
  /// multi walker of spingrads
  std::vector<Value> spingrad_new_local;
  /// mw spin gradients of orbitals, matrix is [nw][norb]
  OffloadMatrix<ComplexType> mw_dspin;
  /// reference to per DDB psiMinvs in a crowd
  RefVector<DualMatrix<Value>> psiMinv_refs;
  ///
  typename UpdateEngine::MultiWalkerResource engine_rsc;
};

/** constructor
 *@param spos the single-particle orbital set
 *@param first index of the first particle
 */
template<PlatformKind PL, typename VT, typename FPVT>
DiracDeterminantBatched<PL, VT, FPVT>::DiracDeterminantBatched(SPOSet& phi,
                                                               int first,
                                                               int last,
                                                               int ndelay,
                                                               DetMatInvertor matrix_inverter_kind)
    : DiracDeterminantBase("DiracDeterminantBatched", phi, first, last),
      det_engine_(NumOrbitals, ndelay),
      ndelay_(ndelay),
      matrix_inverter_kind_(matrix_inverter_kind),
      D2HTimer(createGlobalTimer("DiracDeterminantBatched::D2H", timer_level_fine)),
      H2DTimer(createGlobalTimer("DiracDeterminantBatched::H2D", timer_level_fine))
{
  static_assert(std::is_same<SPOSet::ValueType, typename UpdateEngine::Value>::value);
  resize(NumPtcls, NumPtcls);

#ifndef QMC_COMPLEX
  RotatedSPOs* rot_spo = dynamic_cast<RotatedSPOs*>(&phi_);
  if (rot_spo)
    rot_spo->buildOptVariables(NumPtcls);
#endif
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::invertPsiM(const DualMatrix<Value>& psiM, DualMatrix<Value>& psiMinv)
{
  ScopedTimer inverse_timer(InverseTimer);
  host_inverter_.invert_transpose(psiM, psiMinv, log_value_);
  psiMinv.updateTo();

#ifndef NDEBUG
  // This is easily breakable in that it assumes this function gets psiMinv == det_engine_.psiMinv_
  auto& engine_psiMinv = psiMinv_;
  dummy_vmt.attachReference(engine_psiMinv.data(), engine_psiMinv.rows(), engine_psiMinv.cols());
#endif
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::mw_invertPsiM(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                          const RefVector<const DualMatrix<Value>>& logdetT_list,
                                                          const RefVector<DualMatrix<Value>>& a_inv_list)
{
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<PL, VT, FPVT>>();
  auto& mw_res     = wfc_leader.mw_res_handle_.getResource();
  ScopedTimer inverse_timer(wfc_leader.InverseTimer);
  const auto nw = wfc_list.size();

  mw_res.log_values.resize(nw);

  wfc_leader.accel_inverter_.getResource().mw_invert_transpose(mw_res.engine_rsc.queue, logdetT_list, a_inv_list,
                                                               mw_res.log_values);

  for (int iw = 0; iw < nw; ++iw)
  {
    auto& det      = wfc_list.getCastedElement<DiracDeterminantBatched<PL, VT, FPVT>>(iw);
    det.log_value_ = mw_res.log_values[iw];
  }

#ifndef NDEBUG
  for (int iw = 0; iw < nw; ++iw)
  {
    auto& det            = wfc_list.getCastedElement<DiracDeterminantBatched<PL, VT, FPVT>>(iw);
    auto& engine_psiMinv = det.psiMinv_;
    det.dummy_vmt.attachReference(engine_psiMinv.data(), engine_psiMinv.rows(), engine_psiMinv.cols());
  }
#endif
}

///reset the size: with the number of particles and number of orbtials
template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::resize(int nel, int morb)
{
  int norb = morb;
  if (norb <= 0)
    norb = nel; // for morb == -1 (default)
  psiMinv_.resize(norb, getAlignedSize<Value>(norb));
  psiM_vgl.resize(nel * norb);
  // attach pointers VGL
  psiM_temp.attachReference(psiM_vgl, psiM_vgl.data(0), nel, norb);
  psiM_host.attachReference(psiM_vgl.data(0), nel, norb);
  dpsiM.attachReference(reinterpret_cast<Grad*>(psiM_vgl.data(1)), nel, norb);
  d2psiM.attachReference(psiM_vgl.data(4), nel, norb);

  psiV.resize(NumOrbitals);
  psiV_host_view.attachReference(psiV.data(), NumOrbitals);
  dpsiV.resize(NumOrbitals);
  dpsiV_host_view.attachReference(dpsiV.data(), NumOrbitals);
  d2psiV.resize(NumOrbitals);
  d2psiV_host_view.attachReference(d2psiV.data(), NumOrbitals);
  dspin_psiV.resize(NumOrbitals);
  dspin_psiV_host_view.attachReference(dspin_psiV.data(), NumOrbitals);
}

template<PlatformKind PL, typename VT, typename FPVT>
typename DiracDeterminantBatched<PL, VT, FPVT>::Grad DiracDeterminantBatched<PL, VT, FPVT>::evalGrad(ParticleSet& P,
                                                                                                     int iat)
{
  ScopedTimer local_timer(RatioTimer);
  const int WorkingIndex = iat - FirstIndex;
  Grad g                 = simd::dot(psiMinv_[WorkingIndex], dpsiM[WorkingIndex], NumOrbitals);
  return g;
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::mw_evalGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                        const RefVectorWithLeader<ParticleSet>& p_list,
                                                        int iat,
                                                        std::vector<Grad>& grad_now) const
{
  assert(this == &wfc_list.getLeader());
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<PL, VT, FPVT>>();
  auto& mw_res     = wfc_leader.mw_res_handle_.getResource();
  ScopedTimer local_timer(RatioTimer);

  const int nw = wfc_list.size();
  std::vector<const Value*> dpsiM_row_list(nw, nullptr);
  RefVectorWithLeader<UpdateEngine> engine_list(wfc_leader.det_engine_);
  engine_list.reserve(nw);

  const int WorkingIndex = iat - FirstIndex;
  for (int iw = 0; iw < nw; iw++)
  {
    auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<PL, VT, FPVT>>(iw);
    // capacity is the size of each vector in the VGL so this advances us to the g then makes
    // an offset into the gradients
    dpsiM_row_list[iw] = det.psiM_vgl.device_data() + psiM_vgl.capacity() + NumOrbitals * WorkingIndex * DIM;
    engine_list.push_back(det.det_engine_);
  }

  UpdateEngine::mw_evalGrad(engine_list, wfc_leader.mw_res_handle_.getResource().engine_rsc, mw_res.psiMinv_refs,
                            dpsiM_row_list, WorkingIndex, grad_now);
}

template<PlatformKind PL, typename VT, typename FPVT>
typename DiracDeterminantBatched<PL, VT, FPVT>::Grad DiracDeterminantBatched<PL, VT, FPVT>::evalGradWithSpin(
    ParticleSet& P,
    int iat,
    ComplexType& spingrad)
{
  phi_.evaluate_spin(P, iat, psiV_host_view, dspin_psiV_host_view);
  ScopedTimer local_timer(RatioTimer);
  const int WorkingIndex = iat - FirstIndex;
  Grad g                 = simd::dot(psiMinv_[WorkingIndex], dpsiM[WorkingIndex], NumOrbitals);
  ComplexType spin_g     = simd::dot(psiMinv_[WorkingIndex], dspin_psiV.data(), NumOrbitals);
  spingrad += spin_g;
  return g;
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::mw_evalGradWithSpin(
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
    const RefVectorWithLeader<ParticleSet>& p_list,
    int iat,
    std::vector<Grad>& grad_now,
    std::vector<ComplexType>& spingrad_now) const
{
  assert(this == &wfc_list.getLeader());
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<PL, VT, FPVT>>();
  auto& mw_res     = wfc_leader.mw_res_handle_.getResource();
  auto& mw_dspin   = mw_res.mw_dspin;
  ScopedTimer local_timer(RatioTimer);

  const int nw           = wfc_list.size();
  const int num_orbitals = wfc_leader.phi_.size();
  mw_dspin.resize(nw, num_orbitals);

  //Here, we are just always recomputing the spin gradients from the SPOSet for simplicity.
  //If we stored and modified the accept/reject to include updating stored spin gradients, we could the
  //mw_evaluateVGLWithSpin call below and just use the stored spin gradients.
  //May revisit this in the future.
  RefVectorWithLeader<SPOSet> phi_list(phi_);
  RefVector<SPOSet::ValueVector> psi_v_list, lap_v_list;
  RefVector<SPOSet::GradVector> grad_v_list;
  for (int iw = 0; iw < wfc_list.size(); iw++)
  {
    auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<PL, VT, FPVT>>(iw);
    phi_list.push_back(det.phi_);
    psi_v_list.push_back(det.psiV_host_view);
    grad_v_list.push_back(det.dpsiV_host_view);
    lap_v_list.push_back(det.d2psiV_host_view);
  }

  auto& phi_leader = phi_list.getLeader();
  phi_leader.mw_evaluateVGLWithSpin(phi_list, p_list, iat, psi_v_list, grad_v_list, lap_v_list, mw_dspin);

  std::vector<const Value*> dpsiM_row_list(nw, nullptr);
  RefVectorWithLeader<UpdateEngine> engine_list(wfc_leader.det_engine_);
  engine_list.reserve(nw);

  const int WorkingIndex = iat - FirstIndex;
  for (int iw = 0; iw < nw; iw++)
  {
    auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<PL, VT, FPVT>>(iw);
    // capacity is the size of each vector in the VGL so this advances us to the g then makes
    // an offset into the gradients
    dpsiM_row_list[iw] = det.psiM_vgl.device_data() + psiM_vgl.capacity() + NumOrbitals * WorkingIndex * DIM;
    engine_list.push_back(det.det_engine_);
  }

  UpdateEngine::mw_evalGradWithSpin(engine_list, wfc_leader.mw_res_handle_.getResource().engine_rsc,
                                    mw_res.psiMinv_refs, dpsiM_row_list, mw_dspin, WorkingIndex, grad_now,
                                    spingrad_now);
}

template<PlatformKind PL, typename VT, typename FPVT>
typename DiracDeterminantBatched<PL, VT, FPVT>::PsiValue DiracDeterminantBatched<PL, VT, FPVT>::ratioGrad(
    ParticleSet& P,
    int iat,
    Grad& grad_iat)
{
  UpdateMode = ORB_PBYP_PARTIAL;

  {
    ScopedTimer local_timer(SPOVGLTimer);
    phi_.evaluateVGL(P, iat, psiV_host_view, dpsiV_host_view, d2psiV_host_view);
  }

  {
    ScopedTimer local_timer(RatioTimer);
    const int WorkingIndex = iat - FirstIndex;
    curRatio               = simd::dot(psiMinv_[WorkingIndex], psiV.data(), NumOrbitals);
    grad_iat += static_cast<Value>(static_cast<PsiValue>(1.0) / curRatio) *
        simd::dot(psiMinv_[WorkingIndex], dpsiV.data(), NumOrbitals);
  }
  return curRatio;
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::mw_ratioGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                         const RefVectorWithLeader<ParticleSet>& p_list,
                                                         int iat,
                                                         std::vector<PsiValue>& ratios,
                                                         std::vector<Grad>& grad_new) const
{
  assert(this == &wfc_list.getLeader());
  auto& wfc_leader     = wfc_list.getCastedLeader<DiracDeterminantBatched<PL, VT, FPVT>>();
  auto& mw_res         = wfc_leader.mw_res_handle_.getResource();
  auto& phi_vgl_v      = mw_res.phi_vgl_v;
  auto& ratios_local   = mw_res.ratios_local;
  auto& grad_new_local = mw_res.grad_new_local;
  {
    ScopedTimer local_timer(SPOVGLTimer);
    RefVectorWithLeader<SPOSet> phi_list(phi_);
    phi_list.reserve(wfc_list.size());
    RefVectorWithLeader<UpdateEngine> engine_list(wfc_leader.det_engine_);
    engine_list.reserve(wfc_list.size());
    const int WorkingIndex = iat - FirstIndex;
    for (int iw = 0; iw < wfc_list.size(); iw++)
    {
      auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<PL, VT, FPVT>>(iw);
      phi_list.push_back(det.phi_);
      engine_list.push_back(det.det_engine_);
    }

    auto psiMinv_row_dev_ptr_list =
        UpdateEngine::mw_getInvRow(engine_list, wfc_leader.mw_res_handle_.getResource().engine_rsc, mw_res.psiMinv_refs,
                                   WorkingIndex, !phi_.isOMPoffload());

    phi_vgl_v.resize(SPOSet::DIM_VGL, wfc_list.size(), NumOrbitals);
    ratios_local.resize(wfc_list.size());
    grad_new_local.resize(wfc_list.size());

    wfc_leader.phi_.mw_evaluateVGLandDetRatioGrads(phi_list, p_list, iat, psiMinv_row_dev_ptr_list, phi_vgl_v,
                                                   ratios_local, grad_new_local);
  }

  wfc_leader.UpdateMode = ORB_PBYP_PARTIAL;
  for (int iw = 0; iw < wfc_list.size(); iw++)
  {
    auto& det      = wfc_list.getCastedElement<DiracDeterminantBatched<PL, VT, FPVT>>(iw);
    det.UpdateMode = ORB_PBYP_PARTIAL;
    ratios[iw] = det.curRatio = ratios_local[iw];
    grad_new[iw] += grad_new_local[iw];
  }
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::mw_ratioGradWithSpin(
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
    const RefVectorWithLeader<ParticleSet>& p_list,
    int iat,
    std::vector<PsiValue>& ratios,
    std::vector<Grad>& grad_new,
    std::vector<ComplexType>& spingrad_new) const
{
  assert(this == &wfc_list.getLeader());
  auto& wfc_leader         = wfc_list.getCastedLeader<DiracDeterminantBatched<PL, VT, FPVT>>();
  auto& mw_res             = wfc_leader.mw_res_handle_.getResource();
  auto& phi_vgl_v          = mw_res.phi_vgl_v;
  auto& ratios_local       = mw_res.ratios_local;
  auto& grad_new_local     = mw_res.grad_new_local;
  auto& spingrad_new_local = mw_res.spingrad_new_local;
  {
    ScopedTimer local_timer(SPOVGLTimer);
    RefVectorWithLeader<SPOSet> phi_list(phi_);
    phi_list.reserve(wfc_list.size());
    RefVectorWithLeader<UpdateEngine> engine_list(wfc_leader.det_engine_);
    engine_list.reserve(wfc_list.size());
    const int WorkingIndex = iat - FirstIndex;
    for (int iw = 0; iw < wfc_list.size(); iw++)
    {
      auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<PL, VT, FPVT>>(iw);
      phi_list.push_back(det.phi_);
      engine_list.push_back(det.det_engine_);
    }

    auto psiMinv_row_dev_ptr_list =
        UpdateEngine::mw_getInvRow(engine_list, wfc_leader.mw_res_handle_.getResource().engine_rsc, mw_res.psiMinv_refs,
                                   WorkingIndex, !phi_.isOMPoffload());

    phi_vgl_v.resize(SPOSet::DIM_VGL, wfc_list.size(), NumOrbitals);
    ratios_local.resize(wfc_list.size());
    grad_new_local.resize(wfc_list.size());
    spingrad_new_local.resize(wfc_list.size());

    wfc_leader.phi_.mw_evaluateVGLandDetRatioGradsWithSpin(phi_list, p_list, iat, psiMinv_row_dev_ptr_list, phi_vgl_v,
                                                           ratios_local, grad_new_local, spingrad_new_local);
  }

  wfc_leader.UpdateMode = ORB_PBYP_PARTIAL;
  for (int iw = 0; iw < wfc_list.size(); iw++)
  {
    auto& det      = wfc_list.getCastedElement<DiracDeterminantBatched<PL, VT, FPVT>>(iw);
    det.UpdateMode = ORB_PBYP_PARTIAL;
    ratios[iw] = det.curRatio = ratios_local[iw];
    grad_new[iw] += grad_new_local[iw];
    spingrad_new[iw] += spingrad_new_local[iw];
  }
}

template<PlatformKind PL, typename VT, typename FPVT>
typename DiracDeterminantBatched<PL, VT, FPVT>::PsiValue DiracDeterminantBatched<PL, VT, FPVT>::ratioGradWithSpin(
    ParticleSet& P,
    int iat,
    Grad& grad_iat,
    ComplexType& spingrad_iat)
{
  UpdateMode = ORB_PBYP_PARTIAL;

  {
    ScopedTimer local_timer(SPOVGLTimer);
    phi_.evaluateVGL_spin(P, iat, psiV_host_view, dpsiV_host_view, d2psiV_host_view, dspin_psiV_host_view);
  }

  {
    ScopedTimer local_timer(RatioTimer);
    const int WorkingIndex = iat - FirstIndex;
    curRatio               = simd::dot(psiMinv_[WorkingIndex], psiV.data(), NumOrbitals);
    grad_iat += static_cast<Value>(static_cast<PsiValue>(1.0) / curRatio) *
        simd::dot(psiMinv_[WorkingIndex], dpsiV.data(), NumOrbitals);
    spingrad_iat += static_cast<Value>(static_cast<PsiValue>(1.0) / curRatio) *
        simd::dot(psiMinv_[WorkingIndex], dspin_psiV.data(), NumOrbitals);
  }
  return curRatio;
}


/** move was accepted, update the real container
*/
template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::acceptMove(ParticleSet& P, int iat, bool safe_to_delay)
{
  if (curRatio == PsiValue(0))
  {
    std::ostringstream msg;
    msg << "DiracDeterminantBatched::acceptMove curRatio is " << curRatio << "! Report a bug." << std::endl;
    throw std::runtime_error(msg.str());
  }
  const int WorkingIndex = iat - FirstIndex;
  log_value_ += convertValueToLog(curRatio);
  {
    ScopedTimer local_timer(UpdateTimer);
    psiV.updateTo();
    det_engine_.updateRow(psiMinv_, WorkingIndex, psiV, curRatio);
    if (UpdateMode == ORB_PBYP_PARTIAL)
    {
      simd::copy(dpsiM[WorkingIndex], dpsiV.data(), NumOrbitals);
      simd::copy(d2psiM[WorkingIndex], d2psiV.data(), NumOrbitals);
    }
  }
  curRatio = 1.0;
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::mw_accept_rejectMove(
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
    const RefVectorWithLeader<ParticleSet>& p_list,
    int iat,
    const std::vector<bool>& isAccepted,
    bool safe_to_delay) const
{
  assert(this == &wfc_list.getLeader());
  auto& wfc_leader   = wfc_list.getCastedLeader<DiracDeterminantBatched<PL, VT, FPVT>>();
  auto& mw_res       = wfc_leader.mw_res_handle_.getResource();
  auto& phi_vgl_v    = mw_res.phi_vgl_v;
  auto& ratios_local = mw_res.ratios_local;

  ScopedTimer update(UpdateTimer);

  const int nw = wfc_list.size();
  int count    = 0;
  for (int iw = 0; iw < nw; iw++)
    if (isAccepted[iw])
      count++;
  const int n_accepted = count;

  RefVectorWithLeader<UpdateEngine> engine_list(wfc_leader.det_engine_);
  engine_list.reserve(nw);
  std::vector<Value*> psiM_g_dev_ptr_list(n_accepted, nullptr);
  std::vector<Value*> psiM_l_dev_ptr_list(n_accepted, nullptr);

  const int WorkingIndex = iat - FirstIndex;
  for (int iw = 0, count = 0; iw < nw; iw++)
  {
    // This can be auto but some debuggers can't figure the type out.
    DiracDeterminantBatched<PL, VT, FPVT>& det = wfc_list.getCastedElement<DiracDeterminantBatched<PL, VT, FPVT>>(iw);
    engine_list.push_back(det.det_engine_);
    if (isAccepted[iw])
    {
      psiM_g_dev_ptr_list[count] = det.psiM_vgl.device_data() + psiM_vgl.capacity() + NumOrbitals * WorkingIndex * DIM;
      psiM_l_dev_ptr_list[count] = det.psiM_vgl.device_data() + psiM_vgl.capacity() * 4 + NumOrbitals * WorkingIndex;
      if (det.curRatio == PsiValue(0))

      {
        std::ostringstream msg;
        msg << "DiracDeterminantBatched::mw_accept_rejectMove det.curRatio is " << det.curRatio << "! Report a bug."
            << std::endl;
        throw std::runtime_error(msg.str());
      }
      det.log_value_ += convertValueToLog(det.curRatio);
      count++;
    }
    det.curRatio = 1.0;
  }

  UpdateEngine::mw_accept_rejectRow(engine_list, wfc_leader.mw_res_handle_.getResource().engine_rsc,
                                    mw_res.psiMinv_refs, WorkingIndex, psiM_g_dev_ptr_list, psiM_l_dev_ptr_list,
                                    isAccepted, phi_vgl_v, ratios_local);

  if (!safe_to_delay)
    UpdateEngine::mw_updateInvMat(engine_list, wfc_leader.mw_res_handle_.getResource().engine_rsc, mw_res.psiMinv_refs);
}

/** move was rejected. copy the real container to the temporary to move on
*/
template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::restore(int iat)
{
  curRatio = 1.0;
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::completeUpdates()
{
  ScopedTimer update(UpdateTimer);
  if (UpdateMode == ORB_PBYP_PARTIAL)
  {
    // dpsiM, d2psiM on the device needs to be aligned as host.
    auto* psiM_vgl_ptr = psiM_vgl.data();
    // transfer host to device, total size 4, g(3) + l(1)
    PRAGMA_OFFLOAD("omp target update to(psiM_vgl_ptr[psiM_vgl.capacity():psiM_vgl.capacity()*4])")
  }
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::mw_completeUpdates(
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const
{
  assert(this == &wfc_list.getLeader());
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<PL, VT, FPVT>>();
  auto& mw_res     = wfc_leader.mw_res_handle_.getResource();
  const auto nw    = wfc_list.size();
  RefVectorWithLeader<UpdateEngine> engine_list(wfc_leader.det_engine_);
  engine_list.reserve(nw);
  for (int iw = 0; iw < nw; iw++)
  {
    auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<PL, VT, FPVT>>(iw);
    engine_list.push_back(det.det_engine_);
  }

  {
    ScopedTimer update(UpdateTimer);
    UpdateEngine::mw_updateInvMat(engine_list, mw_res.engine_rsc, mw_res.psiMinv_refs);
  }

  { // transfer dpsiM, d2psiM, psiMinv to host
    ScopedTimer d2h(D2HTimer);

    // this call also completes all the device copying of dpsiM, d2psiM before the target update
    UpdateEngine::mw_transferAinv_D2H(engine_list, mw_res.engine_rsc, mw_res.psiMinv_refs);

    if (UpdateMode == ORB_PBYP_PARTIAL)
    {
      RefVector<DualVGLVector> psiM_vgl_list;
      psiM_vgl_list.reserve(nw);
      for (int iw = 0; iw < nw; iw++)
      {
        auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<PL, VT, FPVT>>(iw);
        psiM_vgl_list.push_back(det.psiM_vgl);
      }

      auto& queue = mw_res.engine_rsc.queue;
      for (DualVGLVector& psiM_vgl : psiM_vgl_list)
      {
        const size_t stride = psiM_vgl.capacity();
        // transfer device to host, total size 4, g(3) + l(1), skipping v
        queue.enqueueD2H(psiM_vgl, stride * 4, stride);
      }
      queue.sync();
    }
  }
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::computeGL(ParticleSet::ParticleGradient& G,
                                                      ParticleSet::ParticleLaplacian& L) const
{
  for (size_t i = 0, iat = FirstIndex; i < NumPtcls; ++i, ++iat)
  {
    Grad rv   = simd::dot(psiMinv_[i], dpsiM[i], NumOrbitals);
    Value lap = simd::dot(psiMinv_[i], d2psiM[i], NumOrbitals);
    G[iat] += rv;
    L[iat] += lap - dot(rv, rv);
  }
}

template<PlatformKind PL, typename VT, typename FPVT>
typename DiracDeterminantBatched<PL, VT, FPVT>::LogValue DiracDeterminantBatched<PL, VT, FPVT>::evaluateGL(
    const ParticleSet& P,
    ParticleSet::ParticleGradient& G,
    ParticleSet::ParticleLaplacian& L,
    bool fromscratch)
{
  if (fromscratch)
    // this updates LogValue
    evaluateLog(P, G, L);
  else
  {
    if (UpdateMode == ORB_PBYP_RATIO)
    { //need to compute dpsiM and d2psiM. Do not touch psiM!
      ScopedTimer spo_timer(SPOVGLTimer);
      phi_.evaluate_notranspose(P, FirstIndex, LastIndex, psiM_host, dpsiM, d2psiM);
    }
    UpdateMode = ORB_WALKER;
    computeGL(G, L);
  }
  return log_value_;
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::mw_evaluateGL(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                          const RefVectorWithLeader<ParticleSet>& p_list,
                                                          const RefVector<ParticleSet::ParticleGradient>& G_list,
                                                          const RefVector<ParticleSet::ParticleLaplacian>& L_list,
                                                          bool fromscratch) const
{
  assert(this == &wfc_list.getLeader());
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<PL, VT, FPVT>>();
  if (fromscratch)
    mw_evaluateLog(wfc_list, p_list, G_list, L_list);
  else
  {
    const auto nw = wfc_list.size();
    RefVectorWithLeader<UpdateEngine> engine_list(wfc_leader.det_engine_);
    engine_list.reserve(nw);

    if (UpdateMode == ORB_PBYP_RATIO)
    { //need to compute dpsiM and d2psiM. psiMinv is not touched!
      ScopedTimer spo_timer(SPOVGLTimer);

      RefVectorWithLeader<SPOSet> phi_list(phi_);
      RefVector<Matrix<Value>> psiM_temp_list;
      RefVector<Matrix<Grad>> dpsiM_list;
      RefVector<Matrix<Value>> d2psiM_list;
      phi_list.reserve(wfc_list.size());
      psiM_temp_list.reserve(nw);
      dpsiM_list.reserve(nw);
      d2psiM_list.reserve(nw);

      for (int iw = 0; iw < nw; iw++)
      {
        auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<PL, VT, FPVT>>(iw);
        engine_list.push_back(det.det_engine_);
        phi_list.push_back(det.phi_);
        psiM_temp_list.push_back(det.psiM_host);
        dpsiM_list.push_back(det.dpsiM);
        d2psiM_list.push_back(det.d2psiM);
      }

      phi_.mw_evaluate_notranspose(phi_list, p_list, FirstIndex, LastIndex, psiM_temp_list, dpsiM_list, d2psiM_list);
    }

    for (int iw = 0; iw < nw; iw++)
    {
      auto& det      = wfc_list.getCastedElement<DiracDeterminantBatched<PL, VT, FPVT>>(iw);
      det.UpdateMode = ORB_WALKER;
      det.computeGL(G_list[iw], L_list[iw]);
    }
  }
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::registerData(ParticleSet& P, WFBufferType& buf)
{
  buf.add(psiMinv_.first_address(), psiMinv_.last_address());
  buf.add(dpsiM.first_address(), dpsiM.last_address());
  buf.add(d2psiM.first_address(), d2psiM.last_address());
  buf.add(log_value_);
}

template<PlatformKind PL, typename VT, typename FPVT>
typename DiracDeterminantBatched<PL, VT, FPVT>::LogValue DiracDeterminantBatched<PL, VT, FPVT>::updateBuffer(
    ParticleSet& P,
    WFBufferType& buf,
    bool fromscratch)
{
  evaluateGL(P, P.G, P.L, fromscratch);
  ScopedTimer local_timer(BufferTimer);
  buf.put(psiMinv_.first_address(), psiMinv_.last_address());
  buf.put(dpsiM.first_address(), dpsiM.last_address());
  buf.put(d2psiM.first_address(), d2psiM.last_address());
  buf.put(log_value_);
  return log_value_;
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  ScopedTimer local_timer(BufferTimer);
  buf.get(psiMinv_.first_address(), psiMinv_.last_address());
  buf.get(dpsiM.first_address(), dpsiM.last_address());
  buf.get(d2psiM.first_address(), d2psiM.last_address());
  psiMinv_.updateTo();
  auto* psiM_vgl_ptr           = psiM_vgl.data();
  const size_t psiM_vgl_stride = psiM_vgl.capacity();
  // transfer host to device, total size 4, g(3) + l(1)
  PRAGMA_OFFLOAD("omp target update to(psiM_vgl_ptr[psiM_vgl_stride:psiM_vgl_stride*4])")
  buf.get(log_value_);
}

/** return the ratio only for the  iat-th partcle move
 * @param P current configuration
 * @param iat the particle thas is being moved
 */
template<PlatformKind PL, typename VT, typename FPVT>
typename DiracDeterminantBatched<PL, VT, FPVT>::PsiValue DiracDeterminantBatched<PL, VT, FPVT>::ratio(ParticleSet& P,
                                                                                                      int iat)
{
  UpdateMode             = ORB_PBYP_RATIO;
  const int WorkingIndex = iat - FirstIndex;
  {
    ScopedTimer local_timer(SPOVTimer);
    phi_.evaluateValue(P, iat, psiV_host_view);
  }
  {
    ScopedTimer local_timer(RatioTimer);
    curRatio = simd::dot(psiMinv_[WorkingIndex], psiV.data(), NumOrbitals);
  }
  return curRatio;
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::mw_calcRatio(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                         const RefVectorWithLeader<ParticleSet>& p_list,
                                                         int iat,
                                                         std::vector<PsiValue>& ratios) const
{
  assert(this == &wfc_list.getLeader());
  auto& wfc_leader     = wfc_list.getCastedLeader<DiracDeterminantBatched<PL, VT, FPVT>>();
  auto& mw_res         = wfc_leader.mw_res_handle_.getResource();
  auto& phi_vgl_v      = mw_res.phi_vgl_v;
  auto& ratios_local   = mw_res.ratios_local;
  auto& grad_new_local = mw_res.grad_new_local;

  {
    ScopedTimer local_timer(SPOVTimer);
    RefVectorWithLeader<SPOSet> phi_list(phi_);
    phi_list.reserve(wfc_list.size());
    RefVectorWithLeader<UpdateEngine> engine_list(wfc_leader.det_engine_);
    engine_list.reserve(wfc_list.size());

    const int WorkingIndex = iat - FirstIndex;
    for (int iw = 0; iw < wfc_list.size(); iw++)
    {
      auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<PL, VT, FPVT>>(iw);
      phi_list.push_back(det.phi_);
      engine_list.push_back(det.det_engine_);
    }

    auto psiMinv_row_dev_ptr_list =
        UpdateEngine::mw_getInvRow(engine_list, wfc_leader.mw_res_handle_.getResource().engine_rsc, mw_res.psiMinv_refs,
                                   WorkingIndex, !phi_.isOMPoffload());

    phi_vgl_v.resize(SPOSet::DIM_VGL, wfc_list.size(), NumOrbitals);
    ratios_local.resize(wfc_list.size());
    grad_new_local.resize(wfc_list.size());

    // calling phi_.mw_evaluateVGLandDetRatioGrads is a temporary workaround.
    // We may implement mw_evaluateVandDetRatio in the future.
    wfc_leader.phi_.mw_evaluateVGLandDetRatioGrads(phi_list, p_list, iat, psiMinv_row_dev_ptr_list, phi_vgl_v,
                                                   ratios_local, grad_new_local);
  }

  wfc_leader.UpdateMode = ORB_PBYP_RATIO;
  for (int iw = 0; iw < wfc_list.size(); iw++)
  {
    auto& det      = wfc_list.getCastedElement<DiracDeterminantBatched<PL, VT, FPVT>>(iw);
    det.UpdateMode = ORB_PBYP_RATIO;
    ratios[iw] = det.curRatio = ratios_local[iw];
  }
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::evaluateRatios(const VirtualParticleSet& VP, std::vector<Value>& ratios)
{
  ScopedTimer local_timer(SPOVTimer);
  Vector<ValueType> inv_row(psiMinv_[VP.refPtcl - FirstIndex], psiV.size());
  phi_.evaluateDetRatios(VP, psiV_host_view, inv_row, ratios);
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::evaluateSpinorRatios(
    const VirtualParticleSet& VP,
    const std::pair<ValueVector, ValueVector>& spinor_multiplier,
    std::vector<Value>& ratios)
{
  ScopedTimer local_timer(SPOVTimer);
  Vector<ValueType> inv_row(psiMinv_[VP.refPtcl - FirstIndex], psiV.size());
  phi_.evaluateDetSpinorRatios(VP, psiV_host_view, spinor_multiplier, inv_row, ratios);
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::mw_evaluateRatios(
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
    const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
    std::vector<std::vector<Value>>& ratios) const
{
  assert(this == &wfc_list.getLeader());
  const size_t nw = wfc_list.size();

  RefVectorWithLeader<SPOSet> phi_list(phi_);
  RefVector<Vector<Value>> psiV_list;
  std::vector<const Value*> invRow_ptr_list;
  phi_list.reserve(nw);
  psiV_list.reserve(nw);
  invRow_ptr_list.reserve(nw);

  {
    ScopedTimer local_timer(RatioTimer);
    for (size_t iw = 0; iw < nw; iw++)
    {
      auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<PL, VT, FPVT>>(iw);
      const VirtualParticleSet& vp(vp_list[iw]);
      const int WorkingIndex = vp.refPtcl - FirstIndex;
      // build lists
      phi_list.push_back(det.phi_);
      psiV_list.push_back(det.psiV_host_view);
      if (phi_.isOMPoffload())
        invRow_ptr_list.push_back(det.psiMinv_.device_data() + WorkingIndex * psiMinv_.cols());
      else
        invRow_ptr_list.push_back(det.psiMinv_[WorkingIndex]);
    }
  }

  {
    ScopedTimer local_timer(SPOVTimer);
    phi_.mw_evaluateDetRatios(phi_list, vp_list, psiV_list, invRow_ptr_list, ratios);
  }
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::mw_evaluateSpinorRatios(
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
    const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
    const RefVector<std::pair<ValueVector, ValueVector>>& spinor_multiplier_list,
    std::vector<std::vector<Value>>& ratios) const
{
  assert(this == &wfc_list.getLeader());
  const size_t nw = wfc_list.size();

  RefVectorWithLeader<SPOSet> phi_list(phi_);
  RefVector<Vector<Value>> psiV_list;
  std::vector<const Value*> invRow_ptr_list;
  phi_list.reserve(nw);
  psiV_list.reserve(nw);
  invRow_ptr_list.reserve(nw);

  {
    ScopedTimer local_timer(RatioTimer);
    for (size_t iw = 0; iw < nw; iw++)
    {
      auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<PL, VT, FPVT>>(iw);
      const VirtualParticleSet& vp(vp_list[iw]);
      const int WorkingIndex = vp.refPtcl - FirstIndex;
      // build lists
      phi_list.push_back(det.phi_);
      psiV_list.push_back(det.psiV_host_view);
      if (phi_.isOMPoffload())
        invRow_ptr_list.push_back(det.psiMinv_.device_data() + WorkingIndex * psiMinv_.cols());
      else
        invRow_ptr_list.push_back(det.psiMinv_[WorkingIndex]);
    }
  }

  {
    ScopedTimer local_timer(SPOVTimer);
    phi_.mw_evaluateDetSpinorRatios(phi_list, vp_list, psiV_list, spinor_multiplier_list, invRow_ptr_list, ratios);
  }
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::evaluateDerivRatios(const VirtualParticleSet& VP,
                                                                const opt_variables_type& optvars,
                                                                std::vector<ValueType>& ratios,
                                                                Matrix<ValueType>& dratios)
{
  const int WorkingIndex = VP.refPtcl - FirstIndex;
  assert(WorkingIndex >= 0);
  std::copy_n(psiMinv_[WorkingIndex], d2psiV.size(), d2psiV.data());
  phi_.evaluateDerivRatios(VP, optvars, psiV_host_view, d2psiV_host_view, ratios, dratios, FirstIndex, LastIndex);
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::evaluateSpinorDerivRatios(
    const VirtualParticleSet& VP,
    const std::pair<ValueVector, ValueVector>& spinor_multiplier,
    const opt_variables_type& optvars,
    std::vector<ValueType>& ratios,
    Matrix<ValueType>& dratios)
{
  const int WorkingIndex = VP.refPtcl - FirstIndex;
  assert(WorkingIndex >= 0);
  std::copy_n(psiMinv_[WorkingIndex], d2psiV.size(), d2psiV.data());
  phi_.evaluateSpinorDerivRatios(VP, spinor_multiplier, optvars, psiV_host_view, d2psiV_host_view, ratios, dratios,
                                 FirstIndex, LastIndex);
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::evaluateRatiosAlltoOne(ParticleSet& P, std::vector<Value>& ratios)
{
  {
    ScopedTimer local_timer(SPOVTimer);
    phi_.evaluateValue(P, -1, psiV_host_view);
  }
  for (int i = 0; i < psiMinv_.rows(); i++)
    ratios[FirstIndex + i] = simd::dot(psiMinv_[i], psiV.data(), NumOrbitals);
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::resizeScratchObjectsForIonDerivs()
{
  grad_source_psiM.resize(NumPtcls, NumOrbitals);
  grad_lapl_source_psiM.resize(NumPtcls, NumOrbitals);
  grad_grad_source_psiM.resize(NumPtcls, NumOrbitals);
  phi_alpha_Minv.resize(NumPtcls, NumOrbitals);
  grad_phi_Minv.resize(NumPtcls, NumOrbitals);
  lapl_phi_Minv.resize(NumPtcls, NumOrbitals);
  grad_phi_alpha_Minv.resize(NumPtcls, NumOrbitals);
}

template<PlatformKind PL, typename VT, typename FPVT>
typename DiracDeterminantBatched<PL, VT, FPVT>::Grad DiracDeterminantBatched<PL, VT, FPVT>::evalGradSource(
    ParticleSet& P,
    ParticleSet& source,
    int iat)
{
  Grad g(0.0);
  if (phi_.hasIonDerivs())
  {
    resizeScratchObjectsForIonDerivs();
    phi_.evaluateGradSource(P, FirstIndex, LastIndex, source, iat, grad_source_psiM);
    // psiMinv columns have padding but grad_source_psiM ones don't
    for (int i = 0; i < psiMinv_.rows(); i++)
      g += simd::dot(psiMinv_[i], grad_source_psiM[i], NumOrbitals);
  }

  return g;
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::evaluateHessian(ParticleSet& P, HessVector& grad_grad_psi)
{
  // Hessian is not often used, so only resize/allocate if used
  grad_grad_source_psiM.resize(psiMinv_.rows(), NumOrbitals);
  //IM A HACK.  Assumes evaluateLog has already been executed.
  Matrix<Value> psiM_temp_host(psiM_temp.data(), psiM_temp.rows(), psiM_temp.cols());
  phi_.evaluate_notranspose(P, FirstIndex, LastIndex, psiM_temp_host, dpsiM, grad_grad_source_psiM);
  invertPsiM(psiM_temp, psiMinv_);

  phi_alpha_Minv      = 0.0;
  grad_phi_Minv       = 0.0;
  lapl_phi_Minv       = 0.0;
  grad_phi_alpha_Minv = 0.0;
  //grad_grad_psi.resize(NumPtcls);

  for (int i = 0, iat = FirstIndex; i < NumPtcls; i++, iat++)
  {
    Grad rv = simd::dot(psiMinv_[i], dpsiM[i], NumOrbitals);
    //  HessType hess_tmp=simd::dot(psiMinv_[i],grad_grad_source_psiM[i],NumOrbitals);
    Hess hess_tmp;
    hess_tmp           = 0.0;
    hess_tmp           = simd::dot(psiMinv_[i], grad_grad_source_psiM[i], NumOrbitals);
    grad_grad_psi[iat] = hess_tmp - outerProduct(rv, rv);
  }
}

template<PlatformKind PL, typename VT, typename FPVT>
typename DiracDeterminantBatched<PL, VT, FPVT>::Grad DiracDeterminantBatched<PL, VT, FPVT>::evalGradSource(
    ParticleSet& P,
    ParticleSet& source,
    int iat,
    TinyVector<ParticleSet::ParticleGradient, OHMMS_DIM>& grad_grad,
    TinyVector<ParticleSet::ParticleLaplacian, OHMMS_DIM>& lapl_grad)
{
  Grad gradPsi(0.0);
  if (phi_.hasIonDerivs())
  {
    resizeScratchObjectsForIonDerivs();
    phi_.evaluateGradSource(P, FirstIndex, LastIndex, source, iat, grad_source_psiM, grad_grad_source_psiM,
                            grad_lapl_source_psiM);

    // Compute matrices
    phi_alpha_Minv      = 0.0;
    grad_phi_Minv       = 0.0;
    lapl_phi_Minv       = 0.0;
    grad_phi_alpha_Minv = 0.0;
    for (int i = 0; i < NumPtcls; i++)
      for (int j = 0; j < NumOrbitals; j++)
      {
        lapl_phi_Minv(i, j) = 0.0;
        for (int k = 0; k < NumOrbitals; k++)
          lapl_phi_Minv(i, j) += d2psiM(i, k) * psiMinv_(j, k);
      }
    for (int dim = 0; dim < OHMMS_DIM; dim++)
    {
      for (int i = 0; i < NumPtcls; i++)
        for (int j = 0; j < NumOrbitals; j++)
        {
          for (int k = 0; k < NumOrbitals; k++)
          {
            phi_alpha_Minv(i, j)[dim] += grad_source_psiM(i, k)[dim] * psiMinv_(j, k);
            grad_phi_Minv(i, j)[dim] += dpsiM(i, k)[dim] * psiMinv_(j, k);
            for (int dim_el = 0; dim_el < OHMMS_DIM; dim_el++)
              grad_phi_alpha_Minv(i, j)(dim, dim_el) += grad_grad_source_psiM(i, k)(dim, dim_el) * psiMinv_(j, k);
          }
        }
    }
    for (int i = 0, iel = FirstIndex; i < NumPtcls; i++, iel++)
    {
      Hess dval(0.0);
      Grad d2val(0.0);
      for (int dim = 0; dim < OHMMS_DIM; dim++)
        for (int dim_el = 0; dim_el < OHMMS_DIM; dim_el++)
          dval(dim, dim_el) = grad_phi_alpha_Minv(i, i)(dim, dim_el);
      for (int j = 0; j < NumOrbitals; j++)
      {
        gradPsi += grad_source_psiM(i, j) * psiMinv_(i, j);
        for (int dim = 0; dim < OHMMS_DIM; dim++)
          for (int k = 0; k < OHMMS_DIM; k++)
            dval(dim, k) -= phi_alpha_Minv(j, i)[dim] * grad_phi_Minv(i, j)[k];
      }
      for (int dim = 0; dim < OHMMS_DIM; dim++)
      {
        for (int k = 0; k < OHMMS_DIM; k++)
          grad_grad[dim][iel][k] += dval(dim, k);
        for (int j = 0; j < NumOrbitals; j++)
        {
          // First term, eq 9
          lapl_grad[dim][iel] += grad_lapl_source_psiM(i, j)[dim] * psiMinv_(i, j);
          // Second term, eq 9
          if (j == i)
            for (int dim_el = 0; dim_el < OHMMS_DIM; dim_el++)
              lapl_grad[dim][iel] -= (Real)2.0 * grad_phi_alpha_Minv(j, i)(dim, dim_el) * grad_phi_Minv(i, j)[dim_el];
          // Third term, eq 9
          // First term, eq 10
          lapl_grad[dim][iel] -= phi_alpha_Minv(j, i)[dim] * lapl_phi_Minv(i, j);
          // Second term, eq 11
          for (int dim_el = 0; dim_el < OHMMS_DIM; dim_el++)
            lapl_grad[dim][iel] +=
                (Real)2.0 * phi_alpha_Minv(j, i)[dim] * grad_phi_Minv(i, i)[dim_el] * grad_phi_Minv(i, j)[dim_el];
        }
      }
    }
  }
  return gradPsi;
}

/** Calculate the log value of the Dirac determinant for particles
 *@param P input configuration containing N particles
 *@param G a vector containing N gradients
 *@param L a vector containing N laplacians
 *\return the complex value of the determinant
 *
 *\f$ (first,first+nel). \f$  Add the gradient and laplacian
 *contribution of the determinant to G(radient) and L(aplacian)
 *for local energy calculations.
 */
template<PlatformKind PL, typename VT, typename FPVT>
typename DiracDeterminantBatched<PL, VT, FPVT>::LogValue DiracDeterminantBatched<PL, VT, FPVT>::evaluateLog(
    const ParticleSet& P,
    ParticleSet::ParticleGradient& G,
    ParticleSet::ParticleLaplacian& L)
{
  recompute(P);
  computeGL(G, L);
  return log_value_;
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::mw_evaluateLog(
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
    const RefVectorWithLeader<ParticleSet>& p_list,
    const RefVector<ParticleSet::ParticleGradient>& G_list,
    const RefVector<ParticleSet::ParticleLaplacian>& L_list) const
{
  assert(this == &wfc_list.getLeader());
  const std::vector<bool> recompute_all(wfc_list.size(), true);
  mw_recompute(wfc_list, p_list, recompute_all);

  for (int iw = 0; iw < wfc_list.size(); iw++)
  {
    auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<PL, VT, FPVT>>(iw);
    det.computeGL(G_list[iw], L_list[iw]);
  }
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::recompute(const ParticleSet& P)
{
  {
    ScopedTimer spo_timer(SPOVGLTimer);

    UpdateMode = ORB_WALKER;
    phi_.evaluate_notranspose(P, FirstIndex, LastIndex, psiM_host, dpsiM, d2psiM);
    auto* psiM_vgl_ptr = psiM_vgl.data();
    // transfer host to device, total size 4, g(3) + l(1)
    PRAGMA_OFFLOAD("omp target update to(psiM_vgl_ptr[psiM_vgl.capacity():psiM_vgl.capacity()*4])")
  }
  invertPsiM(psiM_temp, psiMinv_);
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::mw_recompute(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                         const RefVectorWithLeader<ParticleSet>& p_list,
                                                         const std::vector<bool>& recompute) const
{
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<PL, VT, FPVT>>();
  const auto nw    = wfc_list.size();

  RefVectorWithLeader<WaveFunctionComponent> wfc_filtered_list(wfc_list.getLeader());
  RefVectorWithLeader<ParticleSet> p_filtered_list(p_list.getLeader());
  RefVectorWithLeader<SPOSet> phi_list(wfc_leader.phi_);
  std::vector<Matrix<Value>> psiM_host_views;
  RefVector<const DualMatrix<Value>> psiM_temp_list;
  RefVector<Matrix<Value>> psiM_host_list;
  RefVector<Matrix<Grad>> dpsiM_list;
  RefVector<Matrix<Value>> d2psiM_list;
  RefVector<DualMatrix<Value>> psiMinv_list;

  wfc_filtered_list.reserve(nw);
  p_filtered_list.reserve(nw);
  phi_list.reserve(nw);
  psiM_host_views.reserve(nw);
  psiM_temp_list.reserve(nw);
  dpsiM_list.reserve(nw);
  d2psiM_list.reserve(nw);
  psiMinv_list.reserve(nw);

  for (int iw = 0; iw < nw; iw++)
    if (recompute[iw])
    {
      wfc_filtered_list.push_back(wfc_list[iw]);
      p_filtered_list.push_back(p_list[iw]);

      auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<PL, VT, FPVT>>(iw);
      phi_list.push_back(det.phi_);
      psiM_temp_list.push_back(det.psiM_temp);
      psiM_host_list.push_back(det.psiM_host);
      dpsiM_list.push_back(det.dpsiM);
      d2psiM_list.push_back(det.d2psiM);
      psiMinv_list.push_back(det.psiMinv_);
    }

  if (!wfc_filtered_list.size())
    return;

  {
    ScopedTimer spo_timer(wfc_leader.SPOVGLTimer);
    // I think through the magic of OMPtarget psiM_host_list actually results in psiM_temp being updated
    // on the device. For dspiM_list, d2psiM_list I think they are calculated on CPU and this is not true
    // This is the reason for the strange look omp target update to below.
    wfc_leader.phi_.mw_evaluate_notranspose(phi_list, p_filtered_list, wfc_leader.FirstIndex, wfc_leader.LastIndex,
                                            psiM_host_list, dpsiM_list, d2psiM_list);
  }

  mw_invertPsiM(wfc_filtered_list, psiM_temp_list, psiMinv_list);

  { // transfer dpsiM, d2psiM, psiMinv to device
    ScopedTimer d2h(H2DTimer);

    RefVector<DualVGLVector> psiM_vgl_list;
    psiM_vgl_list.reserve(nw);
    for (int iw = 0; iw < wfc_filtered_list.size(); iw++)
    {
      auto& det = wfc_filtered_list.getCastedElement<DiracDeterminantBatched<PL, VT, FPVT>>(iw);
      psiM_vgl_list.push_back(det.psiM_vgl);
      det.UpdateMode = ORB_WALKER;
    }

    auto& queue = wfc_leader.mw_res_handle_.getResource().engine_rsc.queue;
    for (DualVGLVector& psiM_vgl : psiM_vgl_list)
    {
      const size_t stride = psiM_vgl.capacity();
      // transfer host to device, total size 4, g(3) + l(1), skipping v
      queue.enqueueH2D(psiM_vgl, stride * 4, stride);
    }
    queue.sync();
  }
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::evaluateDerivatives(ParticleSet& P,
                                                                const opt_variables_type& active,
                                                                Vector<Value>& dlogpsi,
                                                                Vector<Value>& dhpsioverpsi)
{
  phi_.evaluateDerivatives(P, active, dlogpsi, dhpsioverpsi, FirstIndex, LastIndex);
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::evaluateDerivativesWF(ParticleSet& P,
                                                                  const opt_variables_type& active,
                                                                  Vector<ValueType>& dlogpsi)
{
  phi_.evaluateDerivativesWF(P, active, dlogpsi, FirstIndex, LastIndex);
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::registerTWFFastDerivWrapper(const ParticleSet& P,
                                                                        TWFFastDerivWrapper& twf) const
{
  twf.addGroup(P, P.getGroupID(FirstIndex), &phi_);
}

template<PlatformKind PL, typename VT, typename FPVT>
std::unique_ptr<DiracDeterminantBase> DiracDeterminantBatched<PL, VT, FPVT>::makeCopy(SPOSet& phi) const
{
  return std::make_unique<DiracDeterminantBatched<PL, VT, FPVT>>(phi, FirstIndex, LastIndex, ndelay_,
                                                                 matrix_inverter_kind_);
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::createResource(ResourceCollection& collection) const
{
  collection.addResource(std::make_unique<DiracDeterminantBatchedMultiWalkerResource>());
  if (matrix_inverter_kind_ == DetMatInvertor::ACCEL)
    collection.addResource(std::make_unique<typename DetInverterSelector<PL, FPVT, VT>::Inverter>());
  else
    collection.addResource(std::make_unique<DiracMatrixInverterOMPTarget<FPVT, VT>>());
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::acquireResource(
    ResourceCollection& collection,
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const
{
  auto& wfc_leader          = wfc_list.getCastedLeader<DiracDeterminantBatched<PL, VT, FPVT>>();
  wfc_leader.mw_res_handle_ = collection.lendResource<DiracDeterminantBatchedMultiWalkerResource>();
  auto& mw_res              = wfc_leader.mw_res_handle_.getResource();
  mw_res.psiMinv_refs.reserve(wfc_list.size());
  for (WaveFunctionComponent& wfc : wfc_list)
  {
    auto& det = static_cast<DiracDeterminantBatched<PL, VT, FPVT>&>(wfc);
    mw_res.psiMinv_refs.push_back(det.psiMinv_);
  }
  wfc_leader.accel_inverter_ = collection.lendResource<DiracMatrixInverter<FPVT, VT>>();
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminantBatched<PL, VT, FPVT>::releaseResource(
    ResourceCollection& collection,
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const
{
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<PL, VT, FPVT>>();
  wfc_leader.mw_res_handle_.getResource().psiMinv_refs.clear();
  collection.takebackResource(wfc_leader.mw_res_handle_);
  collection.takebackResource(wfc_leader.accel_inverter_);
}

template class DiracDeterminantBatched<PlatformKind::OMPTARGET, QMCTraits::ValueType, QMCTraits::QTFull::ValueType>;
#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
template class DiracDeterminantBatched<PlatformKind::CUDA, QMCTraits::ValueType, QMCTraits::QTFull::ValueType>;
#endif
#if defined(ENABLE_SYCL) && defined(ENABLE_OFFLOAD)
template class DiracDeterminantBatched<PlatformKind::SYCL, QMCTraits::ValueType, QMCTraits::QTFull::ValueType>;
#endif

} // namespace qmcplusplus
