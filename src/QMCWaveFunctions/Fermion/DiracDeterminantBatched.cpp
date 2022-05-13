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
#include "Numerics/DeterminantOperators.h"
#include "CPU/BLAS.hpp"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Numerics/MatrixOperators.h"
#include "QMCWaveFunctions/TWFFastDerivWrapper.h"
#include "CPU/SIMD/simd.hpp"
#include <cassert>

namespace qmcplusplus
{
/** constructor
 *@param spos the single-particle orbital set
 *@param first index of the first particle
 */
template<typename DET_ENGINE>
DiracDeterminantBatched<DET_ENGINE>::DiracDeterminantBatched(std::unique_ptr<SPOSet>&& spos,
                                                             int first,
                                                             int last,
                                                             int ndelay,
                                                             DetMatInvertor matrix_inverter_kind)
    : DiracDeterminantBase("DiracDeterminantBatched", std::move(spos), first, last),
      ndelay_(ndelay),
      matrix_inverter_kind_(matrix_inverter_kind),
      D2HTimer(*timer_manager.createTimer("DiracDeterminantBatched::D2H", timer_level_fine)),
      H2DTimer(*timer_manager.createTimer("DiracDeterminantBatched::H2D", timer_level_fine))
{
  static_assert(std::is_same<SPOSet::ValueType, typename DET_ENGINE::Value>::value);
  resize(NumPtcls, NumPtcls);
  if (Optimizable)
    Phi->buildOptVariables(NumPtcls);
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::invertPsiM(const DualMatrix<Value>& psiM, DualMatrix<Value>& psiMinv)
{
  ScopedTimer inverse_timer(InverseTimer);
  host_inverter_.invert_transpose(psiM, psiMinv, log_value_);
  psiMinv.updateTo();

#ifndef NDEBUG
  // This is easily breakable in that it assumes this function gets psiMinv == det_engine_.psiMinv_
  auto& engine_psiMinv = det_engine_.get_ref_psiMinv();
  dummy_vmt.attachReference(engine_psiMinv.data(), engine_psiMinv.rows(), engine_psiMinv.cols());
#endif
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::mw_invertPsiM(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                        const RefVector<const DualMatrix<Value>>& logdetT_list,
                                                        const RefVector<DualMatrix<Value>>& a_inv_list)
{
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<DET_ENGINE>>();
  ScopedTimer inverse_timer(wfc_leader.InverseTimer);
  const auto nw = wfc_list.size();

  if (wfc_leader.matrix_inverter_kind_ == DetMatInvertor::ACCEL)
  {
    RefVectorWithLeader<DET_ENGINE> engine_list(wfc_leader.det_engine_);
    engine_list.reserve(nw);

    wfc_leader.mw_res_->log_values.resize(nw);

    for (int iw = 0; iw < nw; iw++)
    {
      auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE>>(iw);
      engine_list.push_back(det.get_det_engine());
      wfc_leader.mw_res_->log_values[iw] = {0.0, 0.0};
    }

    wfc_leader.accel_inverter_->mw_invertTranspose(wfc_leader.det_engine_.getLAhandles(), logdetT_list, a_inv_list,
                                                   wfc_leader.mw_res_->log_values);

    for (int iw = 0; iw < nw; ++iw)
    {
      auto& det      = wfc_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE>>(iw);
      det.log_value_ = wfc_leader.mw_res_->log_values[iw];
    }
  }
  else
  {
    for (int iw = 0; iw < nw; iw++)
    {
      auto& det     = wfc_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE>>(iw);
      auto& psiMinv = a_inv_list[iw].get();
      det.host_inverter_.invert_transpose(logdetT_list[iw].get(), psiMinv, det.log_value_);
      psiMinv.updateTo();
    }
  }

#ifndef NDEBUG
  for (int iw = 0; iw < nw; ++iw)
  {
    auto& det            = wfc_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE>>(iw);
    auto& engine_psiMinv = det.get_det_engine().get_ref_psiMinv();
    det.dummy_vmt.attachReference(engine_psiMinv.data(), engine_psiMinv.rows(), engine_psiMinv.cols());
  }
#endif
}

///reset the size: with the number of particles and number of orbtials
template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::resize(int nel, int morb)
{
  int norb = morb;
  if (norb <= 0)
    norb = nel; // for morb == -1 (default)
  psiM_vgl.resize(nel * norb);
  // attach pointers VGL
  psiM_temp.attachReference(psiM_vgl, psiM_vgl.data(0), nel, norb);
  psiM_host.attachReference(psiM_vgl.data(0), nel, norb);
  dpsiM.attachReference(reinterpret_cast<Grad*>(psiM_vgl.data(1)), nel, norb);
  d2psiM.attachReference(psiM_vgl.data(4), nel, norb);

  det_engine_.resize(norb, ndelay_);

  psiV.resize(NumOrbitals);
  psiV_host_view.attachReference(psiV.data(), NumOrbitals);
  dpsiV.resize(NumOrbitals);
  dpsiV_host_view.attachReference(dpsiV.data(), NumOrbitals);
  d2psiV.resize(NumOrbitals);
  d2psiV_host_view.attachReference(d2psiV.data(), NumOrbitals);
}

template<typename DET_ENGINE>
typename DiracDeterminantBatched<DET_ENGINE>::Grad DiracDeterminantBatched<DET_ENGINE>::evalGrad(ParticleSet& P,
                                                                                                 int iat)
{
  ScopedTimer local_timer(RatioTimer);
  const int WorkingIndex = iat - FirstIndex;
  Grad g                 = simd::dot(det_engine_.get_psiMinv()[WorkingIndex], dpsiM[WorkingIndex], NumOrbitals);
  assert(checkG(g));
  return g;
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::mw_evalGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                      const RefVectorWithLeader<ParticleSet>& p_list,
                                                      int iat,
                                                      std::vector<Grad>& grad_now) const
{
  assert(this == &wfc_list.getLeader());
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<DET_ENGINE>>();
  ScopedTimer local_timer(RatioTimer);

  const int nw = wfc_list.size();
  std::vector<const Value*> dpsiM_row_list(nw, nullptr);
  RefVectorWithLeader<DET_ENGINE> engine_list(wfc_leader.det_engine_);
  engine_list.reserve(nw);

  const int WorkingIndex = iat - FirstIndex;
  for (int iw = 0; iw < nw; iw++)
  {
    auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE>>(iw);
    // capacity is the size of each vector in the VGL so this advances us to the g then makes
    // an offset into the gradients
    dpsiM_row_list[iw] = det.psiM_vgl.device_data() + psiM_vgl.capacity() + NumOrbitals * WorkingIndex * DIM;
    engine_list.push_back(det.det_engine_);
  }

  DET_ENGINE::mw_evalGrad(engine_list, dpsiM_row_list, WorkingIndex, grad_now);

#ifndef NDEBUG
  for (int iw = 0; iw < nw; iw++)
    checkG(grad_now[iw]);
#endif
}

template<typename DET_ENGINE>
typename DiracDeterminantBatched<DET_ENGINE>::PsiValue DiracDeterminantBatched<DET_ENGINE>::ratioGrad(ParticleSet& P,
                                                                                                      int iat,
                                                                                                      Grad& grad_iat)
{
  UpdateMode = ORB_PBYP_PARTIAL;

  {
    ScopedTimer local_timer(SPOVGLTimer);
    Phi->evaluateVGL(P, iat, psiV_host_view, dpsiV_host_view, d2psiV_host_view);
  }

  {
    ScopedTimer local_timer(RatioTimer);
    auto& psiMinv          = det_engine_.get_psiMinv();
    const int WorkingIndex = iat - FirstIndex;
    curRatio               = simd::dot(psiMinv[WorkingIndex], psiV.data(), NumOrbitals);
    grad_iat += static_cast<Value>(static_cast<PsiValue>(1.0) / curRatio) *
        simd::dot(psiMinv[WorkingIndex], dpsiV.data(), NumOrbitals);
  }
  return curRatio;
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::mw_ratioGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                       const RefVectorWithLeader<ParticleSet>& p_list,
                                                       int iat,
                                                       std::vector<PsiValue>& ratios,
                                                       std::vector<Grad>& grad_new) const
{
  assert(this == &wfc_list.getLeader());
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<DET_ENGINE>>();
  wfc_leader.guardMultiWalkerRes();
  auto& mw_res         = *wfc_leader.mw_res_;
  auto& phi_vgl_v      = mw_res.phi_vgl_v;
  auto& ratios_local   = mw_res.ratios_local;
  auto& grad_new_local = mw_res.grad_new_local;
  {
    ScopedTimer local_timer(SPOVGLTimer);
    RefVectorWithLeader<SPOSet> phi_list(*Phi);
    phi_list.reserve(wfc_list.size());
    RefVectorWithLeader<DET_ENGINE> engine_list(wfc_leader.det_engine_);
    engine_list.reserve(wfc_list.size());
    const int WorkingIndex = iat - FirstIndex;
    for (int iw = 0; iw < wfc_list.size(); iw++)
    {
      auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE>>(iw);
      phi_list.push_back(*det.Phi);
      engine_list.push_back(det.det_engine_);
    }

    auto psiMinv_row_dev_ptr_list = DET_ENGINE::mw_getInvRow(engine_list, WorkingIndex, !Phi->isOMPoffload());

    phi_vgl_v.resize(NumOrbitals * wfc_list.size());
    ratios_local.resize(wfc_list.size());
    grad_new_local.resize(wfc_list.size());

    VectorSoaContainer<Value, DIM + 2> phi_vgl_v_view(phi_vgl_v.data(), NumOrbitals * wfc_list.size(),
                                                      phi_vgl_v.capacity());
    wfc_leader.Phi->mw_evaluateVGLandDetRatioGrads(phi_list, p_list, iat, psiMinv_row_dev_ptr_list, phi_vgl_v_view,
                                                   ratios_local, grad_new_local);
  }

  wfc_leader.UpdateMode = ORB_PBYP_PARTIAL;
  for (int iw = 0; iw < wfc_list.size(); iw++)
  {
    auto& det      = wfc_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE>>(iw);
    det.UpdateMode = ORB_PBYP_PARTIAL;
    ratios[iw] = det.curRatio = ratios_local[iw];
    grad_new[iw] += grad_new_local[iw];
  }
}


/** move was accepted, update the real container
*/
template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::acceptMove(ParticleSet& P, int iat, bool safe_to_delay)
{
  const int WorkingIndex = iat - FirstIndex;
  log_value_ += convertValueToLog(curRatio);
  {
    ScopedTimer local_timer(UpdateTimer);
    det_engine_.updateRow(WorkingIndex, psiV, curRatio);
    if (UpdateMode == ORB_PBYP_PARTIAL)
    {
      simd::copy(dpsiM[WorkingIndex], dpsiV.data(), NumOrbitals);
      simd::copy(d2psiM[WorkingIndex], d2psiV.data(), NumOrbitals);
    }
  }
  curRatio = 1.0;
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::mw_accept_rejectMove(
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
    const RefVectorWithLeader<ParticleSet>& p_list,
    int iat,
    const std::vector<bool>& isAccepted,
    bool safe_to_delay) const
{
  assert(this == &wfc_list.getLeader());
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<DET_ENGINE>>();
  wfc_leader.guardMultiWalkerRes();
  auto& mw_res       = *wfc_leader.mw_res_;
  auto& phi_vgl_v    = mw_res.phi_vgl_v;
  auto& ratios_local = mw_res.ratios_local;

  ScopedTimer update(UpdateTimer);

  const int nw = wfc_list.size();
  int count    = 0;
  for (int iw = 0; iw < nw; iw++)
    if (isAccepted[iw])
      count++;
  const int n_accepted = count;

  RefVectorWithLeader<DET_ENGINE> engine_list(wfc_leader.det_engine_);
  engine_list.reserve(nw);
  std::vector<Value*> psiM_g_dev_ptr_list(n_accepted, nullptr);
  std::vector<Value*> psiM_l_dev_ptr_list(n_accepted, nullptr);

  const int WorkingIndex = iat - FirstIndex;
  for (int iw = 0, count = 0; iw < nw; iw++)
  {
    // This can be auto but some debuggers can't figure the type out.
    DiracDeterminantBatched<DET_ENGINE>& det = wfc_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE>>(iw);
    engine_list.push_back(det.det_engine_);
    if (isAccepted[iw])
    {
      psiM_g_dev_ptr_list[count] = det.psiM_vgl.device_data() + psiM_vgl.capacity() + NumOrbitals * WorkingIndex * DIM;
      psiM_l_dev_ptr_list[count] = det.psiM_vgl.device_data() + psiM_vgl.capacity() * 4 + NumOrbitals * WorkingIndex;
      det.log_value_ += convertValueToLog(det.curRatio);
      count++;
    }
    det.curRatio = 1.0;
  }

  if (!Phi->isOMPoffload() && n_accepted > 0)
  {
    auto* phi_vgl_v_ptr = phi_vgl_v.data();
    // transfer host to device, total size 5, v(1) + g(3) + l(1)
    PRAGMA_OFFLOAD("omp target update to(phi_vgl_v_ptr[:phi_vgl_v.capacity()*5])")
  }

  DET_ENGINE::mw_accept_rejectRow(engine_list, WorkingIndex, psiM_g_dev_ptr_list, psiM_l_dev_ptr_list, isAccepted,
                                  phi_vgl_v.device_data(), phi_vgl_v.capacity(), ratios_local);

  if (!safe_to_delay)
    DET_ENGINE::mw_updateInvMat(engine_list);
}

/** move was rejected. copy the real container to the temporary to move on
*/
template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::restore(int iat)
{
  curRatio = 1.0;
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::completeUpdates()
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

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::mw_completeUpdates(
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const
{
  assert(this == &wfc_list.getLeader());
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<DET_ENGINE>>();
  const auto nw    = wfc_list.size();
  RefVectorWithLeader<DET_ENGINE> engine_list(wfc_leader.det_engine_);
  engine_list.reserve(nw);
  for (int iw = 0; iw < nw; iw++)
  {
    auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE>>(iw);
    engine_list.push_back(det.det_engine_);
  }

  {
    ScopedTimer update(UpdateTimer);
    DET_ENGINE::mw_updateInvMat(engine_list);
  }

  { // transfer dpsiM, d2psiM, psiMinv to host
    ScopedTimer d2h(D2HTimer);

    // this call also completes all the device copying of dpsiM, d2psiM before the target update
    DET_ENGINE::mw_transferAinv_D2H(engine_list);

    if (UpdateMode == ORB_PBYP_PARTIAL)
    {
      RefVector<DualVGLVector> psiM_vgl_list;
      psiM_vgl_list.reserve(nw);
      for (int iw = 0; iw < nw; iw++)
      {
        auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE>>(iw);
        psiM_vgl_list.push_back(det.psiM_vgl);
      }

      // transfer device to host, total size 4, g(3) + l(1), skipping v
      DET_ENGINE::mw_transferVGL_D2H(wfc_leader.det_engine_, psiM_vgl_list, 1, 4);
    }
  }
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::computeGL(ParticleSet::ParticleGradient& G,
                                                    ParticleSet::ParticleLaplacian& L) const
{
  auto& psiMinv = det_engine_.get_psiMinv();
  for (size_t i = 0, iat = FirstIndex; i < NumPtcls; ++i, ++iat)
  {
    Grad rv   = simd::dot(psiMinv[i], dpsiM[i], NumOrbitals);
    Value lap = simd::dot(psiMinv[i], d2psiM[i], NumOrbitals);
    G[iat] += rv;
    L[iat] += lap - dot(rv, rv);
  }
}

template<typename DET_ENGINE>
typename DiracDeterminantBatched<DET_ENGINE>::LogValue DiracDeterminantBatched<DET_ENGINE>::evaluateGL(
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
      Phi->evaluate_notranspose(P, FirstIndex, LastIndex, psiM_host, dpsiM, d2psiM);
    }
    UpdateMode = ORB_WALKER;
    computeGL(G, L);
  }
  return log_value_;
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::mw_evaluateGL(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                        const RefVectorWithLeader<ParticleSet>& p_list,
                                                        const RefVector<ParticleSet::ParticleGradient>& G_list,
                                                        const RefVector<ParticleSet::ParticleLaplacian>& L_list,
                                                        bool fromscratch) const
{
  assert(this == &wfc_list.getLeader());
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<DET_ENGINE>>();
  if (fromscratch)
    mw_evaluateLog(wfc_list, p_list, G_list, L_list);
  else
  {
    const auto nw = wfc_list.size();
    RefVectorWithLeader<DET_ENGINE> engine_list(wfc_leader.get_det_engine());
    engine_list.reserve(nw);

    if (UpdateMode == ORB_PBYP_RATIO)
    { //need to compute dpsiM and d2psiM. psiMinv is not touched!
      ScopedTimer spo_timer(SPOVGLTimer);

      RefVectorWithLeader<SPOSet> phi_list(*Phi);
      RefVector<Matrix<Value>> psiM_temp_list;
      RefVector<Matrix<Grad>> dpsiM_list;
      RefVector<Matrix<Value>> d2psiM_list;
      phi_list.reserve(wfc_list.size());
      psiM_temp_list.reserve(nw);
      dpsiM_list.reserve(nw);
      d2psiM_list.reserve(nw);

      for (int iw = 0; iw < nw; iw++)
      {
        auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE>>(iw);
        engine_list.push_back(det.det_engine_);
        phi_list.push_back(*det.Phi);
        psiM_temp_list.push_back(det.psiM_host);
        dpsiM_list.push_back(det.dpsiM);
        d2psiM_list.push_back(det.d2psiM);
      }

      Phi->mw_evaluate_notranspose(phi_list, p_list, FirstIndex, LastIndex, psiM_temp_list, dpsiM_list, d2psiM_list);
    }

    for (int iw = 0; iw < nw; iw++)
    {
      auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE>>(iw);

#ifndef NDEBUG
      GradMatrix dpsiM_from_device     = det.dpsiM;
      Matrix<Value> d2psiM_from_device = det.d2psiM;

      auto& my_psiM_vgl  = det.psiM_vgl;
      auto* psiM_vgl_ptr = my_psiM_vgl.data();
      // transfer device to host, total size 4, g(3) + l(1)
      PRAGMA_OFFLOAD("omp target update from(psiM_vgl_ptr[my_psiM_vgl.capacity():my_psiM_vgl.capacity()*4])")
      Matrix<Value> psiM_temp_host(det.psiM_temp.data(), det.psiM_temp.rows(), det.psiM_temp.cols());
      det.Phi->evaluate_notranspose(p_list[iw], FirstIndex, LastIndex, psiM_temp_host, det.dpsiM, det.d2psiM);

      assert(dpsiM_from_device == det.dpsiM);
      assert(d2psiM_from_device == det.d2psiM);
#endif

      det.UpdateMode = ORB_WALKER;
      det.computeGL(G_list[iw], L_list[iw]);
    }
  }
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::registerData(ParticleSet& P, WFBufferType& buf)
{
  auto& psiMinv = det_engine_.get_psiMinv();
  buf.add(psiMinv.first_address(), psiMinv.last_address());
  buf.add(dpsiM.first_address(), dpsiM.last_address());
  buf.add(d2psiM.first_address(), d2psiM.last_address());
  buf.add(log_value_);
}

template<typename DET_ENGINE>
typename DiracDeterminantBatched<DET_ENGINE>::LogValue DiracDeterminantBatched<DET_ENGINE>::updateBuffer(
    ParticleSet& P,
    WFBufferType& buf,
    bool fromscratch)
{
  evaluateGL(P, P.G, P.L, fromscratch);
  auto& psiMinv = det_engine_.get_psiMinv();
  ScopedTimer local_timer(BufferTimer);
  buf.put(psiMinv.first_address(), psiMinv.last_address());
  buf.put(dpsiM.first_address(), dpsiM.last_address());
  buf.put(d2psiM.first_address(), d2psiM.last_address());
  buf.put(log_value_);
  return log_value_;
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  ScopedTimer local_timer(BufferTimer);
  auto& psiMinv = det_engine_.get_ref_psiMinv();
  buf.get(psiMinv.first_address(), psiMinv.last_address());
  buf.get(dpsiM.first_address(), dpsiM.last_address());
  buf.get(d2psiM.first_address(), d2psiM.last_address());
  auto* psiMinv_ptr = psiMinv.data();
  PRAGMA_OFFLOAD("omp target update to(psiMinv_ptr[:psiMinv.size()])")
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
template<typename DET_ENGINE>
typename DiracDeterminantBatched<DET_ENGINE>::PsiValue DiracDeterminantBatched<DET_ENGINE>::ratio(ParticleSet& P,
                                                                                                  int iat)
{
  UpdateMode             = ORB_PBYP_RATIO;
  const int WorkingIndex = iat - FirstIndex;
  {
    ScopedTimer local_timer(SPOVTimer);
    Phi->evaluateValue(P, iat, psiV_host_view);
  }
  {
    auto& psiMinv = det_engine_.get_psiMinv();
    ScopedTimer local_timer(RatioTimer);
    curRatio = simd::dot(psiMinv[WorkingIndex], psiV.data(), NumOrbitals);
  }
  return curRatio;
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::mw_calcRatio(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                       const RefVectorWithLeader<ParticleSet>& p_list,
                                                       int iat,
                                                       std::vector<PsiValue>& ratios) const
{
  assert(this == &wfc_list.getLeader());
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<DET_ENGINE>>();
  wfc_leader.guardMultiWalkerRes();
  auto& mw_res         = *wfc_leader.mw_res_;
  auto& phi_vgl_v      = mw_res.phi_vgl_v;
  auto& ratios_local   = mw_res.ratios_local;
  auto& grad_new_local = mw_res.grad_new_local;

  {
    ScopedTimer local_timer(SPOVTimer);
    RefVectorWithLeader<SPOSet> phi_list(*Phi);
    phi_list.reserve(wfc_list.size());
    RefVectorWithLeader<DET_ENGINE> engine_list(wfc_leader.det_engine_);
    engine_list.reserve(wfc_list.size());

    const int WorkingIndex = iat - FirstIndex;
    for (int iw = 0; iw < wfc_list.size(); iw++)
    {
      auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE>>(iw);
      phi_list.push_back(*det.Phi);
      engine_list.push_back(det.det_engine_);
    }

    auto psiMinv_row_dev_ptr_list = DET_ENGINE::mw_getInvRow(engine_list, WorkingIndex, !Phi->isOMPoffload());

    phi_vgl_v.resize(NumOrbitals * wfc_list.size());
    ratios_local.resize(wfc_list.size());
    grad_new_local.resize(wfc_list.size());

    VectorSoaContainer<Value, DIM + 2> phi_vgl_v_view(phi_vgl_v.data(), NumOrbitals * wfc_list.size(),
                                                      phi_vgl_v.capacity());
    // calling Phi->mw_evaluateVGLandDetRatioGrads is a temporary workaround.
    // We may implement mw_evaluateVandDetRatio in the future.
    wfc_leader.Phi->mw_evaluateVGLandDetRatioGrads(phi_list, p_list, iat, psiMinv_row_dev_ptr_list, phi_vgl_v_view,
                                                   ratios_local, grad_new_local);
  }

  wfc_leader.UpdateMode = ORB_PBYP_RATIO;
  for (int iw = 0; iw < wfc_list.size(); iw++)
  {
    auto& det      = wfc_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE>>(iw);
    det.UpdateMode = ORB_PBYP_RATIO;
    ratios[iw] = det.curRatio = ratios_local[iw];
  }
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::evaluateRatios(const VirtualParticleSet& VP, std::vector<Value>& ratios)
{
  {
    ScopedTimer local_timer(RatioTimer);
    const int WorkingIndex = VP.refPtcl - FirstIndex;
    std::copy_n(det_engine_.get_psiMinv()[WorkingIndex], d2psiV.size(), d2psiV.data());
  }
  {
    ScopedTimer local_timer(SPOVTimer);
    Phi->evaluateDetRatios(VP, psiV_host_view, d2psiV_host_view, ratios);
  }
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::mw_evaluateRatios(
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
    const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
    std::vector<std::vector<Value>>& ratios) const
{
  assert(this == &wfc_list.getLeader());
  const size_t nw = wfc_list.size();

  RefVectorWithLeader<SPOSet> phi_list(*Phi);
  RefVector<Vector<Value>> psiV_list;
  std::vector<const Value*> invRow_ptr_list;
  phi_list.reserve(nw);
  psiV_list.reserve(nw);
  invRow_ptr_list.reserve(nw);

  {
    ScopedTimer local_timer(RatioTimer);
    for (size_t iw = 0; iw < nw; iw++)
    {
      auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE>>(iw);
      const VirtualParticleSet& vp(vp_list[iw]);
      const int WorkingIndex = vp.refPtcl - FirstIndex;
      // build lists
      phi_list.push_back(*det.Phi);
      psiV_list.push_back(det.psiV_host_view);
      if (Phi->isOMPoffload())
        invRow_ptr_list.push_back(det.det_engine_.getRow_psiMinv_offload(WorkingIndex));
      else
        invRow_ptr_list.push_back(det.det_engine_.get_psiMinv()[WorkingIndex]);
    }
  }

  {
    ScopedTimer local_timer(SPOVTimer);
    Phi->mw_evaluateDetRatios(phi_list, vp_list, psiV_list, invRow_ptr_list, ratios);
  }
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::evaluateRatiosAlltoOne(ParticleSet& P, std::vector<Value>& ratios)
{
  {
    ScopedTimer local_timer(SPOVTimer);
    Phi->evaluateValue(P, -1, psiV_host_view);
  }
  auto& psiMinv = det_engine_.get_ref_psiMinv();
  for (int i = 0; i < psiMinv.rows(); i++)
    ratios[FirstIndex + i] = simd::dot(psiMinv[i], psiV.data(), NumOrbitals);
}


template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::resizeScratchObjectsForIonDerivs()
{
  grad_source_psiM.resize(NumPtcls, NumOrbitals);
  grad_lapl_source_psiM.resize(NumPtcls, NumOrbitals);
  grad_grad_source_psiM.resize(NumPtcls, NumOrbitals);
  phi_alpha_Minv.resize(NumPtcls, NumOrbitals);
  grad_phi_Minv.resize(NumPtcls, NumOrbitals);
  lapl_phi_Minv.resize(NumPtcls, NumOrbitals);
  grad_phi_alpha_Minv.resize(NumPtcls, NumOrbitals);
}

template<typename DET_ENGINE>
typename DiracDeterminantBatched<DET_ENGINE>::Grad DiracDeterminantBatched<DET_ENGINE>::evalGradSource(
    ParticleSet& P,
    ParticleSet& source,
    int iat)
{
  Grad g(0.0);
  if (Phi->hasIonDerivs())
  {
    resizeScratchObjectsForIonDerivs();
    Phi->evaluateGradSource(P, FirstIndex, LastIndex, source, iat, grad_source_psiM);
    auto& psiMinv = det_engine_.get_psiMinv();
    // psiMinv columns have padding but grad_source_psiM ones don't
    for (int i = 0; i < psiMinv.rows(); i++)
      g += simd::dot(psiMinv[i], grad_source_psiM[i], NumOrbitals);
  }

  return g;
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::evaluateHessian(ParticleSet& P, HessVector& grad_grad_psi)
{
  auto& psiMinv = det_engine_.get_ref_psiMinv();
  // Hessian is not often used, so only resize/allocate if used
  grad_grad_source_psiM.resize(psiMinv.rows(), NumOrbitals);
  //IM A HACK.  Assumes evaluateLog has already been executed.
  Matrix<Value> psiM_temp_host(psiM_temp.data(), psiM_temp.rows(), psiM_temp.cols());
  Phi->evaluate_notranspose(P, FirstIndex, LastIndex, psiM_temp_host, dpsiM, grad_grad_source_psiM);
  invertPsiM(psiM_temp, det_engine_.get_ref_psiMinv());

  phi_alpha_Minv      = 0.0;
  grad_phi_Minv       = 0.0;
  lapl_phi_Minv       = 0.0;
  grad_phi_alpha_Minv = 0.0;
  //grad_grad_psi.resize(NumPtcls);

  for (int i = 0, iat = FirstIndex; i < NumPtcls; i++, iat++)
  {
    Grad rv = simd::dot(psiMinv[i], dpsiM[i], NumOrbitals);
    //  HessType hess_tmp=simd::dot(psiMinv[i],grad_grad_source_psiM[i],NumOrbitals);
    Hess hess_tmp;
    hess_tmp           = 0.0;
    hess_tmp           = simd::dot(psiMinv[i], grad_grad_source_psiM[i], NumOrbitals);
    grad_grad_psi[iat] = hess_tmp - outerProduct(rv, rv);
  }
}

template<typename DET_ENGINE>
typename DiracDeterminantBatched<DET_ENGINE>::Grad DiracDeterminantBatched<DET_ENGINE>::evalGradSource(
    ParticleSet& P,
    ParticleSet& source,
    int iat,
    TinyVector<ParticleSet::ParticleGradient, OHMMS_DIM>& grad_grad,
    TinyVector<ParticleSet::ParticleLaplacian, OHMMS_DIM>& lapl_grad)
{
  Grad gradPsi(0.0);
  if (Phi->hasIonDerivs())
  {
    resizeScratchObjectsForIonDerivs();
    Phi->evaluateGradSource(P, FirstIndex, LastIndex, source, iat, grad_source_psiM, grad_grad_source_psiM,
                            grad_lapl_source_psiM);

    auto& psiMinv = det_engine_.get_psiMinv();

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
          lapl_phi_Minv(i, j) += d2psiM(i, k) * psiMinv(j, k);
      }
    for (int dim = 0; dim < OHMMS_DIM; dim++)
    {
      for (int i = 0; i < NumPtcls; i++)
        for (int j = 0; j < NumOrbitals; j++)
        {
          for (int k = 0; k < NumOrbitals; k++)
          {
            phi_alpha_Minv(i, j)[dim] += grad_source_psiM(i, k)[dim] * psiMinv(j, k);
            grad_phi_Minv(i, j)[dim] += dpsiM(i, k)[dim] * psiMinv(j, k);
            for (int dim_el = 0; dim_el < OHMMS_DIM; dim_el++)
              grad_phi_alpha_Minv(i, j)(dim, dim_el) += grad_grad_source_psiM(i, k)(dim, dim_el) * psiMinv(j, k);
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
        gradPsi += grad_source_psiM(i, j) * psiMinv(i, j);
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
          lapl_grad[dim][iel] += grad_lapl_source_psiM(i, j)[dim] * psiMinv(i, j);
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
template<typename DET_ENGINE>
typename DiracDeterminantBatched<DET_ENGINE>::LogValue DiracDeterminantBatched<DET_ENGINE>::evaluateLog(
    const ParticleSet& P,
    ParticleSet::ParticleGradient& G,
    ParticleSet::ParticleLaplacian& L)
{
  recompute(P);
  computeGL(G, L);
  return log_value_;
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::mw_evaluateLog(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                         const RefVectorWithLeader<ParticleSet>& p_list,
                                                         const RefVector<ParticleSet::ParticleGradient>& G_list,
                                                         const RefVector<ParticleSet::ParticleLaplacian>& L_list) const
{
  assert(this == &wfc_list.getLeader());
  const std::vector<bool> recompute_all(wfc_list.size(), true);
  mw_recompute(wfc_list, p_list, recompute_all);

  for (int iw = 0; iw < wfc_list.size(); iw++)
  {
    auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE>>(iw);
    det.computeGL(G_list[iw], L_list[iw]);
  }
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::recompute(const ParticleSet& P)
{
  {
    ScopedTimer spo_timer(SPOVGLTimer);

    UpdateMode = ORB_WALKER;
    Phi->evaluate_notranspose(P, FirstIndex, LastIndex, psiM_host, dpsiM, d2psiM);
    auto* psiM_vgl_ptr = psiM_vgl.data();
    // transfer host to device, total size 4, g(3) + l(1)
    PRAGMA_OFFLOAD("omp target update to(psiM_vgl_ptr[psiM_vgl.capacity():psiM_vgl.capacity()*4])")
  }
  invertPsiM(psiM_temp, det_engine_.get_ref_psiMinv());
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::mw_recompute(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                       const RefVectorWithLeader<ParticleSet>& p_list,
                                                       const std::vector<bool>& recompute) const
{
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<DET_ENGINE>>();
  const auto nw    = wfc_list.size();

  RefVectorWithLeader<WaveFunctionComponent> wfc_filtered_list(wfc_list.getLeader());
  RefVectorWithLeader<ParticleSet> p_filtered_list(p_list.getLeader());
  RefVectorWithLeader<SPOSet> phi_list(*wfc_leader.Phi);
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

      auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE>>(iw);
      phi_list.push_back(*det.Phi);
      psiM_temp_list.push_back(det.psiM_temp);
      psiM_host_list.push_back(det.psiM_host);
      dpsiM_list.push_back(det.dpsiM);
      d2psiM_list.push_back(det.d2psiM);
      // We need get_ref_psiMinv because C++ can't deduce the correct overload from
      // return type.
      psiMinv_list.push_back(det.get_det_engine().get_ref_psiMinv());
    }

  if (!wfc_filtered_list.size())
    return;

  {
    ScopedTimer spo_timer(wfc_leader.SPOVGLTimer);
    // I think through the magic of OMPtarget psiM_host_list actually results in psiM_temp being updated
    // on the device. For dspiM_list, d2psiM_list I think they are calculated on CPU and this is not true
    // This is the reason for the strange look omp target update to below.
    wfc_leader.Phi->mw_evaluate_notranspose(phi_list, p_filtered_list, wfc_leader.FirstIndex, wfc_leader.LastIndex,
                                            psiM_host_list, dpsiM_list, d2psiM_list);
  }

  mw_invertPsiM(wfc_filtered_list, psiM_temp_list, psiMinv_list);

  { // transfer dpsiM, d2psiM, psiMinv to device
    ScopedTimer d2h(H2DTimer);

    RefVector<DualVGLVector> psiM_vgl_list;
    psiM_vgl_list.reserve(nw);
    for (int iw = 0; iw < wfc_filtered_list.size(); iw++)
    {
      auto& det = wfc_filtered_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE>>(iw);
      psiM_vgl_list.push_back(det.psiM_vgl);
      det.UpdateMode = ORB_WALKER;
    }

    // transfer host to device, total size 4, g(3) + l(1), skipping v
    DET_ENGINE::mw_transferVGL_H2D(wfc_leader.det_engine_, psiM_vgl_list, 1, 4);
  }
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::evaluateDerivatives(ParticleSet& P,
                                                              const opt_variables_type& active,
                                                              std::vector<Value>& dlogpsi,
                                                              std::vector<Value>& dhpsioverpsi)
{
  Phi->evaluateDerivatives(P, active, dlogpsi, dhpsioverpsi, FirstIndex, LastIndex);
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::registerTWFFastDerivWrapper(const ParticleSet& P,
                                                                      TWFFastDerivWrapper& twf) const
{
  twf.addGroup(P, P.getGroupID(FirstIndex), Phi.get());
}

template<typename DET_ENGINE>
std::unique_ptr<DiracDeterminantBase> DiracDeterminantBatched<DET_ENGINE>::makeCopy(std::unique_ptr<SPOSet>&& spo) const
{
  return std::make_unique<DiracDeterminantBatched<DET_ENGINE>>(std::move(spo), FirstIndex, LastIndex, ndelay_,
                                                               matrix_inverter_kind_);
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::createResource(ResourceCollection& collection) const
{
  collection.addResource(std::make_unique<DiracDeterminantBatchedMultiWalkerResource>());
  Phi->createResource(collection);
  det_engine_.createResource(collection);
  collection.addResource(std::make_unique<typename DET_ENGINE::DetInverter>());
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::acquireResource(
    ResourceCollection& collection,
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const
{
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<DET_ENGINE>>();
  auto res_ptr     = dynamic_cast<DiracDeterminantBatchedMultiWalkerResource*>(collection.lendResource().release());
  if (!res_ptr)
    throw std::runtime_error("DiracDeterminantBatched::acquireResource dynamic_cast failed");
  wfc_leader.mw_res_.reset(res_ptr);

  RefVectorWithLeader<SPOSet> phi_list(*wfc_leader.Phi);
  for (WaveFunctionComponent& wfc : wfc_list)
  {
    auto& det = static_cast<DiracDeterminantBatched<DET_ENGINE>&>(wfc);
    phi_list.push_back(*det.Phi);
  }
  wfc_leader.Phi->acquireResource(collection, phi_list);

  wfc_leader.det_engine_.acquireResource(collection);

  auto det_eng_ptr = dynamic_cast<typename DET_ENGINE::DetInverter*>(collection.lendResource().release());
  if (!det_eng_ptr)
    throw std::runtime_error(
        "DiracDeterminantBatched::acquireResource dynamic_cast to DET_ENGINE::DetInverter* failed");
  wfc_leader.accel_inverter_.reset(det_eng_ptr);
}

template<typename DET_ENGINE>
void DiracDeterminantBatched<DET_ENGINE>::releaseResource(
    ResourceCollection& collection,
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const
{
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<DET_ENGINE>>();
  collection.takebackResource(std::move(wfc_leader.mw_res_));
  RefVectorWithLeader<SPOSet> phi_list(*wfc_leader.Phi);
  for (WaveFunctionComponent& wfc : wfc_list)
  {
    auto& det = static_cast<DiracDeterminantBatched<DET_ENGINE>&>(wfc);
    phi_list.push_back(*det.Phi);
  }
  wfc_leader.Phi->releaseResource(collection, phi_list);
  wfc_leader.det_engine_.releaseResource(collection);
  collection.takebackResource(std::move(wfc_leader.accel_inverter_));
}

template class DiracDeterminantBatched<>;
#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
template class DiracDeterminantBatched<MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
#endif

} // namespace qmcplusplus
