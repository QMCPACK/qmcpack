
#include "DiracDeterminantBatched.h"

namespace qmcplusplus
{
template<>
inline void DiracDeterminantBatched::mw_recompute<MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                   const RefVectorWithLeader<ParticleSet>& p_list,
                                                   const std::vector<bool>& recompute)
{
  using DetEngine = MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>;
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<DetEngine>>();
  const auto nw    = wfc_list.size();

  {
    ScopedTimer spo_timer(wfc_leader.SPOVGLTimer);

    RefVectorWithLeader<SPOSet> phi_list(*wfc_leader.Phi);
    RefVector<GradMatrix_t> dpsiM_list;
    RefVector<ValueMatrix_t> d2psiM_list;
    phi_list.reserve(wfc_list.size());
    dpsiM_list.reserve(nw);
    d2psiM_list.reserve(nw);
    std::vector<DDBT::ValueMatrix_t> psiM_temp_host_list;
    psiM_temp_host_list.reserve(nw);

    for (int iw = 0; iw < nw; iw++)
    {
      auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<DetEngine>>(iw);
      phi_list.push_back(*det.Phi);
      psiM_temp_host_list.emplace_back(det.psiM_temp.data(), det.psiM_temp.rows(), det.psiM_temp.cols());
      dpsiM_list.push_back(det.dpsiM);
      d2psiM_list.push_back(det.d2psiM);
    }

    wfc_leader.Phi->mw_evaluate_notranspose(phi_list, p_list, wfc_leader.FirstIndex, wfc_leader.LastIndex,
                                            makeRefVector<DDBT::ValueMatrix_t>(psiM_temp_host_list), dpsiM_list, d2psiM_list);
  }
  RefVector<DDBT::OffloadPinnedValueMatrix_t> psiM_temp_list;
  psiM_temp_list.reserve(nw);

  for (int iw = 0; iw < nw; iw++)
  {
    auto& det          = wfc_list.getCastedElement<DiracDeterminantBatched<DetEngine>>(iw);
    auto* psiM_vgl_ptr = det.psiM_vgl.data();
    size_t stride      = wfc_leader.psiM_vgl.capacity();
    PRAGMA_OFFLOAD("omp target update to(psiM_vgl_ptr[stride:stride*4]) nowait")
    psiM_temp_list.push_back(det.psiM_temp);
  }
  DDBT::mw_invertPsiM(*(wfc_leader.mw_res_), wfc_list, psiM_temp_list);
  PRAGMA_OFFLOAD("omp taskwait")
}
  
template<>
void DiracDeterminantBatched<MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>::mw_invertPsiM(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                             const RefVector<const ValueMatrix_t>& logdetT_list)
{
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<DET_ENGINE_TYPE>>();
  ScopedTimer inverse_timer(&wfc_leader.InverseTimer);
  const auto nw = wfc_list.size();

  RefVector<LogValueType> log_value_list;

  for (int iw = 0; iw < nw; iw++)
  {
    auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE_TYPE>>(iw);
    log_value_list.push_back(det.LogValue);
  }
  wfc_leader.det_inverter_.mw_invert_transpose(*cuda_handles_,logdetT_list, log_value_list);
}

template<>
inline void DiracDeterminantBatched<MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>::mw_evaluateLog(
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
    const RefVectorWithLeader<ParticleSet>& p_list,
    const RefVector<ParticleSet::ParticleGradient_t>& G_list,
    const RefVector<ParticleSet::ParticleLaplacian_t>& L_list) const
{
  using DetEngine = MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>;
  assert(this == &wfc_list.getLeader());
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<DetEngine>>();
  wfc_leader.guardMultiWalkerRes();
  auto& mw_res       = *wfc_leader.mw_res_;
  mw_res.log_values.resize(wfc_list.size());
  const std::vector<bool> recompute_all(wfc_list.size(), true);
  mw_recompute(wfc_list, p_list, recompute_all);

  for (int iw = 0; iw < wfc_list.size(); iw++)
  {
    auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<DetEngine>>(iw);
    det.computeGL(G_list[iw], L_list[iw]);
  }
}

  
}
