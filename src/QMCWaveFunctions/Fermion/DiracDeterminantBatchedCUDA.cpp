
/**@file 
 * @brief CUDA specific specialization of DiracDeterminantBatched methods.
 */

#include "DiracDeterminantBatched.h"

namespace qmcplusplus
{
// template<>
// void DiracDeterminantBatched<MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>::mw_recompute_impl(
//     const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
//     const RefVectorWithLeader<ParticleSet>& p_list,
//     const std::vector<bool>& recompute)
// {
//   using DetEngine  = MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>;
//   auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<DetEngine>>();
//   const auto nw    = wfc_list.size();

// std:
//   cout << "You have the CUDA mw_recompute!\n";
//   {
//     ScopedTimer spo_timer(wfc_leader.SPOVGLTimer);

//     RefVectorWithLeader<SPOSet> phi_list(*wfc_leader.Phi);
//     RefVector<GradMatrix_t> dpsiM_list;
//     RefVector<ValueMatrix_t> d2psiM_list;
//     phi_list.reserve(wfc_list.size());
//     dpsiM_list.reserve(nw);
//     d2psiM_list.reserve(nw);
//     std::vector<DDBT::ValueMatrix_t> psiM_temp_host_list;
//     psiM_temp_host_list.reserve(nw);

//     for (int iw = 0; iw < nw; iw++)
//     {
//       auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<DetEngine>>(iw);
//       phi_list.push_back(*det.Phi);
//       psiM_temp_host_list.emplace_back(det.psiM_temp.data(), det.psiM_temp.rows(), det.psiM_temp.cols());
//       dpsiM_list.push_back(det.dpsiM);
//       d2psiM_list.push_back(det.d2psiM);
//     }

//     wfc_leader.Phi->mw_evaluate_notranspose(phi_list, p_list, wfc_leader.FirstIndex, wfc_leader.LastIndex,
//                                             makeRefVector<DDBT::ValueMatrix_t>(psiM_temp_host_list), dpsiM_list,
//                                             d2psiM_list);
//   }
//   RefVector<DDBT::OffloadPinnedValueMatrix_t> psiM_temp_list;
//   psiM_temp_list.reserve(nw);

//   for (int iw = 0; iw < nw; iw++)
//   {
//     auto& det          = wfc_list.getCastedElement<DiracDeterminantBatched<DetEngine>>(iw);
//     auto* psiM_vgl_ptr = det.psiM_vgl.data();
//     size_t stride      = wfc_leader.psiM_vgl.capacity();
//     PRAGMA_OFFLOAD("omp target update to(psiM_vgl_ptr[stride:stride*4]) nowait")
//     psiM_temp_list.push_back(det.psiM_temp);
//   }
//   DDBT::mw_invertPsiM(*(wfc_leader.mw_res_), wfc_list, psiM_temp_list);
//   PRAGMA_OFFLOAD("omp taskwait")
// }

class SPOSet;

template<>
void DiracDeterminantDetails<DiracDeterminantBatched<MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>,
    MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>::mw_recomputeDispatch(
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
    const RefVectorWithLeader<ParticleSet>& p_list,
    const std::vector<bool>& recompute)
{
  using DetEngine  = MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>;
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<DetEngine>>();
  const auto nw    = wfc_list.size();
  using DDBT       = std::decay_t<decltype(wfc_leader)>;
  std::cout << "wpppp CUDA\n";
  {
    ScopedTimer spo_timer(wfc_leader.SPOVGLTimer);

    RefVectorWithLeader<SPOSet> phi_list(*wfc_leader.Phi);
    RefVector<DDBT::GradMatrix_t> dpsiM_list;
    RefVector<DDBT::ValueMatrix_t> d2psiM_list;
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
                                            makeRefVector<DDBT::ValueMatrix_t>(psiM_temp_host_list), dpsiM_list,
                                            d2psiM_list);
  }
  RefVector<const DDBT::OffloadPinnedValueMatrix_t> psiM_temp_list;
  psiM_temp_list.reserve(nw);

  for (int iw = 0; iw < nw; iw++)
  {
    auto& det          = wfc_list.getCastedElement<DiracDeterminantBatched<DetEngine>>(iw);
    auto* psiM_vgl_ptr = det.psiM_vgl.data();
    size_t stride      = wfc_leader.psiM_vgl.capacity();
    PRAGMA_OFFLOAD("omp target update to(psiM_vgl_ptr[stride:stride*4]) nowait")
    psiM_temp_list.push_back(det.psiM_temp);
  }
  DiracDeterminantBatched<DetEngine>::mw_invertPsiM(wfc_list, psiM_temp_list);
  PRAGMA_OFFLOAD("omp taskwait")
}

} // namespace qmcplusplus
