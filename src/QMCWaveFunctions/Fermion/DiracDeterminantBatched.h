/////////////////////////////////////////////////////////////////////////////////////
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


/**@file DiracDeterminantBatched.h
 * @brief Declaration of DiracDeterminantBatched with a S(ingle)P(article)O(rbital)Set
 */
#ifndef QMCPLUSPLUS_DIRACDETERMINANTBATCHED_H
#define QMCPLUSPLUS_DIRACDETERMINANTBATCHED_H

#include "Configuration.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#include "QMCWaveFunctions/Fermion/MatrixUpdateOMPTarget.h"
#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
#include "QMCWaveFunctions/Fermion/MatrixDelayedUpdateCUDA.h"
#endif
#include "Platforms/PinnedAllocator.h"
#include "OMPTarget/OMPallocator.hpp"

namespace qmcplusplus
{
namespace testing
{
class DiracDeterminantBatchedTest;
struct SetupDiracDetResources;
}

template<typename DET_ENGINE>
struct DiracDeterminantBatchedMultiWalkerResource : public Resource
{
  using ValueType = QMCTraits::ValueType;
  using GradType  = QMCTraits::GradType;
  using Real      = QMCTraits::RealType;
  template<typename DT>
  using OffloadPinnedAllocator = OMPallocator<DT, PinnedAlignedAllocator<DT>>;
  using OffloadVGLVector_t     = VectorSoaContainer<ValueType, QMCTraits::DIM + 2, OffloadPinnedAllocator<ValueType>>;
  // I don't think its a good idea create a hard dependency all the way back to WaveFunctionComponent for this.
  using LogValue                    = std::complex<Real>;
  using OffloadPinnedLogValueVector = Vector<LogValue, OffloadPinnedAllocator<LogValue>>;

  DiracDeterminantBatchedMultiWalkerResource() : Resource("DiracDeterminantBatched") {}

  DiracDeterminantBatchedMultiWalkerResource(const DiracDeterminantBatchedMultiWalkerResource&)
      : DiracDeterminantBatchedMultiWalkerResource()
  {}

  Resource* makeClone() const override { return new DiracDeterminantBatchedMultiWalkerResource(*this); }

  OffloadPinnedLogValueVector log_values;
  /// value, grads, laplacian of single-particle orbital for particle-by-particle update and multi walker [5][nw*norb]
  OffloadVGLVector_t phi_vgl_v;
  // multi walker of ratio
  std::vector<ValueType> ratios_local;
  // multi walker of grads
  std::vector<GradType> grad_new_local;
};


// class DiracDeterminantDetails
// {
// public:
// template<class DET_ENGINE>
// static void mw_recomputeDispatch(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
//                           const RefVectorWithLeader<ParticleSet>& p_list,
//                           const std::vector<bool>& recompute);
// };

template<typename DET_ENGINE>
class DiracDeterminantBatched : public DiracDeterminantBase
{
public:
  using ValueVector_t = SPOSet::ValueVector_t;
  using ValueMatrix_t = SPOSet::ValueMatrix_t;
  using GradVector_t  = SPOSet::GradVector_t;
  using GradMatrix_t  = SPOSet::GradMatrix_t;
  using HessMatrix_t  = SPOSet::HessMatrix_t;
  using HessVector_t  = SPOSet::HessVector_t;
  using HessType      = SPOSet::HessType;
  using Real          = QMCTraits::RealType;
  using mValueType    = QMCTraits::QTFull::ValueType;
  using mGradType     = TinyVector<mValueType, DIM>;
  using DetEngine_t   = DET_ENGINE;

  template<typename DT>
  using OffloadPinnedAllocator        = OMPallocator<DT, PinnedAlignedAllocator<DT>>;
  using OffloadPinnedValueVector_t    = Vector<ValueType, OffloadPinnedAllocator<ValueType>>;
  using OffloadPinnedLogValueVector_t = Vector<LogValueType, OffloadPinnedAllocator<LogValueType>>;
  using OffloadPinnedValueMatrix_t    = Matrix<ValueType, OffloadPinnedAllocator<ValueType>>;
  using OffloadPinnedPsiValueVector_t = Vector<PsiValueType, OffloadPinnedAllocator<PsiValueType>>;
  using OffloadVGLVector_t            = VectorSoaContainer<ValueType, DIM + 2, OffloadPinnedAllocator<ValueType>>;

  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   */
  DiracDeterminantBatched(std::shared_ptr<SPOSet>&& spos, int first = 0);

  // copy constructor and assign operator disabled
  DiracDeterminantBatched(const DiracDeterminantBatched& s) = delete;
  DiracDeterminantBatched& operator=(const DiracDeterminantBatched& s) = delete;

  /** set the index of the first particle in the determinant and reset the size of the determinant
   *@param first index of first particle
   *@param nel number of particles in the determinant
   */
  void set(int first, int nel, int delay = 1) override final;

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           std::vector<ValueType>& dlogpsi,
                           std::vector<ValueType>& dhpsioverpsi) override;

  ///reset the size: with the number of particles and number of orbtials
  void resize(int nel, int morb);

  void registerData(ParticleSet& P, WFBufferType& buf) override;

  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;

  /** return the ratio only for the  iat-th partcle move
   * @param P current configuration
   * @param iat the particle thas is being moved
   */
  PsiValueType ratio(ParticleSet& P, int iat) override;

  void mw_calcRatio(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    int iat,
                    std::vector<PsiValueType>& ratios) const override;

  /** compute multiple ratios for a particle move
   */
  void evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios) override;

  void mw_evaluateRatios(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                         const RefVector<const VirtualParticleSet>& vp_list,
                         std::vector<std::vector<ValueType>>& ratios) const override;

  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;

  void mw_ratioGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    int iat,
                    std::vector<PsiValueType>& ratios,
                    std::vector<GradType>& grad_new) const override;

  GradType evalGrad(ParticleSet& P, int iat) override;

  void mw_evalGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                   const RefVectorWithLeader<ParticleSet>& p_list,
                   int iat,
                   std::vector<GradType>& grad_now) const override;

  GradType evalGradSource(ParticleSet& P, ParticleSet& source, int iat) override;

  GradType evalGradSource(ParticleSet& P,
                          ParticleSet& source,
                          int iat,
                          TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM>& grad_grad,
                          TinyVector<ParticleSet::ParticleLaplacian_t, OHMMS_DIM>& lapl_grad) override;

  /** move was accepted, update the real container
   */
  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override;

  void mw_accept_rejectMove(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                            const RefVectorWithLeader<ParticleSet>& p_list,
                            int iat,
                            const std::vector<bool>& isAccepted,
                            bool safe_to_delay = false) const override;

  /** complete any left over determinant matrix updates.
   * Usually this is the end of pbyp moves for a given spin of electrons
   * The psiM, dpsiM, d2psiM should be up-to-date on both device and host sides.
   */
  void completeUpdates() override;

  void mw_completeUpdates(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override;

  /** move was rejected. copy the real container to the temporary to move on
   */
  void restore(int iat) override;

  /** evaluate log of a determinant for a particle set
   * This is the most defensive call. The psiM, dpsiM, d2psiM should be up-to-date on both device and host sides.
   */
  LogValueType evaluateLog(const ParticleSet& P,
                           ParticleSet::ParticleGradient_t& G,
                           ParticleSet::ParticleLaplacian_t& L) override;

  void mw_evaluateLog(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                      const RefVectorWithLeader<ParticleSet>& p_list,
                      const RefVector<ParticleSet::ParticleGradient_t>& G_list,
                      const RefVector<ParticleSet::ParticleLaplacian_t>& L_list) const override;

  void recompute(DiracDeterminantBatchedMultiWalkerResource<DET_ENGINE>& mw_res, const ParticleSet& P);

  void mw_recompute(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                          const RefVectorWithLeader<ParticleSet>& p_list,
                    const std::vector<bool>& recompute) const override;
  
  LogValueType evaluateGL(const ParticleSet& P,
                          ParticleSet::ParticleGradient_t& G,
                          ParticleSet::ParticleLaplacian_t& L,
                          bool fromscratch) override;

  void mw_evaluateGL(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                     const RefVectorWithLeader<ParticleSet>& p_list,
                     const RefVector<ParticleSet::ParticleGradient_t>& G_list,
                     const RefVector<ParticleSet::ParticleLaplacian_t>& L_list,
                     bool fromscratch) const override;

  void evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi) override;

  void createResource(ResourceCollection& collection) const override;
  void acquireResource(ResourceCollection& collection) override;
  void releaseResource(ResourceCollection& collection) override;

  /** cloning function
   * @param tqp target particleset
   * @param spo spo set
   *
   * This interface is exposed only to SlaterDet and its derived classes
   * can overwrite to clone itself correctly.
   */
  DiracDeterminantBatched* makeCopy(std::shared_ptr<SPOSet>&& spo) const override;

  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios) override;

  auto& get_det_engine() { return det_engine_; }

  /// memory for psiM, dpsiM and d2psiM. [5][norb*norb]
  OffloadVGLVector_t psiM_vgl;
  /// psiM(j,i) \f$= \psi_j({\bf r}_i)\f$. partial memory view of psiM_vgl
  OffloadPinnedValueMatrix_t psiM_temp;
  /// dpsiM(i,j) \f$= \nabla_i \psi_j({\bf r}_i)\f$. partial memory view of psiM_vgl
  GradMatrix_t dpsiM;
  /// d2psiM(i,j) \f$= \nabla_i^2 \psi_j({\bf r}_i)\f$. partial memory view of psiM_vgl
  ValueMatrix_t d2psiM;

  /// Used for force computations
  GradMatrix_t grad_source_psiM, grad_lapl_source_psiM;
  HessMatrix_t grad_grad_source_psiM;

  GradMatrix_t phi_alpha_Minv, grad_phi_Minv;

  ValueMatrix_t lapl_phi_Minv;
  HessMatrix_t grad_phi_alpha_Minv;

  /// value of single-particle orbital for particle-by-particle update
  OffloadPinnedValueVector_t psiV;
  ValueVector_t psiV_host_view;
  GradVector_t dpsiV;
  ValueVector_t d2psiV;


  /// Log values for invert_transpose results
  //OffloadPinnedLogValueVector_t log_values;

  /// delayed update engine
  DET_ENGINE det_engine_;

  // psi(r')/psi(r) during a PbyP move
  PsiValueType curRatio;

  std::unique_ptr<DiracDeterminantBatchedMultiWalkerResource<DET_ENGINE>> mw_res_;

  LogValueType get_log_value() const { return LogValue; }

private:

  /// compute G adn L assuming psiMinv, dpsiM, d2psiM are ready for use
  void computeGL(ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L) const;

  /// invert logdetT(psiM), result is in the engine.
  void invertPsiM(DiracDeterminantBatchedMultiWalkerResource<DET_ENGINE>& mw_res, OffloadPinnedValueMatrix_t& logdetT);

  // static void mw_invertPsiM(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
  //                           const RefVector<const ValueMatrix_t>& logdetT_list);

  static void mw_invertPsiM(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                            const RefVector<const OffloadPinnedValueMatrix_t>& logdetT_list);

  /// Resize all temporary arrays required for force computation.
  void resizeScratchObjectsForIonDerivs();

  /// maximal number of delayed updates
  int ndelay;

  /// timers
  NewTimer &D2HTimer, &H2DTimer;

  // make this class unit tests friendly without the need of setup resources.
  void guardMultiWalkerRes()
  {
    if (!mw_res_)
    {
      std::cerr
          << "WARNING DiracDeterminantBatched : This message should not be seen in production (performance bug) runs "
             "but only unit tests (expected)."
          << std::endl;
      mw_res_ = std::make_unique<DiracDeterminantBatchedMultiWalkerResource<DET_ENGINE>>();
      mw_res_->log_values.resize(1);
    }
  }

  //  friend class qmcplusplus::DiracDeterminantDetails;
  friend struct qmcplusplus::testing::SetupDiracDetResources;
  friend class qmcplusplus::testing::DiracDeterminantBatchedTest;
};


// template<class DET_ENGINE>
// void DiracDeterminantDetails::mw_recomputeDispatch(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
//                           const RefVectorWithLeader<ParticleSet>& p_list,
//                           const std::vector<bool>& recompute)
// {
//   auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<DET_ENGINE>>();
//   const auto nw    = wfc_list.size();
//   using DDBT = decltype(wfc_leader);
//   RefVectorWithLeader<WaveFunctionComponent> wfc_filtered_list(wfc_list.getLeader());
//   RefVectorWithLeader<ParticleSet> p_filtered_list(p_list.getLeader());
//   RefVectorWithLeader<SPOSet> phi_list(*wfc_leader.Phi);
//   RefVector<ValueMatrix_t> psiM_temp_list;
//   RefVector<GradMatrix_t> dpsiM_list;
//   RefVector<ValueMatrix_t> d2psiM_list;

//   wfc_filtered_list.reserve(nw);
//   p_filtered_list.reserve(nw);
//   phi_list.reserve(nw);
//   psiM_temp_list.reserve(nw);
//   dpsiM_list.reserve(nw);
//   d2psiM_list.reserve(nw);

//   for (int iw = 0; iw < nw; iw++)
//     if (recompute[iw])
//     {
//       wfc_filtered_list.push_back(wfc_list[iw]);
//       p_filtered_list.push_back(p_list[iw]);

//       auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE>>(iw);
//       phi_list.push_back(*det.Phi);
//       psiM_temp_list.push_back(det.psiM_temp);
//       dpsiM_list.push_back(det.dpsiM);
//       d2psiM_list.push_back(det.d2psiM);
//     }

//   if (!wfc_filtered_list.size())
//     return;

//   {
//     ScopedTimer spo_timer(wfc_leader.SPOVGLTimer);
//     wfc_leader.Phi->mw_evaluate_notranspose(phi_list, p_filtered_list, wfc_leader.FirstIndex, wfc_leader.LastIndex,
//                                             psiM_temp_list, dpsiM_list, d2psiM_list);
//   }

//   { // transfer dpsiM, d2psiM, psiMinv to device
//     ScopedTimer d2h(DDBT::H2DTimer);

//     RefVector<const typename DDBT::ValueMatrix_t> const_psiM_temp_list;
//     for (int iw = 0; iw < wfc_filtered_list.size(); iw++)
//     {
//       auto& det          = wfc_filtered_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE>>(iw);
//       auto* psiM_vgl_ptr = det.psiM_vgl.data();
//       size_t stride      = wfc_leader.psiM_vgl.capacity();
//       PRAGMA_OFFLOAD("omp target update to(psiM_vgl_ptr[stride:stride*4]) nowait")
//       const_psiM_temp_list.push_back(det.psiM_temp);
//     }
//     DDBT::mw_invertPsiM(wfc_filtered_list, const_psiM_temp_list);
//     PRAGMA_OFFLOAD("omp taskwait")
//   }
// }

// template<>
// inline void DiracDeterminantDetails::mw_recomputeDispatch<MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
//                                                    const RefVectorWithLeader<ParticleSet>& p_list,
//                                                    const std::vector<bool>& recompute)
// {
//   using DetEngine = MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>;
//   auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<DetEngine>>();
//   const auto nw    = wfc_list.size();
//   using DDBT = std::decay_t<decltype(wfc_leader)>;

//   {
//     ScopedTimer spo_timer(wfc_leader.SPOVGLTimer);

//     RefVectorWithLeader<SPOSet> phi_list(*wfc_leader.Phi);
//     RefVector<DDBT::GradMatrix_t> dpsiM_list;
//     RefVector<DDBT::ValueMatrix_t> d2psiM_list;
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
//                                             makeRefVector<DDBT::ValueMatrix_t>(psiM_temp_host_list), dpsiM_list, d2psiM_list);
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


extern template struct DiracDeterminantBatchedMultiWalkerResource<MatrixUpdateOMPTarget<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
extern template class DiracDeterminantBatched<MatrixUpdateOMPTarget<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
extern template struct DiracDeterminantBatchedMultiWalkerResource<
    MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
extern template class DiracDeterminantBatched<
    MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
#endif

} // namespace qmcplusplus
#endif
