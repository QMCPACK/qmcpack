//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
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
}

template<typename DET_ENGINE = MatrixUpdateOMPTarget<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>
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
  using OffloadPinnedLogValueVector = Vector<LogValue, typename DET_ENGINE::template OffloadPinnedAllocator<LogValue>>;

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

template<typename DET_ENGINE = MatrixUpdateOMPTarget<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>
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
  using LogValue      = std::complex<QTFull::RealType>;
  using DetEngine_t   = DET_ENGINE;
  
  template<typename DT>
  using OffloadPinnedAllocator        = OMPallocator<DT, PinnedAlignedAllocator<DT>>;
  using OffloadPinnedValueVector_t    = Vector<ValueType, OffloadPinnedAllocator<ValueType>>;
  using OffloadPinnedLogValueVector_t = Vector<LogValue, OffloadPinnedAllocator<LogValue>>;
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
  LogValueType evaluateLog(ParticleSet& P,
                           ParticleSet::ParticleGradient_t& G,
                           ParticleSet::ParticleLaplacian_t& L) override;

  void mw_evaluateLog(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                      const RefVectorWithLeader<ParticleSet>& p_list,
                      const RefVector<ParticleSet::ParticleGradient_t>& G_list,
                      const RefVector<ParticleSet::ParticleLaplacian_t>& L_list) const override;

  void recompute(DiracDeterminantBatchedMultiWalkerResource<DET_ENGINE>& mw_res, ParticleSet& P);

  LogValueType evaluateGL(ParticleSet& P,
                          ParticleSet::ParticleGradient_t& G,
                          ParticleSet::ParticleLaplacian_t& L,
                          bool fromscratch) override;

  void mw_evaluateGL(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                     const RefVectorWithLeader<ParticleSet>& p_list,
                     const RefVector<ParticleSet::ParticleGradient_t>& G_list,
                     const RefVector<ParticleSet::ParticleLaplacian_t>& L_list,
                     bool fromscratch) const override;

  void evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi) override;

  void createResource(ResourceCollection& collection) override;
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

  /// return  for testing
  auto& getPsiMinv() const { return psiMinv; }

  auto& get_det_engine() { return det_engine_; }

  /// inverse transpose of psiM(j,i) \f$= \psi_j({\bf r}_i)\f$, actual memory owned by det_inverter_
  ValueMatrix_t psiMinv;

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

  LogValue get_log_value() const { return log_value_; }

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

private:
  /// Smelly second source of truth for this DDB's logvalue in the mw_res_;
  LogValue log_value_;

  /// compute G adn L assuming psiMinv, dpsiM, d2psiM are ready for use
  void computeGL(ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L) const;

  /// invert logdetT(psiM), result is in the engine.
  void invertPsiM(DiracDeterminantBatchedMultiWalkerResource<DET_ENGINE>& mw_res, OffloadPinnedValueMatrix_t& logdetT);

  static void mw_invertPsiM(DiracDeterminantBatchedMultiWalkerResource<DET_ENGINE>& mw_res,
                            const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                            RefVector<OffloadPinnedValueMatrix_t>& logdetT_list);

  static void mw_recompute(DiracDeterminantBatchedMultiWalkerResource<DET_ENGINE>& mw_res,
                           const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                           const RefVectorWithLeader<ParticleSet>& p_list);

  /// Resize all temporary arrays required for force computation.
  void resizeScratchObjectsForIonDerivs();

  /// maximal number of delayed updates
  int ndelay;

  /// timers
  NewTimer &D2HTimer, &H2DTimer;

  friend class qmcplusplus::testing::DiracDeterminantBatchedTest;
};


extern template struct DiracDeterminantBatchedMultiWalkerResource<>;
extern template class DiracDeterminantBatched<>;
#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
extern template struct DiracDeterminantBatchedMultiWalkerResource<
    MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
extern template class DiracDeterminantBatched<
    MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
#endif

} // namespace qmcplusplus
#endif
