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


/**@file DiracDeterminantBatched.h
 * @brief Declaration of DiracDeterminantBatched with a S(ingle)P(article)O(rbital)Set
 */
#ifndef QMCPLUSPLUS_DIRACDETERMINANTBATCHED_H
#define QMCPLUSPLUS_DIRACDETERMINANTBATCHED_H

#include "Configuration.h"
#include "DeterminantAllocators.hpp"
#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#if defined(ENABLE_CUDA)
#include "QMCWaveFunctions/Fermion/MatrixDelayedUpdateCUDA.h"
#endif
#if defined(ENABLE_OFFLOAD)
#include "QMCWaveFunctions/Fermion/MatrixUpdateOMPTarget.h"
#endif
#include "Resource.h"

namespace qmcplusplus
{
namespace testing
{
class DiracDeterminantBatchedTest;
struct SetupDiracDetResources;
} // namespace testing

/** Should contain the Multiwalker data.
 *  in truth this is not used in many cases and instead single walker data stashed in the "extra" copies of
 *  DiracDeterminantBatched are used..
 */
struct DiracDeterminantBatchedMultiWalkerResource : public Resource
{
  using ValueType = QMCTraits::ValueType;
  using GradType  = QMCTraits::GradType;
  using Real      = QMCTraits::RealType;
  using FullPrecReal = QMCTraits::QTFull::RealType;
  using DualVGLVector_t     = VectorSoaContainer<ValueType, QMCTraits::DIM + 2, PinnedDualAllocator<ValueType>>;
  // I don't think its a good idea create a hard dependency all the way back to WaveFunctionComponent for this.
  using LogValue                    = std::complex<FullPrecReal>;
  using DualPinnedLogValueVector = Vector<LogValue, PinnedDualAllocator<LogValue>>;

  DiracDeterminantBatchedMultiWalkerResource() : Resource("DiracDeterminantBatched") {}

  DiracDeterminantBatchedMultiWalkerResource(const DiracDeterminantBatchedMultiWalkerResource&)
      : DiracDeterminantBatchedMultiWalkerResource()
  {}

  Resource* makeClone() const override { return new DiracDeterminantBatchedMultiWalkerResource(*this); }

  DualPinnedLogValueVector log_values;
  /// value, grads, laplacian of single-particle orbital for particle-by-particle update and multi walker [5][nw*norb]
  DualVGLVector_t phi_vgl_v;
  // multi walker of ratio
  std::vector<ValueType> ratios_local;
  // multi walker of grads
  std::vector<GradType> grad_new_local;
};

/** Currently DiracDeterminantBatched (DDB) is an object managing batched computation and update of DiracDeterminants,
 *  the multiwalker resources needed to accomplish this, and is also the container element class for all the 
 *  determinant single walker data members and as such is derived from DiracDeterminantBase : WaveFunctionComponent. 
 *  One of these perwalker members is the DET_ENGINE. The DET_ENGINE is specialized for
 *  portability framework and/or accelerator. The per walker psiMinv is owned by the DET_ENGINE. The psiM its derivative and
 *  temp elements are still ownded by the DDB.
 *
 */
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

  using DualPinnedValueVector_t    = Vector<ValueType, PinnedDualAllocator<ValueType>>;
  using DualPinnedLogValueVector_t = Vector<LogValueType, PinnedDualAllocator<LogValueType>>;
  using DualPinnedValueMatrix_t    = Matrix<ValueType, PinnedDualAllocator<ValueType>>;
  using DualPinnedPsiValueVector_t = Vector<PsiValueType, PinnedDualAllocator<PsiValueType>>;
  using DualPinnedGradVector       = Vector<OrbitalSetTraits<SPOSet::ValueType>::GradType, PinnedDualAllocator<OrbitalSetTraits<SPOSet::ValueType>::GradType>>;
  using DualPinnedGradMatrix       = Matrix<OrbitalSetTraits<SPOSet::ValueType>::GradType, PinnedDualAllocator<OrbitalSetTraits<SPOSet::ValueType>::GradType>>;

  using DualVGLVector_t            = VectorSoaContainer<ValueType, DIM + 2, PinnedDualAllocator<ValueType>>;

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
  void set(int first, int nel, int delay = 1) final;

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           std::vector<ValueType>& dlogpsi,
                           std::vector<ValueType>& dhpsioverpsi) override;

  ///reset the size: with the number of particles and number of orbtials
  void resize(int nel, int morb);

  void mw_resize(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list, int nel, int morb);

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
                         const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                         std::vector<std::vector<ValueType>>& ratios) const override;

  /** Legacy single method */
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

  void recompute(DiracDeterminantBatchedMultiWalkerResource& mw_res, const ParticleSet& P);

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
  // This sucks for lower level testing
  void acquireResource(ResourceCollection& collection, const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override;
  void releaseResource(ResourceCollection& collection, const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override;

  /** cloning function
   * @param tqp target particleset
   * @param spo spo set
   *
   * This interface is exposed only to SlaterDet and its derived classes
   * can overwrite to clone itself correctly.
   */
  DiracDeterminantBatched* makeCopy(std::shared_ptr<SPOSet>&& spo) const override;

  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios) override;

  DET_ENGINE& get_det_engine() { return det_engine_; }

  /** @defgroup LegacySingleData Single Walker Data Members of Legacy OO design
   *  @brief    Deprecated as high throughput of walkers requires a division between
   *            walker data which should be "SoA" and traditional OO design which is generally AoS with
   *            single "structure" functions bundled.
   *  
   *  @ingroup LegacySingleData
   *  @{
   */
  /// Legacy single determinant values of single-particle orbital for particle-by-particle update
  /// Ideally DDB should use the mw_res resources and not these duplicative values.

  /// memory for psiM, dpsiM and d2psiM. [5][norb*norb]
  DualVGLVector_t psiM_vgl;
  /// psiM(j,i) \f$= \psi_j({\bf r}_i)\f$. partial memory view of psiM_vgl
  DualPinnedValueMatrix_t psiM_temp;
  /// dpsiM(i,j) \f$= \nabla_i \psi_j({\bf r}_i)\f$. partial memory view of psiM_vgl
  DualPinnedGradMatrix dpsiM;
  /// d2psiM(i,j) \f$= \nabla_i^2 \psi_j({\bf r}_i)\f$. partial memory view of psiM_vgl
  DualPinnedValueMatrix_t d2psiM;

  /// Used for force computations
  GradMatrix_t grad_source_psiM, grad_lapl_source_psiM;
  HessMatrix_t grad_grad_source_psiM;

  GradMatrix_t phi_alpha_Minv, grad_phi_Minv;
  ValueMatrix_t lapl_phi_Minv;
  HessMatrix_t grad_phi_alpha_Minv;

  DualPinnedValueVector_t psiV;
  ValueVector_t psiV_host_view;
  DualPinnedGradVector dpsiV;
  GradVector_t dpsiV_host_view;
  DualPinnedValueVector_t d2psiV;
  ValueVector_t d2psiV_host_view;

  /// Delayed update engine 1 per walker.
  DET_ENGINE det_engine_;
  /**@}*/
  
  std::unique_ptr<DiracDeterminantBatchedMultiWalkerResource> mw_res_;

  LogValueType get_log_value() const { return LogValue; }

  static void mw_invertPsiM(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                            RefVector<DualPinnedValueMatrix_t>& logdetT_list,
                            RefVector<DualPinnedValueMatrix_t>& a_inv_list,
                            const std::vector<bool>& compute_mask);

  /// maximal number of delayed updates
  int ndelay;

  /// timers
  NewTimer &D2HTimer, &H2DTimer;

private:
  /// compute G adn L assuming psiMinv, dpsiM, d2psiM are ready for use
  void computeGL(ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L) const;

  /// Legacy single invert logdetT(psiM), result is stored in DDB object.
  void invertPsiM(DiracDeterminantBatchedMultiWalkerResource& mw_res, DualPinnedValueMatrix_t& logdetT, DualPinnedValueMatrix_t& a_inv);

  /// Resize all temporary arrays required for force computation.
  void resizeScratchObjectsForIonDerivs();

  //  friend class qmcplusplus::DiracDeterminantDetails;
  friend struct qmcplusplus::testing::SetupDiracDetResources;
  friend class qmcplusplus::testing::DiracDeterminantBatchedTest;
};



} // namespace qmcplusplus
#endif
