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

#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#include "QMCWaveFunctions/Fermion/MatrixUpdateOMPTarget.h"
#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
#include "QMCWaveFunctions/Fermion/MatrixDelayedUpdateCUDA.h"
#endif
#include "DualAllocatorAliases.hpp"
#include "WaveFunctionTypes.hpp"
#include "type_traits/complex_help.hpp"

namespace qmcplusplus
{

template<typename DET_ENGINE = MatrixUpdateOMPTarget<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>
class DiracDeterminantBatched : public DiracDeterminantBase
{
public:
  using WFT           = typename DET_ENGINE::WFT;
  using Value         = typename WFT::Value;
  using FullPrecValue = typename WFT::FullPrecValue;
  using PsiValue      = typename WFT::PsiValue;
  using LogValue      = typename WFT::LogValue;
  using Grad          = typename WFT::Grad;
  using Hess          = typename WFT::Hess;
  using Real          = typename WFT::Real;
  using FullPrecGrad  = TinyVector<FullPrecValue, DIM>;
  template<typename DT>
  using PinnedVector = Vector<DT, PinnedDualAllocator<DT>>;
  template<typename DT>
  using DualMatrix    = Matrix<DT, PinnedDualAllocator<DT>>;
  using DualVGLVector = VectorSoaContainer<Value, DIM + 2, PinnedDualAllocator<Value>>;

  struct DiracDeterminantBatchedMultiWalkerResource : public Resource
  {
    DiracDeterminantBatchedMultiWalkerResource() : Resource("DiracDeterminantBatched") {}
    DiracDeterminantBatchedMultiWalkerResource(const DiracDeterminantBatchedMultiWalkerResource&)
        : DiracDeterminantBatchedMultiWalkerResource()
    {}

    Resource* makeClone() const override { return new DiracDeterminantBatchedMultiWalkerResource(*this); }
    PinnedVector<LogValue> log_values;
    /// value, grads, laplacian of single-particle orbital for particle-by-particle update and multi walker [5][nw*norb]
    DualVGLVector phi_vgl_v;
    // multi walker of ratio
    std::vector<Value> ratios_local;
    // multi walker of grads
    std::vector<Grad> grad_new_local;
  };

  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   *@param last index of last particle
   *@param ndelay delayed update rank
   */
  DiracDeterminantBatched(std::shared_ptr<SPOSet>&& spos, int first, int last, int ndelay = 1);

  // copy constructor and assign operator disabled
  DiracDeterminantBatched(const DiracDeterminantBatched& s) = delete;
  DiracDeterminantBatched& operator=(const DiracDeterminantBatched& s) = delete;

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           std::vector<Value>& dlogpsi,
                           std::vector<Value>& dhpsioverpsi) override;

  void registerData(ParticleSet& P, WFBufferType& buf) override;

  LogValue updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;

  /** return the ratio only for the  iat-th partcle move
   * @param P current configuration
   * @param iat the particle thas is being moved
   */
  PsiValue ratio(ParticleSet& P, int iat) override;

  void mw_calcRatio(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    int iat,
                    std::vector<PsiValue>& ratios) const override;

  /** compute multiple ratios for a particle move
   */
  void evaluateRatios(const VirtualParticleSet& VP, std::vector<Value>& ratios) override;

  void mw_evaluateRatios(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                         const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                         std::vector<std::vector<Value>>& ratios) const override;

  PsiValue ratioGrad(ParticleSet& P, int iat, Grad& grad_iat) override;

  void mw_ratioGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    int iat,
                    std::vector<PsiValue>& ratios,
                    std::vector<Grad>& grad_new) const override;

  GradType evalGrad(ParticleSet& P, int iat) override;

  void mw_evalGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                   const RefVectorWithLeader<ParticleSet>& p_list,
                   int iat,
                   std::vector<Grad>& grad_now) const override;

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
  LogValue evaluateLog(const ParticleSet& P,
                       ParticleSet::ParticleGradient_t& G,
                       ParticleSet::ParticleLaplacian_t& L) override;

  void mw_evaluateLog(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                      const RefVectorWithLeader<ParticleSet>& p_list,
                      const RefVector<ParticleSet::ParticleGradient_t>& G_list,
                      const RefVector<ParticleSet::ParticleLaplacian_t>& L_list) const override;

  void recompute(const ParticleSet& P) override;

  LogValue evaluateGL(const ParticleSet& P,
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
  void acquireResource(ResourceCollection& collection,
                       const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override;
  void releaseResource(ResourceCollection& collection,
                       const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override;

  /** cloning function
   * @param tqp target particleset
   * @param spo spo set
   *
   * This interface is exposed only to SlaterDet and its derived classes
   * can overwrite to clone itself correctly.
   */
  DiracDeterminantBatched* makeCopy(std::shared_ptr<SPOSet>&& spo) const override;

  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<Value>& ratios) override;

  /// return  for testing
  auto& getPsiMinv() const { return psiMinv; }

  /** inverse transpose of psiM(j,i) \f$= \psi_j({\bf r}_i)\f$ memory view
   * Actual memory is owned by det_engine_
   * Only NumOrbitals x NumOrbitals subblock has meaningful data
   * The number of rows is equal to NumOrbitals
   * The number of columns in each row is padded to a multiple of QMC_SIMD_ALIGNMENT
   */
  Matrix<Value> psiMinv;

  /// memory for psiM, dpsiM and d2psiM. [5][norb*norb]
  DualVGLVector psiM_vgl;
  /// psiM(j,i) \f$= \psi_j({\bf r}_i)\f$. partial memory view of psiM_vgl
  Matrix<Value> psiM_temp;
  /// dpsiM(i,j) \f$= \nabla_i \psi_j({\bf r}_i)\f$. partial memory view of psiM_vgl
  Matrix<Grad> dpsiM;
  /// d2psiM(i,j) \f$= \nabla_i^2 \psi_j({\bf r}_i)\f$. partial memory view of psiM_vgl
  Matrix<Value> d2psiM;

  /// Used for force computations
  Matrix<Grad> grad_source_psiM, grad_lapl_source_psiM;
  Matrix<Hess> grad_grad_source_psiM;

  Matrix<Grad> phi_alpha_Minv, grad_phi_Minv;
  Matrix<Value> lapl_phi_Minv;
  Matrix<Hess> grad_phi_alpha_Minv;

  /// value of single-particle orbital for particle-by-particle update
  PinnedVector<Value> psiV;
  Vector<Value> psiV_host_view;
  Vector<Grad> dpsiV;
  Vector<Value> d2psiV;

  /// delayed update engine
  DET_ENGINE det_engine_;

  // psi(r')/psi(r) during a PbyP move
  PsiValue curRatio;

  std::unique_ptr<DiracDeterminantBatchedMultiWalkerResource> mw_res_;

private:
  ///reset the size: with the number of particles and number of orbtials
  void resize(int nel, int morb);

  /// compute G adn L assuming psiMinv, dpsiM, d2psiM are ready for use
  void computeGL(ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L) const;

  /// invert logdetT(psiM), result is in the engine.
  void invertPsiM(const Matrix<Value>& logdetT);

  static void mw_invertPsiM(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                            const RefVector<const Matrix<Value>>& logdetT_list);

  void mw_recompute(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    const std::vector<bool>& recompute) const override;

  // make this class unit tests friendly without the need of setup resources.
  void guardMultiWalkerRes()
  {
    if (!mw_res_)
    {
      std::cerr
          << "WARNING DiracDeterminantBatched : This message should not be seen in production (performance bug) runs "
             "but only unit tests (expected)."
          << std::endl;
      mw_res_ = std::make_unique<DiracDeterminantBatchedMultiWalkerResource>();
    }
  }

  /// Resize all temporary arrays required for force computation.
  void resizeScratchObjectsForIonDerivs();

  /// maximal number of delayed updates
  const int ndelay_;

  /// timers
  NewTimer &D2HTimer, &H2DTimer;
};

extern template class DiracDeterminantBatched<>;
#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
extern template class DiracDeterminantBatched<
    MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
#endif

} // namespace qmcplusplus
#endif
