//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file DiracDeterminant.h
 * @brief Declaration of DiracDeterminant with a S(ingle)P(article)O(rbital)Set
 */
#ifndef QMCPLUSPLUS_DIRACDETERMINANT_H
#define QMCPLUSPLUS_DIRACDETERMINANT_H

#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#include "QMCWaveFunctions/Fermion/DelayedUpdate.h"
#if defined(ENABLE_CUDA)
#include "QMCWaveFunctions/Fermion/DelayedUpdateCUDA.h"
#endif

namespace qmcplusplus
{
template<typename DU_TYPE = DelayedUpdate<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>
class DiracDeterminant : public DiracDeterminantBase
{
protected:
  int ndelay;

public:
  using ValueVector_t = SPOSet::ValueVector_t;
  using ValueMatrix_t = SPOSet::ValueMatrix_t;
  using GradVector_t  = SPOSet::GradVector_t;
  using GradMatrix_t  = SPOSet::GradMatrix_t;
  using HessMatrix_t  = SPOSet::HessMatrix_t;
  using HessVector_t  = SPOSet::HessVector_t;
  using HessType      = SPOSet::HessType;

  using mValueType = QMCTraits::QTFull::ValueType;
  using mGradType  = TinyVector<mValueType, DIM>;

  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   */
  DiracDeterminant(std::shared_ptr<SPOSet>&& spos, int first = 0);

  // copy constructor and assign operator disabled
  DiracDeterminant(const DiracDeterminant& s) = delete;
  DiracDeterminant& operator=(const DiracDeterminant& s) = delete;

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

  void updateAfterSweep(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L);

  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;

  /** return the ratio only for the  iat-th partcle move
   * @param P current configuration
   * @param iat the particle thas is being moved
   */
  PsiValueType ratio(ParticleSet& P, int iat) override;

  //Ye: TODO, good performance needs batched SPO evaluation.
  //void mw_calcRatio(const std::vector<WaveFunctionComponent*>& wfc_list,
  //                  const std::vector<ParticleSet*>& p_list,
  //                  int iat,
  //                  std::vector<PsiValueType>& ratios) override;

  /** compute multiple ratios for a particle move
   */
  void evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios) override;

  void mw_evaluateRatios(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                         const RefVector<const VirtualParticleSet>& vp_list,
                         std::vector<std::vector<ValueType>>& ratios) const override;

  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;

  PsiValueType ratioGradWithSpin(ParticleSet& P, int iat, GradType& grad_iat, ComplexType& spingrad) override final;

  void mw_ratioGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    int iat,
                    std::vector<PsiValueType>& ratios,
                    std::vector<GradType>& grad_new) const override;

  GradType evalGrad(ParticleSet& P, int iat) override;

  GradType evalGradWithSpin(ParticleSet& P, int iat, ComplexType& spingrad) override final;

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
                            bool safe_to_delay = false) const override
  {
    for (int iw = 0; iw < wfc_list.size(); iw++)
      if (isAccepted[iw])
        wfc_list[iw].acceptMove(p_list[iw], iat, safe_to_delay);
      else
        wfc_list[iw].restore(iat);
  }

  void completeUpdates() override;

  void mw_completeUpdates(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override
  {
    for (int iw = 0; iw < wfc_list.size(); iw++)
      wfc_list[iw].completeUpdates();
  }

  /** move was rejected. copy the real container to the temporary to move on
   */
  void restore(int iat) override;

  ///evaluate log of a determinant for a particle set
  LogValueType evaluateLog(ParticleSet& P,
                           ParticleSet::ParticleGradient_t& G,
                           ParticleSet::ParticleLaplacian_t& L) override;

  //Ye: TODO, good performance needs batched SPO evaluation.
  //void mw_evaluateLog(const std::vector<WaveFunctionComponent*>& wfc_list,
  //                    const std::vector<ParticleSet*>& p_list,
  //                    const std::vector<ParticleSet::ParticleGradient_t*>& G_list,
  //                    const std::vector<ParticleSet::ParticleLaplacian_t*>& L_list) override;

  void recompute(ParticleSet& P) override;

  LogValueType evaluateGL(ParticleSet& P,
                          ParticleSet::ParticleGradient_t& G,
                          ParticleSet::ParticleLaplacian_t& L,
                          bool fromscratch) override;

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
  DiracDeterminant* makeCopy(std::shared_ptr<SPOSet>&& spo) const override;

  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios) override;

#ifndef NDEBUG
  /// return  for testing
  ValueMatrix_t& getPsiMinv() override { return psiM; }
#else
  ValueMatrix_t& getPsiMinv() { return psiM; }
#endif

  /// psiM(j,i) \f$= \psi_j({\bf r}_i)\f$
  ValueMatrix_t psiM_temp;

  /// inverse transpose of psiM(j,i) \f$= \psi_j({\bf r}_i)\f$
  ValueMatrix_t psiM;

  /// temporary container for testing
  ValueMatrix_t psiMinv;

  /// dpsiM(i,j) \f$= \nabla_i \psi_j({\bf r}_i)\f$
  GradMatrix_t dpsiM;

  /// d2psiM(i,j) \f$= \nabla_i^2 \psi_j({\bf r}_i)\f$
  ValueMatrix_t d2psiM;

  /// Used for force computations
  GradMatrix_t grad_source_psiM, grad_lapl_source_psiM;
  HessMatrix_t grad_grad_source_psiM;

  GradMatrix_t phi_alpha_Minv, grad_phi_Minv;
  ValueMatrix_t lapl_phi_Minv;
  HessMatrix_t grad_phi_alpha_Minv;

  /// value of single-particle orbital for particle-by-particle update
  ValueVector_t psiV;
  ValueVector_t dspin_psiV;
  GradVector_t dpsiV;
  ValueVector_t d2psiV;

  /// delayed update engine
  DU_TYPE updateEng;

  /// the row of up-to-date inverse matrix
  ValueVector_t invRow;

  /** row id correspond to the up-to-date invRow. [0 norb), invRow is ready; -1, invRow is not valid.
   *  This id is set after calling getInvRow indicating invRow has been prepared for the invRow_id row
   *  ratioGrad checks if invRow_id is consistent. If not, invRow needs to be recomputed.
   *  acceptMove and completeUpdates mark invRow invalid by setting invRow_id to -1
   */
  int invRow_id;

  PsiValueType curRatio;
  ValueType* FirstAddressOfdV;
  ValueType* LastAddressOfdV;

private:
  /// invert psiM or its copies
  void invertPsiM(const ValueMatrix_t& logdetT, ValueMatrix_t& invMat);

  /// Resize all temporary arrays required for force computation.
  void resizeScratchObjectsForIonDerivs();

  /// internal function computing ratio and gradients after computing the SPOs, used by ratioGrad.
  PsiValueType ratioGrad_compute(int iat, GradType& grad_iat);
};

extern template class DiracDeterminant<>;
#if defined(ENABLE_CUDA)
extern template class DiracDeterminant<DelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
#endif

} // namespace qmcplusplus
#endif
