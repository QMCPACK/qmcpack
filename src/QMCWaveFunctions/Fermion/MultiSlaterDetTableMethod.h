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


#ifndef QMCPLUSPLUS_MULTISLATERDETERMINANTFAST_ORBITAL_H
#define QMCPLUSPLUS_MULTISLATERDETERMINANTFAST_ORBITAL_H
#include <Configuration.h>
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Fermion/MultiDiracDeterminant.h"
#include "Utilities/TimerManager.h"
#include "Platforms/PinnedAllocator.h"
#include "OMPTarget/OMPallocator.hpp"

namespace qmcplusplus
{
/// CSF related dataset
struct CSFData
{
  using RealType  = WaveFunctionComponent::RealType;
  using ValueType = WaveFunctionComponent::ValueType;

  // coefficients of csfs, these are only used during optm
  std::vector<ValueType> coeffs;
  // number of dets per csf
  std::vector<size_t> dets_per_csf;
  // coefficient of csf expansion (smaller dimension)
  std::vector<RealType> expansion;
};

/** @ingroup WaveFunctionComponent
 *  @brief An AntiSymmetric WaveFunctionComponent composed of a linear combination of SlaterDeterminants.
 *
 *\f[
 *MS({\bf R}) = \sum_n c_n S_n({\bf R})
 *\f].
 *
 *The (S)ingle(P)article(O)rbitalSet template parameter is an
 *engine which fills in single-particle orbital terms.
 *
 \f[
 \frac{\nabla_i\Psi}{\Psi} = \frac{\sum_{n=1}^M c_n \nabla_i D_n}
 {\sum_{n=1}^M c_n D_n}
 \f]
 \f[
 \frac{\sum_{n=1}^M c_n S_n(\sum_{j=1}^N(\nabla_i
 S^{ij}_n({\bf r_i}))(S^{-1})^{ji}_n}{\sum_{n=1}^M c_n S_n}
 \f]
 The Laplacian
 \f[
 \frac{\nabla^2_i\Psi}{\Psi} = \frac{\sum_{n=1}^M c_n S_n(\sum_{j=1}^N
 (\nabla_i^2S^{ij}_n({\bf r_i}))(S^{-1})^{ji}_n}{\sum_{n=1}^M c_n S_n}
 \f]
 */
class MultiSlaterDetTableMethod : public WaveFunctionComponent
{
public:
  void registerTimers();
  NewTimer &RatioTimer, &MWRatioTimer, &OffloadRatioTimer, &OffloadGradTimer;
  NewTimer &EvalGradTimer, &MWEvalGradTimer, &RatioGradTimer, &MWRatioGradTimer;
  NewTimer &PrepareGroupTimer, &UpdateTimer, &AccRejTimer, &EvaluateTimer;

  template<typename DT>
  using OffloadPinnedAllocator = OMPallocator<DT, PinnedAlignedAllocator<DT>>;

  template<typename DT>
  using OffloadVector = Vector<DT, OffloadPinnedAllocator<DT>>;

  using SPOSetPtr   = SPOSet*;
  using IndexVector = OrbitalSetTraits<ValueType>::IndexVector;
  using ValueVector = OrbitalSetTraits<ValueType>::ValueVector;
  using GradVector  = OrbitalSetTraits<ValueType>::GradVector;
  using HessMatrix  = OrbitalSetTraits<ValueType>::HessMatrix;
  using ValueMatrix = OrbitalSetTraits<ValueType>::ValueMatrix;
  using HessType    = OrbitalSetTraits<ValueType>::HessType;
  using HessArray   = Array<HessType, 3>;
  using GGGType     = TinyVector<HessType, OHMMS_DIM>;
  using GGGVector   = Vector<GGGType>;
  using GGGMatrix   = Matrix<GGGType>;
  using Walker_t    = ParticleSet::Walker_t;


  ///constructor
  MultiSlaterDetTableMethod(ParticleSet& targetPtcl,
                            std::vector<std::unique_ptr<MultiDiracDeterminant>>&& dets,
                            bool use_pre_computing);

  ///destructor
  ~MultiSlaterDetTableMethod() override;

  void checkInVariables(opt_variables_type& active) override;
  void checkOutVariables(const opt_variables_type& active) override;
  void resetParameters(const opt_variables_type& active) override;
  void reportStatus(std::ostream& os) override;

  //builds orbital rotation parameters using MultiSlater member variables
  void buildOptVariables();

  LogValueType evaluate_vgl_impl(const ParticleSet& P,
                                 ParticleSet::ParticleGradient& g_tmp,
                                 ParticleSet::ParticleLaplacian& l_tmp);

  LogValueType evaluateLog(const ParticleSet& P,
                           ParticleSet::ParticleGradient& G,
                           ParticleSet::ParticleLaplacian& L) override;

  void prepareGroup(ParticleSet& P, int ig) override;

  GradType evalGrad(ParticleSet& P, int iat) override;
  //evalGrad, but returns the spin gradient as well
  GradType evalGradWithSpin(ParticleSet& P, int iat, ComplexType& spingrad) override;

  void mw_evalGrad(const RefVectorWithLeader<WaveFunctionComponent>& WFC_list,
                   const RefVectorWithLeader<ParticleSet>& P_list,
                   int iat,
                   std::vector<GradType>& grad_now) const override;

  void mw_ratioGrad(const RefVectorWithLeader<WaveFunctionComponent>& WFC_list,
                    const RefVectorWithLeader<ParticleSet>& P_list,
                    int iat,
                    std::vector<PsiValueType>& ratios,
                    std::vector<GradType>& grad_new) const override;

  void mw_calcRatio(const RefVectorWithLeader<WaveFunctionComponent>& WFC_list,
                    const RefVectorWithLeader<ParticleSet>& P_list,
                    int iat,
                    std::vector<PsiValueType>& ratios) const override;

  PsiValueType ratio(ParticleSet& P, int iat) override;
  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;
  //ratioGradWithSpin, but includes tthe spin gradient info
  PsiValueType ratioGradWithSpin(ParticleSet& P, int iat, GradType& grad_iat, ComplexType& spingrad_iat) override;

  void evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios) override;

  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios) override
  {
    // the base class routine may probably work, just never tested.
    // it can also be highly optimized with a specialized implementation.
    throw std::runtime_error(" Need to implement MultiSlaterDetTableMethod::evaluateRatiosAlltoOne. \n");
  }

  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override;
  void restore(int iat) override;

  void registerData(ParticleSet& P, WFBufferType& buf) override;
  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;
  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;

  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tqp) const override;
  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           std::vector<ValueType>& dlogpsi,
                           std::vector<ValueType>& dhpsioverpsi) override;

  void evaluateDerivativesWF(ParticleSet& P,
                             const opt_variables_type& optvars,
                             std::vector<ValueType>& dlogpsi) override;

  /** initialize a few objects and states by the builder
   * YL: it should be part of the constructor. It cannot be added to the constructor
   * because the constructor is used by makeClone. The right way of fix needs:
   *  1. implement a copy constructor and reroute makeClone to it.
   *  2. merge initialize() into the constructor.
   */
  void initialize(std::unique_ptr<std::vector<std::vector<size_t>>> C2node_in,
                  std::unique_ptr<std::vector<ValueType>> C_in,
                  std::unique_ptr<opt_variables_type> myVars_in,
                  std::unique_ptr<CSFData> csf_data_in,
                  bool optimizable,
                  bool CI_optimizable);

private:
  //get Det ID. It should be consistent with particle group id within the particle set.
  inline int getDetID(const int iat) const
  {
    int id = 0;
    while (iat > Last[id])
      id++;
    return id;
  }

  /** an implementation shared by evalGrad and ratioGrad. Use precomputed data
   * @param newpos to distinguish evalGrad(false) ratioGrad(true)
   */
  PsiValueType evalGrad_impl(ParticleSet& P, int iat, bool newpos, GradType& g_at);
  /// multi walker version of evalGrad_impl
  static void mw_evalGrad_impl(const RefVectorWithLeader<WaveFunctionComponent>& WFC_list,
                               const RefVectorWithLeader<ParticleSet>& P_list,
                               int iat,
                               bool newpos,
                               std::vector<GradType>& grad_now,
                               std::vector<PsiValueType>& psi_list);

  /** an implementation shared by evalGrad and ratioGrad. No use of precomputed data
   * @param newpos to distinguish evalGrad(false) ratioGrad(true)
   */
  PsiValueType evalGrad_impl_no_precompute(ParticleSet& P, int iat, bool newpos, GradType& g_at);

  //implemtation for evalGradWithSpin
  PsiValueType evalGradWithSpin_impl(ParticleSet& P, int iat, bool newpos, GradType& g_at, ComplexType& sg_at);
  //implemtation for evalGradWithSpin with no precomputation
  PsiValueType evalGradWithSpin_impl_no_precompute(ParticleSet& P,
                                                   int iat,
                                                   bool newpos,
                                                   GradType& g_at,
                                                   ComplexType& sg_at);

  // an implementation of ratio. Use precomputed data
  PsiValueType ratio_impl(ParticleSet& P, int iat);
  // an implementation of ratio. No use of precomputed data
  PsiValueType ratio_impl_no_precompute(ParticleSet& P, int iat);

  /** precompute C_otherDs for a given particle group
   * @param P a particle set
   * @param ig group id
   */
  void precomputeC_otherDs(const ParticleSet& P, int ig);

  void evaluateMultiDiracDeterminantDerivatives(ParticleSet& P,
                                                const opt_variables_type& optvars,
                                                std::vector<ValueType>& dlogpsi,
                                                std::vector<ValueType>& dhpsioverpsi);

  void evaluateMultiDiracDeterminantDerivativesWF(ParticleSet& P,
                                                  const opt_variables_type& optvars,
                                                  std::vector<ValueType>& dlogpsi);

  /// determinant collection
  std::vector<std::unique_ptr<MultiDiracDeterminant>> Dets;

  /** map determinant in linear combination to unique det list
   * map global det id to unique det id. [spin, global det id] = unique det id
   */
  std::shared_ptr<std::vector<std::vector<size_t>>> C2node;
  /// CI coefficients
  std::shared_ptr<std::vector<ValueType>> C;
  /// if true, the CI coefficients are optimized
  bool CI_Optimizable;
  //optimizable variable is shared with the clones
  std::shared_ptr<opt_variables_type> myVars;

  /// CSF data set. If nullptr, not using CSF
  std::shared_ptr<CSFData> csf_data_;

  ///the last particle of each group
  std::vector<int> Last;
  ///use pre-compute (fast) algorithm
  const bool use_pre_computing_;

  /// current psi over ref single det
  PsiValueType psi_ratio_to_ref_det_;
  /// new psi over new ref single det when one particle is moved
  PsiValueType new_psi_ratio_to_new_ref_det_;

  size_t ActiveSpin;
  PsiValueType curRatio;

  /// C_n x D^1_n x D^2_n ... D^3_n with one D removed. Summed by group. [spin, unique det id]
  //std::vector<Vector<ValueType, OffloadPinnedAllocator<ValueType>>> C_otherDs;
  std::vector<OffloadVector<ValueType>> C_otherDs;
  /// a collection of device pointers of multiple walkers fused for fast H2D transfer.
  OffloadVector<const ValueType*> C_otherDs_ptr_list;
  OffloadVector<const ValueType*> det_value_ptr_list;

  ParticleSet::ParticleGradient myG, myG_temp;
  ParticleSet::ParticleLaplacian myL, myL_temp;
  std::vector<ValueVector> laplSum;

  // temporary storage for evaluateDerivatives
  ParticleSet::ParticleGradient gmPG;
  std::vector<Matrix<RealType>> dpsia, dLa;
  std::vector<Array<GradType, OHMMS_DIM>> dGa;
};

} // namespace qmcplusplus
#endif
