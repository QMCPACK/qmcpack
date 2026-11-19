//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SLATERDETERMINANT_WITHBASE_H
#define QMCPLUSPLUS_SLATERDETERMINANT_WITHBASE_H

#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/SPOSet.h"

namespace qmcplusplus
{
class DiracDeterminantBase;
class TWFFastDerivWrapper;

class SlaterDet : public WaveFunctionComponent
{
public:
  /**  constructor
   * @param targetPtcl target Particleset
   */
  SlaterDet(ParticleSet& targetPtcl,
            std::vector<std::unique_ptr<SPOSet>>&& sposets,
            std::vector<std::unique_ptr<DiracDeterminantBase>>&& dets,
            const std::string& class_name = "SlaterDet");

  ///destructor
  ~SlaterDet() override;

  std::string getClassName() const override { return "SlaterDet"; }

  bool isFermionic() const final { return true; }
  bool isOptimizable() const override;

  void extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs) override;

  void checkOutVariables(const OptVariables& active) override;

  void registerTWFFastDerivWrapper(const ParticleSet& P, TWFFastDerivWrapper& twf) const override;

  LogValue evaluateLog(const ParticleSet& P,
                       ParticleSet::ParticleGradient& G,
                       ParticleSet::ParticleLaplacian& L) override;

  void mw_evaluateLog(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                      const RefVectorWithLeader<ParticleSet>& p_list,
                      const RefVector<ParticleSet::ParticleGradient>& G_list,
                      const RefVector<ParticleSet::ParticleLaplacian>& L_list) const override;

  LogValue evaluateGL(const ParticleSet& P,
                      ParticleSet::ParticleGradient& G,
                      ParticleSet::ParticleLaplacian& L,
                      bool fromscratch) override;

  void mw_evaluateGL(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                     const RefVectorWithLeader<ParticleSet>& p_list,
                     const RefVector<ParticleSet::ParticleGradient>& G_list,
                     const RefVector<ParticleSet::ParticleLaplacian>& L_list,
                     bool fromscratch) const override;

  void recompute(const ParticleSet& P) override;

  void mw_recompute(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    const std::vector<bool>& recompute) const override;

  void evaluateHessian(ParticleSet& P, HessVector& grad_grad_psi) override;

  void registerData(ParticleSet& P, WFBufferType& buf) override;

  LogValue updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;

  void createResource(ResourceCollection& collection) const override;

  void acquireResource(ResourceCollection& collection,
                       const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override;

  void releaseResource(ResourceCollection& collection,
                       const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override;

  void evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios) override;

  void evaluateSpinorRatios(const VirtualParticleSet& VP,
                            const std::pair<ValueVector, ValueVector>& spinor_multiplier,
                            std::vector<ValueType>& ratios) override;

  void evaluateDerivRatios(const VirtualParticleSet& VP,
                           const OptVariables& optvars,
                           std::vector<ValueType>& ratios,
                           Matrix<ValueType>& dratios) override;

  void evaluateSpinorDerivRatios(const VirtualParticleSet& VP,
                                 const std::pair<ValueVector, ValueVector>& spinor_multiplier,
                                 const OptVariables& optvars,
                                 std::vector<ValueType>& ratios,
                                 Matrix<ValueType>& dratios) override;

  void mw_evaluateRatios(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                         const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                         std::vector<std::vector<ValueType>>& ratios) const override;

  void mw_evaluateSpinorRatios(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                               const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                               const RefVector<std::pair<ValueVector, ValueVector>>& spinor_multiplier_list,
                               std::vector<std::vector<ValueType>>& ratios) const final;

  PsiValue ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;

  PsiValue ratioGradWithSpin(ParticleSet& P, int iat, GradType& grad_iat, ComplexType& spingrad_iat) override;

  void mw_ratioGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    int iat,
                    std::vector<PsiValue>& ratios,
                    std::vector<GradType>& grad_now) const override;

  void mw_ratioGradWithSpin(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                            const RefVectorWithLeader<ParticleSet>& p_list,
                            int iat,
                            std::vector<PsiValue>& ratios,
                            std::vector<GradType>& grad_now,
                            std::vector<ComplexType>& spingrad_now) const override;

  GradType evalGrad(ParticleSet& P, int iat) override;

  GradType evalGradWithSpin(ParticleSet& P, int iat, ComplexType& spingrad) override;

  void mw_evalGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                   const RefVectorWithLeader<ParticleSet>& p_list,
                   int iat,
                   std::vector<GradType>& grad_now) const override;

  void mw_evalGradWithSpin(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                           const RefVectorWithLeader<ParticleSet>& p_list,
                           int iat,
                           std::vector<GradType>& grad_now,
                           std::vector<ComplexType>& spingrad_now) const override;

  GradType evalGradSource(ParticleSet& P, ParticleSet& src, int iat) override;

  GradType evalGradSource(ParticleSet& P,
                          ParticleSet& src,
                          int iat,
                          TinyVector<ParticleSet::ParticleGradient, OHMMS_DIM>& grad_grad,
                          TinyVector<ParticleSet::ParticleLaplacian, OHMMS_DIM>& lapl_grad) override;

  void restore(int iat) override;

  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override;

  void mw_accept_rejectMove(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                            const RefVectorWithLeader<ParticleSet>& p_list,
                            int iat,
                            const std::vector<bool>& isAccepted,
                            bool safe_to_delay = false) const override;

  void completeUpdates() override;

  void mw_completeUpdates(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override;

  PsiValue ratio(ParticleSet& P, int iat) override;

  void mw_calcRatio(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    int iat,
                    std::vector<PsiValue>& ratios) const override;

  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tqp) const override;

  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios) override;

  void evaluateDerivatives(ParticleSet& P,
                           const OptVariables& active,
                           Vector<ValueType>& dlogpsi,
                           Vector<ValueType>& dhpsioverpsi) override;

  void evaluateDerivativesWF(ParticleSet& P, const OptVariables& active, Vector<ValueType>& dlogpsi) override;

  ///return the total number of Dirac determinants
  inline int getNumDets() const { return Dets.size(); }
  ///return the i-th determinant
  DiracDeterminantBase& getDet(int i);
  const DiracDeterminantBase& getDet(int i) const;
  ///return the sposet of the i-th determinant
  SPOSet& getPhi(int i = 0);

private:
  //get Det ID
  inline int getDetID(const int iat) const
  {
    int id = 0;
    while (iat > Last[id])
      id++;
    return id;
  }

  // helper function for extracting a list of WaveFunctionComponent from a list of TrialWaveFunction
  RefVectorWithLeader<WaveFunctionComponent> extract_DetRef_list(
      const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
      int det_id) const;

  ///the last particle of each group
  std::vector<int> Last;

  ///container for the unique SPOSets
  const std::vector<std::unique_ptr<SPOSet>> sposets_;

  ///container for the DiracDeterminants
  const std::vector<std::unique_ptr<DiracDeterminantBase>> Dets;
};
} // namespace qmcplusplus
#endif
