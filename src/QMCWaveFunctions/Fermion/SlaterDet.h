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
#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#include <map>

namespace qmcplusplus
{
// NOTE NOTE NOTE
// template<bool backflow>
//  class SlaterDet: public WaveFunctionComponent {}
//     then change SlaterDet to SlaterDet<false>
//     and SlaterDeterminantWithBackflow to SlaterDet<true>
//     and remove all virtuals and inline them

class SlaterDet : public WaveFunctionComponent
{
public:
  using Determinant_t = DiracDeterminantBase;
  ///container for the DiracDeterminants
  const std::vector<std::unique_ptr<Determinant_t>> Dets;

  /**  constructor
   * @param targetPtcl target Particleset
   */
  SlaterDet(ParticleSet& targetPtcl, std::vector<std::unique_ptr<Determinant_t>> dets, const std::string& class_name = "SlaterDet");

  ///destructor
  ~SlaterDet() override;

  void checkInVariables(opt_variables_type& active) override;

  void checkOutVariables(const opt_variables_type& active) override;

  ///reset all the Dirac determinants, Optimizable is true
  void resetParameters(const opt_variables_type& optVariables) override;

  void reportStatus(std::ostream& os) override;

  LogValueType evaluateLog(const ParticleSet& P,
                           ParticleSet::ParticleGradient_t& G,
                           ParticleSet::ParticleLaplacian_t& L) override;

  void mw_evaluateLog(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                      const RefVectorWithLeader<ParticleSet>& p_list,
                      const RefVector<ParticleSet::ParticleGradient_t>& G_list,
                      const RefVector<ParticleSet::ParticleLaplacian_t>& L_list) const override;

  LogValueType evaluateGL(const ParticleSet& P,
                          ParticleSet::ParticleGradient_t& G,
                          ParticleSet::ParticleLaplacian_t& L,
                          bool fromscratch) override;

  void mw_evaluateGL(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                     const RefVectorWithLeader<ParticleSet>& p_list,
                     const RefVector<ParticleSet::ParticleGradient_t>& G_list,
                     const RefVector<ParticleSet::ParticleLaplacian_t>& L_list,
                     bool fromscratch) const override;

  void recompute(const ParticleSet& P) override;

  void mw_recompute(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    const std::vector<bool>& recompute) const override;

  void evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi) override;

  ///return the total number of Dirac determinants
  inline int size() const { return Dets.size(); }

  void registerData(ParticleSet& P, WFBufferType& buf) override;

  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;

  void createResource(ResourceCollection& collection) const override;

  void acquireResource(ResourceCollection& collection, const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override;

  void releaseResource(ResourceCollection& collection, const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override;

  inline void evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios) override
  {
    return Dets[getDetID(VP.refPtcl)]->evaluateRatios(VP, ratios);
  }

  inline void mw_evaluateRatios(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                                std::vector<std::vector<ValueType>>& ratios) const override
  {
    if (wfc_list.size())
    {
      // assuming all the VP.refPtcl are identical
      const int det_id = getDetID(vp_list[0].refPtcl);
      Dets[det_id]->mw_evaluateRatios(extract_DetRef_list(wfc_list, det_id), vp_list, ratios);
    }
  }

  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;

  PsiValueType ratioGradWithSpin(ParticleSet& P, int iat, GradType& grad_iat, ComplexType& spingrad_iat) override;

  void mw_ratioGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    int iat,
                    std::vector<PsiValueType>& ratios,
                    std::vector<GradType>& grad_now) const override;

  GradType evalGrad(ParticleSet& P, int iat) override { return Dets[getDetID(iat)]->evalGrad(P, iat); }

  GradType evalGradWithSpin(ParticleSet& P, int iat, ComplexType& spingrad) override
  {
    return Dets[getDetID(iat)]->evalGradWithSpin(P, iat, spingrad);
  }

  void mw_evalGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                   const RefVectorWithLeader<ParticleSet>& p_list,
                   int iat,
                   std::vector<GradType>& grad_now) const override
  {
    const int det_id = getDetID(iat);
    Dets[det_id]->mw_evalGrad(extract_DetRef_list(wfc_list, det_id), p_list, iat, grad_now);
  }

  GradType evalGradSource(ParticleSet& P, ParticleSet& src, int iat) override
  {
    GradType G = GradType();
    for (int iz = 0; iz < size(); iz++)
      G += Dets[iz]->evalGradSource(P, src, iat);
    return G;
  }

  GradType evalGradSource(ParticleSet& P,
                          ParticleSet& src,
                          int iat,
                          TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM>& grad_grad,
                          TinyVector<ParticleSet::ParticleLaplacian_t, OHMMS_DIM>& lapl_grad) override
  {
    GradType G = GradType();
    for (int iz = 0; iz < size(); iz++)
      G += Dets[iz]->evalGradSource(P, src, iat, grad_grad, lapl_grad);
    return G;
  }

  inline void restore(int iat) override { return Dets[getDetID(iat)]->restore(iat); }

  inline void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override
  {
    Dets[getDetID(iat)]->acceptMove(P, iat, safe_to_delay);

    log_value_ = 0.0;
    for (int i = 0; i < Dets.size(); ++i)
      log_value_ += Dets[i]->get_log_value();
  }

  void mw_accept_rejectMove(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                            const RefVectorWithLeader<ParticleSet>& p_list,
                            int iat,
                            const std::vector<bool>& isAccepted,
                            bool safe_to_delay = false) const override
  {
    constexpr LogValueType czero(0);

    // This log_value_ is in the slater determinant, it's still around but not consistent anymore with the
    // sum of the log_values in its determinants.  Caching the state seems like a bad call, but the wfc base class
    // having log_value_ as a data member asks for this sort of consistency issue when wfc can contain wfc.
    for (int iw = 0; iw < wfc_list.size(); iw++)
      if (isAccepted[iw])
        wfc_list.getCastedElement<SlaterDet>(iw).log_value_ = czero;

    for (int i = 0; i < Dets.size(); ++i)
    {
      const auto Det_list(extract_DetRef_list(wfc_list, i));

      if (i == getDetID(iat))
        Dets[i]->mw_accept_rejectMove(Det_list, p_list, iat, isAccepted, safe_to_delay);

      for (int iw = 0; iw < wfc_list.size(); iw++)
        if (isAccepted[iw])
          wfc_list.getCastedElement<SlaterDet>(iw).log_value_ += Det_list[iw].get_log_value();
    }
  }

  void completeUpdates() override
  {
    for (int i = 0; i < Dets.size(); i++)
      Dets[i]->completeUpdates();
  }

  void mw_completeUpdates(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override
  {
    for (int i = 0; i < Dets.size(); i++)
      Dets[i]->mw_completeUpdates(extract_DetRef_list(wfc_list, i));
  }

  inline PsiValueType ratio(ParticleSet& P, int iat) override { return Dets[getDetID(iat)]->ratio(P, iat); }

  void mw_calcRatio(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    int iat,
                    std::vector<PsiValueType>& ratios) const override
  {
    const int det_id = getDetID(iat);
    Dets[det_id]->mw_calcRatio(extract_DetRef_list(wfc_list, det_id), p_list, iat, ratios);
  }

  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tqp) const override;

  virtual SPOSetPtr getPhi(int i = 0) { return Dets[i]->getPhi(); }

  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios) override;

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           std::vector<ValueType>& dlogpsi,
                           std::vector<ValueType>& dhpsioverpsi) override
  {
    // First zero out values, since each determinant only adds on
    // its contribution (i.e. +=) , rather than setting the value
    // (i.e. =)
    for (int k = 0; k < myVars.size(); ++k)
    {
      int kk = myVars.where(k);
      if (kk >= 0)
        dlogpsi[kk] = dhpsioverpsi[kk] = 0.0;
    }
    // Now add on contribution from each determinant to the derivatives
    for (int i = 0; i < Dets.size(); i++)
      Dets[i]->evaluateDerivatives(P, active, dlogpsi, dhpsioverpsi);
  }

  void evaluateGradDerivatives(const ParticleSet::ParticleGradient_t& G_in,
                               std::vector<ValueType>& dgradlogpsi) override
  {
    for (int i = 0; i < Dets.size(); i++)
      Dets[i]->evaluateGradDerivatives(G_in, dgradlogpsi);
  }

#ifdef QMC_CUDA
  /////////////////////////////////////////////////////
  // Functions for vectorized evaluation and updates //
  /////////////////////////////////////////////////////
  void recompute(MCWalkerConfiguration& W, bool firstTime) override
  {
    for (int id = 0; id < Dets.size(); id++)
      Dets[id]->recompute(W, firstTime);
  }

  void reserve(PointerPool<gpu::device_vector<CTS::ValueType>>& pool, int kblocksize = 1) override
  {
    for (int id = 0; id < Dets.size(); id++)
      Dets[id]->reserve(pool, kblocksize);
  }

  void addLog(MCWalkerConfiguration& W, std::vector<RealType>& logPsi) override
  {
    for (int id = 0; id < Dets.size(); id++)
      Dets[id]->addLog(W, logPsi);
  }

  void ratio(MCWalkerConfiguration& W,
             int iat,
             std::vector<ValueType>& psi_ratios,
             std::vector<GradType>& grad,
             std::vector<ValueType>& lapl) override
  {
    Dets[getDetID(iat)]->ratio(W, iat, psi_ratios, grad, lapl);
  }

  void det_lookahead(MCWalkerConfiguration& W,
                     std::vector<ValueType>& psi_ratios,
                     std::vector<GradType>& grad,
                     std::vector<ValueType>& lapl,
                     int iat,
                     int k,
                     int kd,
                     int nw) override
  {
    Dets[getDetID(iat)]->det_lookahead(W, psi_ratios, grad, lapl, iat, k, kd, nw);
  }
  void calcRatio(MCWalkerConfiguration& W,
                 int iat,
                 std::vector<ValueType>& psi_ratios,
                 std::vector<GradType>& grad,
                 std::vector<ValueType>& lapl) override
  {
    Dets[getDetID(iat)]->calcRatio(W, iat, psi_ratios, grad, lapl);
  }

  void addRatio(MCWalkerConfiguration& W,
                int iat,
                int k,
                std::vector<ValueType>& psi_ratios,
                std::vector<GradType>& grad,
                std::vector<ValueType>& lapl) override
  {
    Dets[getDetID(iat)]->addRatio(W, iat, k, psi_ratios, grad, lapl);
  }

  void ratio(std::vector<Walker_t*>& walkers,
             std::vector<int>& iatList,
             std::vector<PosType>& rNew,
             std::vector<ValueType>& psi_ratios,
             std::vector<GradType>& grad,
             std::vector<ValueType>& lapl) override;

  void calcGradient(MCWalkerConfiguration& W, int iat, int k, std::vector<GradType>& grad) override
  {
    Dets[getDetID(iat)]->calcGradient(W, iat, k, grad);
  }

  void addGradient(MCWalkerConfiguration& W, int iat, std::vector<GradType>& grad) override
  {
    Dets[getDetID(iat)]->addGradient(W, iat, grad);
  }

  void update(MCWalkerConfiguration* W,
              std::vector<Walker_t*>& walkers,
              int iat,
              std::vector<bool>* acc,
              int k) override
  {
    Dets[getDetID(iat)]->update(W, walkers, iat, acc, k);
  }

  void update(const std::vector<Walker_t*>& walkers, const std::vector<int>& iatList) override;

  void gradLapl(MCWalkerConfiguration& W, GradMatrix_t& grads, ValueMatrix_t& lapl) override
  {
    for (int id = 0; id < Dets.size(); id++)
      Dets[id]->gradLapl(W, grads, lapl);
  }

  void NLratios(MCWalkerConfiguration& W,
                std::vector<NLjob>& jobList,
                std::vector<PosType>& quadPoints,
                std::vector<ValueType>& psi_ratios) override
  {
    for (int id = 0; id < Dets.size(); id++)
      Dets[id]->NLratios(W, jobList, quadPoints, psi_ratios);
  }
#endif

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
      int det_id) const
  {
    RefVectorWithLeader<WaveFunctionComponent> Det_list(*wfc_list.getCastedLeader<SlaterDet>().Dets[det_id]);
    Det_list.reserve(wfc_list.size());
    for (WaveFunctionComponent& wfc : wfc_list)
      Det_list.push_back(*static_cast<SlaterDet&>(wfc).Dets[det_id]);
    return Det_list;
  }

  ///the last particle of each group
  std::vector<int> Last;
};
} // namespace qmcplusplus
#endif
