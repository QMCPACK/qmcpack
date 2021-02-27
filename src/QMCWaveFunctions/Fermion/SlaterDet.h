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
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
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
  typedef DiracDeterminantBase Determinant_t;
  ///container for the DiracDeterminants
  std::vector<std::unique_ptr<Determinant_t>> Dets;

  /**  constructor
   * @param targetPtcl target Particleset
   */
  SlaterDet(ParticleSet& targetPtcl, const std::string& class_name = "SlaterDet");

  ///destructor
  ~SlaterDet();

  ///add a new DiracDeterminant to the list of determinants
  virtual void add(Determinant_t* det, int ispin);

  ///set BF pointers
  virtual void setBF(BackflowTransformation* BFTrans) {}

  virtual void checkInVariables(opt_variables_type& active) override;

  virtual void checkOutVariables(const opt_variables_type& active) override;

  ///reset all the Dirac determinants, Optimizable is true
  virtual void resetParameters(const opt_variables_type& optVariables) override;

  void reportStatus(std::ostream& os) override;

  virtual LogValueType evaluateLog(ParticleSet& P,
                                   ParticleSet::ParticleGradient_t& G,
                                   ParticleSet::ParticleLaplacian_t& L) override;

  virtual void mw_evaluateLog(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                              const RefVectorWithLeader<ParticleSet>& p_list,
                              const RefVector<ParticleSet::ParticleGradient_t>& G_list,
                              const RefVector<ParticleSet::ParticleLaplacian_t>& L_list) const override;

  virtual LogValueType evaluateGL(ParticleSet& P,
                                  ParticleSet::ParticleGradient_t& G,
                                  ParticleSet::ParticleLaplacian_t& L,
                                  bool fromscratch) override;

  virtual void mw_evaluateGL(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                             const RefVectorWithLeader<ParticleSet>& p_list,
                             const RefVector<ParticleSet::ParticleGradient_t>& G_list,
                             const RefVector<ParticleSet::ParticleLaplacian_t>& L_list,
                             bool fromscratch) const override;

  virtual void recompute(ParticleSet& P) override;

  virtual void evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi) override;

  ///return the total number of Dirac determinants
  inline int size() const { return Dets.size(); }

  virtual void registerData(ParticleSet& P, WFBufferType& buf) override;

  virtual LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;

  virtual void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;

  void createResource(ResourceCollection& collection) const override;

  void acquireResource(ResourceCollection& collection) override;

  void releaseResource(ResourceCollection& collection) override;

  virtual inline void evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios) override
  {
    return Dets[getDetID(VP.refPtcl)]->evaluateRatios(VP, ratios);
  }

  virtual inline void mw_evaluateRatios(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                        const RefVector<const VirtualParticleSet>& vp_list,
                                        std::vector<std::vector<ValueType>>& ratios) const override
  {
    // assuming all the VP.refPtcl are identical
    const int det_id = getDetID(vp_list[0].get().refPtcl);
    return Dets[det_id]->mw_evaluateRatios(extract_DetRef_list(wfc_list, det_id), vp_list, ratios);
  }

  virtual PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;
  virtual void ratioGradAsync(ParticleSet& P, int iat, PsiValueType& ratio, GradType& grad_iat) override;

  virtual PsiValueType ratioGradWithSpin(ParticleSet& P,
                                         int iat,
                                         GradType& grad_iat,
                                         ComplexType& spingrad_iat) override;

  virtual void mw_ratioGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                            const RefVectorWithLeader<ParticleSet>& p_list,
                            int iat,
                            std::vector<PsiValueType>& ratios,
                            std::vector<GradType>& grad_now) const override;

  void mw_ratioGradAsync(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                         const RefVectorWithLeader<ParticleSet>& p_list,
                         int iat,
                         std::vector<PsiValueType>& ratios,
                         std::vector<GradType>& grad_now) const override;

  virtual GradType evalGrad(ParticleSet& P, int iat) override { return Dets[getDetID(iat)]->evalGrad(P, iat); }

  virtual GradType evalGradWithSpin(ParticleSet& P, int iat, ComplexType& spingrad) override
  {
    return Dets[getDetID(iat)]->evalGradWithSpin(P, iat, spingrad);
  }

  virtual void mw_evalGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                           const RefVectorWithLeader<ParticleSet>& p_list,
                           int iat,
                           std::vector<GradType>& grad_now) const override
  {
    const int det_id = getDetID(iat);
    Dets[det_id]->mw_evalGrad(extract_DetRef_list(wfc_list, det_id), p_list, iat, grad_now);
  }

  virtual GradType evalGradSource(ParticleSet& P, ParticleSet& src, int iat) override
  {
    GradType G = GradType();
    for (int iz = 0; iz < size(); iz++)
      G += Dets[iz]->evalGradSource(P, src, iat);
    return G;
  }

  virtual GradType evalGradSource(ParticleSet& P,
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

  virtual inline void restore(int iat) override { return Dets[getDetID(iat)]->restore(iat); }

  virtual inline void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override
  {
    Dets[getDetID(iat)]->acceptMove(P, iat, safe_to_delay);

    LogValue = 0.0;
    for (int i = 0; i < Dets.size(); ++i)
      LogValue += Dets[i]->LogValue;
  }

  virtual void mw_accept_rejectMove(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                    const RefVectorWithLeader<ParticleSet>& p_list,
                                    int iat,
                                    const std::vector<bool>& isAccepted,
                                    bool safe_to_delay = false) const override
  {
    constexpr LogValueType czero(0);

    for (int iw = 0; iw < wfc_list.size(); iw++)
      if (isAccepted[iw])
        wfc_list[iw].LogValue = czero;

    for (int i = 0; i < Dets.size(); ++i)
    {
      const auto Det_list(extract_DetRef_list(wfc_list, i));

      if (i == getDetID(iat))
        Dets[i]->mw_accept_rejectMove(Det_list, p_list, iat, isAccepted, safe_to_delay);

      for (int iw = 0; iw < wfc_list.size(); iw++)
        if (isAccepted[iw])
          wfc_list[iw].LogValue += Det_list[iw].LogValue;
    }
  }

  virtual void completeUpdates() override
  {
    for (int i = 0; i < Dets.size(); i++)
      Dets[i]->completeUpdates();
  }

  virtual void mw_completeUpdates(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override
  {
    for (int i = 0; i < Dets.size(); i++)
      Dets[i]->mw_completeUpdates(extract_DetRef_list(wfc_list, i));
  }

  virtual inline PsiValueType ratio(ParticleSet& P, int iat) override { return Dets[getDetID(iat)]->ratio(P, iat); }

  virtual void mw_calcRatio(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                            const RefVectorWithLeader<ParticleSet>& p_list,
                            int iat,
                            std::vector<PsiValueType>& ratios) const override
  {
    const int det_id = getDetID(iat);
    Dets[det_id]->mw_calcRatio(extract_DetRef_list(wfc_list, det_id), p_list, iat, ratios);
  }

  virtual WaveFunctionComponentPtr makeClone(ParticleSet& tqp) const override;

  virtual SPOSetPtr getPhi(int i = 0) { return Dets[i]->getPhi(); }

  virtual void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios) override;

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
  void recompute(MCWalkerConfiguration& W, bool firstTime)
  {
    for (int id = 0; id < Dets.size(); id++)
      Dets[id]->recompute(W, firstTime);
  }

  void reserve(PointerPool<gpu::device_vector<CTS::ValueType>>& pool, int kblocksize = 1)
  {
    for (int id = 0; id < Dets.size(); id++)
      Dets[id]->reserve(pool, kblocksize);
  }

  void addLog(MCWalkerConfiguration& W, std::vector<RealType>& logPsi)
  {
    for (int id = 0; id < Dets.size(); id++)
      Dets[id]->addLog(W, logPsi);
  }

  void ratio(MCWalkerConfiguration& W,
             int iat,
             std::vector<ValueType>& psi_ratios,
             std::vector<GradType>& grad,
             std::vector<ValueType>& lapl)
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
                     int nw)
  {
    Dets[getDetID(iat)]->det_lookahead(W, psi_ratios, grad, lapl, iat, k, kd, nw);
  }
  void calcRatio(MCWalkerConfiguration& W,
                 int iat,
                 std::vector<ValueType>& psi_ratios,
                 std::vector<GradType>& grad,
                 std::vector<ValueType>& lapl)
  {
    Dets[getDetID(iat)]->calcRatio(W, iat, psi_ratios, grad, lapl);
  }

  void addRatio(MCWalkerConfiguration& W,
                int iat,
                int k,
                std::vector<ValueType>& psi_ratios,
                std::vector<GradType>& grad,
                std::vector<ValueType>& lapl)
  {
    Dets[getDetID(iat)]->addRatio(W, iat, k, psi_ratios, grad, lapl);
  }

  void ratio(std::vector<Walker_t*>& walkers,
             std::vector<int>& iatList,
             std::vector<PosType>& rNew,
             std::vector<ValueType>& psi_ratios,
             std::vector<GradType>& grad,
             std::vector<ValueType>& lapl);

  void calcGradient(MCWalkerConfiguration& W, int iat, int k, std::vector<GradType>& grad)
  {
    Dets[getDetID(iat)]->calcGradient(W, iat, k, grad);
  }

  void addGradient(MCWalkerConfiguration& W, int iat, std::vector<GradType>& grad)
  {
    Dets[getDetID(iat)]->addGradient(W, iat, grad);
  }

  void update(MCWalkerConfiguration* W, std::vector<Walker_t*>& walkers, int iat, std::vector<bool>* acc, int k)
  {
    Dets[getDetID(iat)]->update(W, walkers, iat, acc, k);
  }

  void update(const std::vector<Walker_t*>& walkers, const std::vector<int>& iatList);

  void gradLapl(MCWalkerConfiguration& W, GradMatrix_t& grads, ValueMatrix_t& lapl)
  {
    for (int id = 0; id < Dets.size(); id++)
      Dets[id]->gradLapl(W, grads, lapl);
  }

  void NLratios(MCWalkerConfiguration& W,
                std::vector<NLjob>& jobList,
                std::vector<PosType>& quadPoints,
                std::vector<ValueType>& psi_ratios)
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
