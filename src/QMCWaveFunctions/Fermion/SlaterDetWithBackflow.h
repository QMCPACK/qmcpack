//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SLATERDETERMINANT_WITHBACKFLOW_H
#define QMCPLUSPLUS_SLATERDETERMINANT_WITHBACKFLOW_H
#include "QMCWaveFunctions/Fermion/DiracDeterminantWithBackflow.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include <cmath>

namespace qmcplusplus
{
class SlaterDetWithBackflow : public WaveFunctionComponent
{
public:
  using Determinant_t = DiracDeterminantWithBackflow;
  /**  constructor
   * @param targetPtcl target Particleset
   * @param rn release node
   */
  SlaterDetWithBackflow(ParticleSet& targetPtcl,
                        std::vector<std::unique_ptr<Determinant_t>> dets,
                        std::unique_ptr<BackflowTransformation> BF);

  ///destructor
  ~SlaterDetWithBackflow() override;

  void checkInVariables(opt_variables_type& active) override
  {
    //if(Optimizable) {
    if (BFTrans->isOptimizable())
    {
      BFTrans->checkInVariables(active);
      for (int i = 0; i < Dets.size(); i++)
        Dets[i]->checkInVariables(active);
    }
  }

  void checkOutVariables(const opt_variables_type& active) override
  {
    //if(Optimizable) {
    if (BFTrans->isOptimizable())
    {
      BFTrans->checkOutVariables(active);
      for (int i = 0; i < Dets.size(); i++)
        Dets[i]->checkOutVariables(active);
    }
  }

  ///reset all the Dirac determinants, Optimizable is true
  void resetParameters(const opt_variables_type& active) override
  {
    //if(Optimizable) {
    if (BFTrans->isOptimizable())
    {
      BFTrans->resetParameters(active);
      for (int i = 0; i < Dets.size(); i++)
        Dets[i]->resetParameters(active);
    }
  }

  void reportStatus(std::ostream& os) override {}
  LogValueType evaluateLog(const ParticleSet& P,
                           ParticleSet::ParticleGradient& G,
                           ParticleSet::ParticleLaplacian& L) override;

  void registerData(ParticleSet& P, WFBufferType& buf) override;
  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;
  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;

  inline PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override
  {
    BFTrans->evaluatePbyPWithGrad(P, iat);
    //BFTrans->evaluate(P);
    PsiValueType psi = 1.0;
    for (int i = 0; i < Dets.size(); ++i)
      psi *= Dets[i]->ratioGrad(P, iat, grad_iat);
    return psi;
  }

  GradType evalGrad(ParticleSet& P, int iat) override
  {
    QMCTraits::GradType g;
    for (int i = 0; i < Dets.size(); ++i)
      g += Dets[i]->evalGrad(P, iat);
    return g;
  }

  GradType evalGradSource(ParticleSet& P, ParticleSet& src, int iat) override
  {
    APP_ABORT("Need to implement SlaterDetWithBackflow::evalGradSource() \n");
    return ValueType();
  }

  GradType evalGradSource(ParticleSet& P,
                          ParticleSet& src,
                          int iat,
                          TinyVector<ParticleSet::ParticleGradient, OHMMS_DIM>& grad_grad,
                          TinyVector<ParticleSet::ParticleLaplacian, OHMMS_DIM>& lapl_grad) override
  {
    APP_ABORT("Need to implement SlaterDetWithBackflow::evalGradSource() \n");
    return GradType();
  }

  inline void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override
  {
    BFTrans->acceptMove(P, iat);
    for (int i = 0; i < Dets.size(); i++)
      Dets[i]->acceptMove(P, iat);
  }

  inline void restore(int iat) override
  {
    BFTrans->restore(iat);
    for (int i = 0; i < Dets.size(); i++)
      Dets[i]->restore(iat);
  }


  inline PsiValueType ratio(ParticleSet& P, int iat) override
  {
    BFTrans->evaluatePbyP(P, iat);
    //BFTrans->evaluate(P);
    PsiValueType ratio = 1.0;
    for (int i = 0; i < Dets.size(); ++i)
      ratio *= Dets[i]->ratio(P, iat);
    return ratio;
  }

  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tqp) const override;

  SPOSetPtr getPhi(int i = 0) const { return Dets[i]->getPhi(); }

  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios) override;

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           std::vector<ValueType>& dlogpsi,
                           std::vector<ValueType>& dhpsioverpsi) override;

  void testDerivGL(ParticleSet& P);

private:
  ///container for the DiracDeterminants
  const std::vector<std::unique_ptr<Determinant_t>> Dets;
  /// backflow transformation
  const std::unique_ptr<BackflowTransformation> BFTrans;
};
} // namespace qmcplusplus
#endif
