//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_CONSTANTORBITAL_H
#define QMCPLUSPLUS_CONSTANTORBITAL_H

#include "QMCWaveFunctions/WaveFunctionComponent.h"

namespace qmcplusplus
{
class ConstantOrbital : public WaveFunctionComponent
{
public:
  void checkInVariables(opt_variables_type& active) override {}
  void checkOutVariables(const opt_variables_type& active) override {}
  void resetParameters(const opt_variables_type& active) override {}
  void reportStatus(std::ostream& os) override {}

  PsiValueType FakeGradRatio;

  ConstantOrbital() : WaveFunctionComponent("ConstantOrbital"), FakeGradRatio(1.0) {}

  LogValueType evaluateLog(const ParticleSet& P,
                           ParticleSet::ParticleGradient& G,
                           ParticleSet::ParticleLaplacian& L) override
  {
    G = 0.0;
    L = 0.0;
    return 0.0;
  }

  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override {}

  void restore(int iat) override {}

  PsiValueType ratio(ParticleSet& P, int iat) override { return 1.0; }

  GradType evalGrad(ParticleSet& P, int iat) override { return GradType(0.0); }

  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override { return FakeGradRatio; }

  void registerData(ParticleSet& P, WFBufferType& buf) override {}

  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override { return 0.0; }

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override {}

  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tpq) const override
  {
    return std::make_unique<ConstantOrbital>();
  }
};


} // namespace qmcplusplus
#endif
