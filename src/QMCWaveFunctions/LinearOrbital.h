//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_LINEARORBITAL_H
#define QMCPLUSPLUS_LINEARORBITAL_H

#include "QMCWaveFunctions/WaveFunctionComponent.h"

namespace qmcplusplus
{
// Only for testing - has constant gradient and zero laplacian
// phi = coeff[0]*x + coeff[1]*y + coeff[2]*z


class LinearOrbital : public WaveFunctionComponent
{
public:
  TinyVector<ValueType, 3> coeff;

  LinearOrbital()
  {
    coeff[0] = 1.0;
    coeff[1] = 2.0;
    coeff[2] = 3.0;
  }

  std::string getClassName() const override { return "LinearOrbital"; }

  LogValueType evaluateLog(const ParticleSet& P,
                           ParticleSet::ParticleGradient& G,
                           ParticleSet::ParticleLaplacian& L) override
  {
    ValueType v = 0.0;
    for (int i = 0; i < P.R.size(); i++)
    {
      for (int d = 0; d < OHMMS_DIM; d++)
      {
        v += coeff[d] * P.R[i][d];
      }
      G[i] = coeff;
    }
    L          = 0.0;
    log_value_ = convertValueToLog(v);
    return log_value_;
  }

  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override {}

  void restore(int iat) override {}

  PsiValueType ratio(ParticleSet& P, int iat) override { return 1.0; }

  GradType evalGrad(ParticleSet& P, int iat) override { return GradType(coeff); }

  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override { return 1.0; }

  void registerData(ParticleSet& P, WFBufferType& buf) override {}

  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override
  {
    return evaluateLog(P, P.G, P.L);
  }

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override {}

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           Vector<ValueType>& dlogpsi,
                           Vector<ValueType>& dhpsioverpsi) override
  {}
};


} // namespace qmcplusplus
#endif
