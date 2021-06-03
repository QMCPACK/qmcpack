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
  virtual void checkInVariables(opt_variables_type& active) {}
  virtual void checkOutVariables(const opt_variables_type& active) {}
  virtual void resetParameters(const opt_variables_type& active) {}
  virtual void reportStatus(std::ostream& os) {}

  TinyVector<ValueType, 3> coeff;

  LinearOrbital() : WaveFunctionComponent("LinearOrbital")
  {
    coeff[0] = 1.0;
    coeff[1] = 2.0;
    coeff[2] = 3.0;
  }

  virtual LogValueType evaluateLog(const ParticleSet& P,
                                   ParticleSet::ParticleGradient_t& G,
                                   ParticleSet::ParticleLaplacian_t& L)
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
    L        = 0.0;
    LogValue = convertValueToLog(v);
    return LogValue;
  }

  virtual void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) {}

  virtual void restore(int iat) {}

  virtual PsiValueType ratio(ParticleSet& P, int iat) { return 1.0; }

  virtual GradType evalGrad(ParticleSet& P, int iat) { return GradType(coeff); }

  virtual PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) { return 1.0; }

  virtual void registerData(ParticleSet& P, WFBufferType& buf) {}

  virtual LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) { return evaluateLog(P, P.G, P.L); }

  virtual void copyFromBuffer(ParticleSet& P, WFBufferType& buf) {}
};


} // namespace qmcplusplus
#endif
