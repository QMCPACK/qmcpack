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

#include "QMCWaveFunctions/OrbitalBase.h"

namespace qmcplusplus
{


class ConstantOrbital: public OrbitalBase
{
public:
  virtual void checkInVariables(opt_variables_type &active) {}
  virtual void checkOutVariables(const opt_variables_type &active) {}
  virtual void resetParameters(const opt_variables_type &active) {}
  virtual void reportStatus(std::ostream& os) {}
  virtual void resetTargetParticleSet(ParticleSet& P) {}

  ValueType FakeGradRatio;

  ConstantOrbital() : FakeGradRatio(1.0) {}

  virtual ValueType
  evaluate(ParticleSet& P,
           ParticleSet::ParticleGradient_t& G,
           ParticleSet::ParticleLaplacian_t& L)
  {
    G = 0.0;
    L = 0.0;
    return 1.0;
  }

  virtual RealType
  evaluateLog(ParticleSet& P,
              ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
  {
    G = 0.0;
    L = 0.0;
    return 0.0;
  }

  virtual ValueType ratio(ParticleSet& P, int iat,
                          ParticleSet::ParticleGradient_t& dG,
                          ParticleSet::ParticleLaplacian_t& dL)
  {
    dG = 0.0;
    dL = 0.0;
    return 1.0;
  }

  virtual void acceptMove(ParticleSet& P, int iat) {}

  virtual void restore(int iat) {}

  virtual ValueType ratio(ParticleSet& P, int iat)
  {
    return 1.0;
  }

  virtual GradType evalGrad(ParticleSet &P, int iat)
  {
    return GradType(0.0);
  }

  virtual ValueType ratioGrad(ParticleSet &P, int iat, GradType& grad_iat)
  {
    return FakeGradRatio;
  }

  virtual void registerData(ParticleSet& P, WFBufferType& buf) {}

  virtual RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false) {return 0.0;}

  virtual void copyFromBuffer(ParticleSet& P, WFBufferType& buf) {}

};


}
#endif
