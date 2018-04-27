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
  virtual void checkInVariables(opt_variables_type &active) override {}
  virtual void checkOutVariables(const opt_variables_type &active) override {}
  virtual void resetParameters(const opt_variables_type &active) override {}
  virtual void reportStatus(std::ostream& os) override {}
  virtual void resetTargetParticleSet(ParticleSet& P) override {}

  ValueType FakeGradRatio;

  ConstantOrbital() : FakeGradRatio(1.0) {}

  virtual RealType
  evaluateLog(ParticleSet& P,
              ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L) override
  {
    G = 0.0;
    L = 0.0;
    return 0.0;
  }

  virtual void acceptMove(ParticleSet& P, int iat) override {}

  virtual void restore(int iat) override {}

  virtual ValueType ratio(ParticleSet& P, int iat) override
  {
    return 1.0;
  }

  virtual GradType evalGrad(ParticleSet &P, int iat) override
  {
    return GradType(0.0);
  }

  virtual ValueType ratioGrad(ParticleSet &P, int iat, GradType& grad_iat) override
  {
    return FakeGradRatio;
  }

  virtual void registerData(ParticleSet& P, WFBufferType& buf) override {}

  virtual RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false) override {return 0.0;}

  virtual void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override {}

  virtual OrbitalBasePtr makeClone(ParticleSet& tpq) const override {
    return new ConstantOrbital();
  }

};


}
#endif
