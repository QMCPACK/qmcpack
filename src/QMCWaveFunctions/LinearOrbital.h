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

#include "QMCWaveFunctions/OrbitalBase.h"

namespace qmcplusplus
{

// Only for testing - has constant gradient and zero laplacian
// phi = coeff[0]*x + coeff[1]*y + coeff[2]*z


class LinearOrbital: public OrbitalBase
{
public:
  virtual void checkInVariables(opt_variables_type &active) {}
  virtual void checkOutVariables(const opt_variables_type &active) {}
  virtual void resetParameters(const opt_variables_type &active) {}
  virtual void reportStatus(std::ostream& os) {}
  virtual void resetTargetParticleSet(ParticleSet& P) {}

  TinyVector<ValueType,3> coeff;

  LinearOrbital() {coeff[0]=1.0; coeff[1]= 2.0; coeff[2]=3.0;}

  virtual ValueType
  evaluate(ParticleSet& P,
           ParticleSet::ParticleGradient_t& G,
           ParticleSet::ParticleLaplacian_t& L)
  {
    APP_ABORT("LinearOrbital. evaluate");
    ValueType v = 0.0;
    for (int i = 0; i < P.R.size(); i++) {
      for (int d = 0; d < OHMMS_DIM; d++) {
        v += coeff[d]*P.R[i][d];
      }
      G[i] = coeff;
    }
    L = 0.0;
    return v;
  }

  virtual RealType
  evaluateLog(ParticleSet& P,
              ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
  {
    APP_ABORT("LinearOrbital. evaluateLog");
    G = 0.0;
    L = 0.0;
    return 0.0;
  }

  virtual void acceptMove(ParticleSet& P, int iat) {}

  virtual void restore(int iat) {}

  virtual ValueType ratio(ParticleSet& P, int iat)
  {
    return 1.0;
  }

  virtual GradType evalGrad(ParticleSet &P, int iat)
  {
    return GradType(coeff);
  }

  virtual ValueType ratioGrad(ParticleSet &P, int iat, GradType& grad_iat)
  {
    return 1.0;
  }

  virtual void registerData(ParticleSet& P, WFBufferType& buf) {}

  virtual RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false) {return 0.0;}

  virtual void copyFromBuffer(ParticleSet& P, WFBufferType& buf) {}

};


}
#endif
