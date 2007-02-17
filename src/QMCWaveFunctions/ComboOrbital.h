//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file ComboOrbital.h
 * @brief Declaration of a Composite Orbital
 */
#ifndef QMCPLUSPLUS_GENERIC_COMBO_WITHCONSTRAINTS_H
#define QMCPLUSPLUS_GENERIC_COMBO_WITHCONSTRAINTS_H
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/OrbitalConstraintsBase.h"

namespace qmcplusplus {

  /** A composite Orbital
   */
  struct ComboOrbital: public OrbitalBase 
  {

    ///A list of OrbitalBase* 
    vector<OrbitalBase*> Psi;
    /** Contraints on Psi
     *
     * Constraints reset optimizable variables that may be shared by Psi. 
     */
    OrbitalConstraintsBase* Constraints;

    ComboOrbital(OrbitalConstraintsBase* control):
      Constraints(control) {
        Optimizable=true;
      }

    ~ComboOrbital();

    void setContraints(OrbitalConstraintsBase* control) {
      Constraints=control;
    }

    /** reset parameters
     */
    void resetParameters(OptimizableSetType& optVariables);

    void resetTargetParticleSet(ParticleSet& P);

    ValueType
      evaluate(ParticleSet& P, 
          ParticleSet::ParticleGradient_t& G, 
          ParticleSet::ParticleLaplacian_t& L) {
        return std::exp(evaluateLog(P,G,L));
      }

    ValueType
      evaluateLog(ParticleSet& P, 
          ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L);

    ValueType 
      ratio(ParticleSet& P, int iat,
          ParticleSet::ParticleGradient_t& dG,
          ParticleSet::ParticleLaplacian_t& dL);

    ValueType 
      ratio(ParticleSet& P, int iat);

    ValueType 
      logRatio(ParticleSet& P, int iat,
          ParticleSet::ParticleGradient_t& dG,
          ParticleSet::ParticleLaplacian_t& dL);

    void acceptMove(ParticleSet& P, int iat);

    void restore(int iat);

    void update(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& dG, 
        ParticleSet::ParticleLaplacian_t& dL,
        int iat);

    ValueType 
      registerData(ParticleSet& P, BufferType& buf);

    ValueType 
      updateBuffer(ParticleSet& P, BufferType& buf);

    void 
      copyFromBuffer(ParticleSet& P, BufferType& buf);

    ValueType 
      evaluate(ParticleSet& P,BufferType& buf);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
