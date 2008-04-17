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
/** @file ComboOrbital.cpp
 * @brief Definitions of ComboOrbital
 */
#include "QMCWaveFunctions/ComboOrbital.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus {

  ComboOrbital::~ComboOrbital() 
  {
    delete_iter(Psi.begin(), Psi.end());
    delete Constraints;
  }

  void ComboOrbital::resetParameters(OptimizableSetType& optVariables) 
  {
    Constraints->resetParameters(optVariables);
  }

  void ComboOrbital::resetTargetParticleSet(ParticleSet& P) {
    for(int i=0; i<Psi.size(); i++) Psi[i]->resetTargetParticleSet(P);
  }

  ComboOrbital::ValueType
    ComboOrbital::evaluateLog(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L) {
      LogValue=0.0;
      for(int i=0; i<Psi.size(); i++) 
        LogValue += Psi[i]->evaluateLog(P,G,L);
      return LogValue;
    }

  ComboOrbital::ValueType 
    ComboOrbital::ratio(ParticleSet& P, int iat,
        ParticleSet::ParticleGradient_t& dG,
        ParticleSet::ParticleLaplacian_t& dL) {
      ValueType r(1.0);
      for(int i=0; i<Psi.size(); i++) 
        r *= Psi[i]->ratio(P,iat,dG,dL);
      return r;
    }

  ComboOrbital::ValueType 
    ComboOrbital::ratio(ParticleSet& P, int iat) {
      ValueType r(1.0);
      for(int i=0; i<Psi.size(); i++) 
        r *= Psi[i]->ratio(P,iat);
      return r;
    }

  ComboOrbital::ValueType 
    ComboOrbital::logRatio(ParticleSet& P, int iat,
        ParticleSet::ParticleGradient_t& dG,
        ParticleSet::ParticleLaplacian_t& dL) {
      ValueType r(0.0);
      for(int i=0; i<Psi.size(); i++) 
        r += Psi[i]->logRatio(P,iat,dG,dL);
      return r;
    }

  void ComboOrbital::acceptMove(ParticleSet& P, int iat) {
    for(int i=0; i<Psi.size(); i++) 
      Psi[i]->acceptMove(P,iat);
  }

  void ComboOrbital::restore(int iat) {
    for(int i=0; i<Psi.size(); i++) 
      Psi[i]->restore(iat);
  }

  void ComboOrbital::update(ParticleSet& P, 
      ParticleSet::ParticleGradient_t& dG, 
      ParticleSet::ParticleLaplacian_t& dL,
      int iat) {
    for(int i=0; i<Psi.size(); i++) 
      Psi[i]->update(P,dG,dL,iat);
  }

  ComboOrbital::ValueType 
    ComboOrbital::registerData(ParticleSet& P, BufferType& buf) {
      LogValue=0.0;
      for(int i=0; i<Psi.size(); i++) 
        LogValue += Psi[i]->registerData(P,buf);
      return LogValue;
    }

  ComboOrbital::ValueType 
    ComboOrbital::updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch) {
      LogValue=0.0;
      for(int i=0; i<Psi.size(); i++) 
        LogValue += Psi[i]->updateBuffer(P,buf,fromscratch);
      return LogValue;
    }

  void 
    ComboOrbital::copyFromBuffer(ParticleSet& P, BufferType& buf){
      for(int i=0; i<Psi.size(); i++) 
        Psi[i]->copyFromBuffer(P,buf);
    }

  ComboOrbital::ValueType 
    ComboOrbital::evaluate(ParticleSet& P,BufferType& buf) {
      LogValue=0.0;
      for(int i=0; i<Psi.size(); i++) 
        LogValue += Psi[i]->evaluate(P,buf);
      return LogValue;
    }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
