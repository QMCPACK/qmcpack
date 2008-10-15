//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
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

  void ComboOrbital::resetParameters(const opt_variables_type& active) 
  {
    //APP_ABORT("ComboOrbital::resetParameters is incomplete");
    //Constraints->resetParameters(active);
    for(int i=0; i<Psi.size(); i++) Psi[i]->resetParameters(active);
  }

  void ComboOrbital::checkOutVariables(const opt_variables_type& o)
  {
    for(int i=0; i<Psi.size(); i++) Psi[i]->checkOutVariables(o);
  }

  void ComboOrbital::checkInVariables(opt_variables_type& o)
  {
    for(int i=0; i<Psi.size(); i++) Psi[i]->checkInVariables(o);
  }

  void ComboOrbital::reportStatus(ostream& os)
  {
    for(int i=0; i<Psi.size(); i++) Psi[i]->reportStatus(os);
  }

  void ComboOrbital::resetTargetParticleSet(ParticleSet& P) {
    for(int i=0; i<Psi.size(); i++) Psi[i]->resetTargetParticleSet(P);
  }

  ComboOrbital::RealType
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

  //ComboOrbital::ValueType 
  //  ComboOrbital::logRatio(ParticleSet& P, int iat,
  //      ParticleSet::ParticleGradient_t& dG,
  //      ParticleSet::ParticleLaplacian_t& dL) {
  //    ValueType r(0.0);
  //    for(int i=0; i<Psi.size(); i++) 
  //      r += Psi[i]->logRatio(P,iat,dG,dL);
  //    return r;
  //  }

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

  ComboOrbital::RealType 
    ComboOrbital::registerData(ParticleSet& P, BufferType& buf) {
      LogValue=0.0;
      for(int i=0; i<Psi.size(); i++) 
        LogValue += Psi[i]->registerData(P,buf);
      return LogValue;
    }

  ComboOrbital::RealType 
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

  ComboOrbital::RealType 
    ComboOrbital::evaluateLog(ParticleSet& P,BufferType& buf) {
      LogValue=0.0;
      for(int i=0; i<Psi.size(); i++) 
        LogValue += Psi[i]->evaluateLog(P,buf);
      return LogValue;
    }

  OrbitalBase* ComboOrbital::makeClone(ParticleSet& tpq) const
  {
    app_warning() << "  ComboOrbital::makeClone for long-range breakup stuff won't work." << endl;
    ComboOrbital* myclone=new ComboOrbital(*this);
    for(int i=0; i<Psi.size(); ++i)
    {
      myclone->Psi[i]=Psi[i]->makeClone(tpq);
    }
    return myclone;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
