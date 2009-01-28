//////////////////////////////////////////////////////////////////
// (c) Copyright 2009-  by Ken Esler
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: esler@.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-

#ifndef QMCPLUSPLUS_PULAY_FOCE_H
#define QMCPLUSPLUS_PULAY_FOCE_H

#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCHamiltonians/ForceBase.h"

namespace qmcplusplus {

  struct PulayForce : public public QMCHamiltonianBase, public ForceBase
  {
    ParticleSet& Ions;
    ParticleSet& Electrons;

    vector<RealType> WarpNorm;

    ParticleSet::ParticlePos_t GradLogPsi, EGradLogPsi;
    
    CoulombPBCABTemp(ParticleSet& ions, ParticleSet& elns);

    void resetTargetParticleSet(ParticleSet& P);
    Return_t evaluate(ParticleSet& P);

    inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy) 
    { return evaluate(P); }

    bool put(xmlNodePtr cur) 
    { return true; }

    bool get(std::ostream& os) const 
    { return true;    }

    inline RealType
    WarpFunction (RealType r)
    { return 1.0/(r*r*r*r); }

    void addObservables(PropertySetType& plist);

    void setObservables(PropertySetType& plist);

    void setParticlePropertyList(PropertySetType& plist, int offset);
  };

