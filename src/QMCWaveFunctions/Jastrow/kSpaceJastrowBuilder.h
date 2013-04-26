//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Ken Esler
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: kesler@ciw.edu
//   Tel:
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_KSPACE_JASTROW_BUILDER_H
#define QMCPLUSPLUS_KSPACE_JASTROW_BUILDER_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/Jastrow/kSpaceJastrow.h"

namespace qmcplusplus
{
//forward declaration
class ParticleSet;

struct kSpaceJastrowBuilder: public OrbitalBuilderBase
{
  ParticleSet sourcePtcl;
  std::map<string,kSpaceJastrow::SymmetryType> SymmMap;
  // One-body constructor
  kSpaceJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi,
                       ParticleSet& source) :
    OrbitalBuilderBase(target,psi), sourcePtcl(source)
  {
    // nothing for now
    SymmMap["crystal"]   = kSpaceJastrow::CRYSTAL;
    SymmMap["isotropic"] = kSpaceJastrow::ISOTROPIC;
    SymmMap["none"]      = kSpaceJastrow::NOSYMM;
  }

  bool put(xmlNodePtr cur);
};

}
#endif
