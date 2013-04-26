//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Ken Esler
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
#ifndef QMCPLUSPLUS_EEI_JASTROW_BUILDER_H
#define QMCPLUSPLUS_EEI_JASTROW_BUILDER_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"

namespace qmcplusplus
{
//forward declaration
class ParticleSet;

struct eeI_JastrowBuilder: public OrbitalBuilderBase
{
  ParticleSet *sourcePtcl;
  // Two-body constructor
  eeI_JastrowBuilder(ParticleSet& target, TrialWaveFunction& psi,
                     ParticleSet& source) :
    OrbitalBuilderBase(target,psi), sourcePtcl(&source)
  {
    ClassName="eeI_JastrowBuilder";
  }

  bool put(xmlNodePtr cur);
  template<typename J3type>  bool putkids (xmlNodePtr kids, J3type &J3);
};

}
#endif
