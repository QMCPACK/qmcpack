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
#ifndef QMCPLUSPLUS_BSPLINE_JASTROW_BUILDER_H
#define QMCPLUSPLUS_BSPLINE_JASTROW_BUILDER_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"

namespace qmcplusplus
{
//forward declaration
class ParticleSet;

struct BsplineJastrowBuilder: public OrbitalBuilderBase
{
  ParticleSet *sourcePtcl;
  // One-body constructor
  BsplineJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi,
                        ParticleSet& source) :
    OrbitalBuilderBase(target,psi), sourcePtcl(&source)
  {
    ClassName="BsplineJastrowBuilder";
  }
  // Two-body constructor
  BsplineJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi) :
    OrbitalBuilderBase(target,psi), sourcePtcl(NULL)
  {
    ClassName="BsplineJastrowBuilder";
  }

  /** create one-body Jastrow
   * @tparm OBJT one-body jastrow orbital class
   */
  template<typename OBJT, typename DOBJT>
  bool createOneBodyJastrow(xmlNodePtr cur);

  bool put(xmlNodePtr cur);
};

}
#endif
