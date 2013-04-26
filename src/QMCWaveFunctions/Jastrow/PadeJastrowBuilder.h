//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_PADE_JASTROW_BUILDER_H
#define QMCPLUSPLUS_PADE_JASTROW_BUILDER_H
#include "QMCWaveFunctions/Jastrow/PadeFunctors.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"

namespace qmcplusplus
{

/** JastrowBuilder using Pade functor
 *
 * To replace PadeConstraints
 */
struct PadeJastrowBuilder: public OrbitalBuilderBase
{

  PadeJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi, PtclPoolType& psets);
  bool put(xmlNodePtr cur);
  PtclPoolType& ptclPool;
};
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2930 $   $Date: 2008-07-31 10:30:42 -0500 (Thu, 31 Jul 2008) $
 * $Id: PadeConstraints.h 2930 2008-07-31 15:30:42Z jnkim $
 ***************************************************************************/
