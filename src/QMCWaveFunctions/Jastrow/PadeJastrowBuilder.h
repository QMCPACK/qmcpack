//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
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
