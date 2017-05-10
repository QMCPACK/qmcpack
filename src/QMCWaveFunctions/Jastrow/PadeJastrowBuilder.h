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
