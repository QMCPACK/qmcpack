//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_ORIGINAL_JASTROW_AB_BUILDER_H
#define QMCPLUSPLUS_ORIGINAL_JASTROW_AB_BUILDER_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"

namespace qmcplusplus
{

//forward declaration
class ParticleSet;

/** Generic two-body Jastrow builder
 *
 * Replacement of JastrowBuilder::createTwoBodySpin and JastrowBuilder::createTwoBodyNoSpin
 */
struct JABBuilder: public OrbitalBuilderBase
{

  JABBuilder(ParticleSet& p, TrialWaveFunction& psi,
             PtclPoolType& psets):OrbitalBuilderBase(p,psi), ptclPool(psets) {}

  bool put(xmlNodePtr cur);

  template<class FN> bool createJAB(xmlNodePtr cur, const std::string& jname);


  PtclPoolType& ptclPool;
};

}
#endif
