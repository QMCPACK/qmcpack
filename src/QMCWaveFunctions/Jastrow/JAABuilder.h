//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    D.C. Yang, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_ORIGINAL_JASTROW_AA_BUILDER_H
#define QMCPLUSPLUS_ORIGINAL_JASTROW_AA_BUILDER_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"

namespace qmcplusplus
{

//forward declaration
class ParticleSet;

/** Generic two-body Jastrow builder
 *
 * Replacement of JastrowBuilder::createTwoBodySpin and JastrowBuilder::createTwoBodyNoSpin
 */
struct JAABuilder: public OrbitalBuilderBase
{

  JAABuilder(ParticleSet& p, TrialWaveFunction& psi);

  bool put(xmlNodePtr cur);

  template <class FN> TwoBodyJastrowOrbital<FN>* createJAA(xmlNodePtr cur, const std::string& jname);

  bool IgnoreSpin;

  DistanceTableData* d_table;
};

}
#endif
