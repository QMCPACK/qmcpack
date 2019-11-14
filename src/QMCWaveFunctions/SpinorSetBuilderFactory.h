//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SPINORBASISSETFACTORY_H
#define QMCPLUSPLUS_SPINORBASISSETFACTORY_H

#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"
#include "QMCWaveFunctions/SPOSetBuilder.h"
#include "QMCWaveFunctions/SPOSetBuilderFactory.h"

namespace qmcplusplus
{

class SpinorSetBuilderFactory  : public SPOSetBuilderFactory
{
public:
  typedef std::map<std::string, ParticleSet*> PtclPoolType;

  ///set of basis set builders resolved by type
//  static std::map<std::string, SPOSetBuilder*> spo_builders;

  /// Reset the map and last_builder pointers.  Mostly for unit tests.
//  static void clear();

  SpinorSetBuilderFactory(Communicate* comm, ParticleSet& els, PtclPoolType& psets);

  ~SpinorSetBuilderFactory(){};
  SPOSetBuilder* createSPOSetBuilder(xmlNodePtr rootNode);

//  void loadBasisSetFromXML(xmlNodePtr cur) { last_builder->loadBasisSetFromXML(cur); }
//
//  SPOSet* createSPOSet(xmlNodePtr cur);

//  void build_sposet_collection(xmlNodePtr cur);
//
//private:
  ///store the last builder, use if type not provided
//  static SPOSetBuilder* last_builder;

  ///reference to the target particle
 // ParticleSet& targetPtcl;

  ///reference to the particle pool
//  PtclPoolType& ptclPool;

//  static std::string basisset_tag; 
};
} // namespace qmcplusplus
#endif
