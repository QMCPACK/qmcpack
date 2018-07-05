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
    
    
#ifndef QMCPLUSPLUS_BASISSETFACTORY_H
#define QMCPLUSPLUS_BASISSETFACTORY_H

#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/SPOSetBuilder.h"

namespace qmcplusplus
{

  ///writes info about contained sposets to stdout
  void write_spo_builders(const std::string& pad="");

  /**returns a named sposet from the global pool
   *  only use in serial portion of execution
   *  ie during initialization prior to threaded code
   */
  SPOSetBase* get_sposet(const std::string& name);



/** derived class from OrbitalBuilderBase
 */
class SPOSetBuilderFactory: public OrbitalBuilderBase
{

public:

  ///set of basis set builders resolved by type
  static std::map<std::string,SPOSetBuilder*> spo_builders;

  /// Reset the map and last_builder pointers.  Mostly for unit tests.
  static void clear();

  /** constructor
   * \param els reference to the electrons
   * \param psi reference to the wavefunction
   * \param ions reference to the ions
   */
  SPOSetBuilderFactory(ParticleSet& els, TrialWaveFunction& psi, PtclPoolType& psets);

  ~SPOSetBuilderFactory();
  bool put(xmlNodePtr cur);

  SPOSetBuilder* createSPOSetBuilder(xmlNodePtr rootNode);

  void loadBasisSetFromXML(xmlNodePtr cur)
  {
    last_builder->loadBasisSetFromXML(cur);
  }

  SPOSetBase* createSPOSet(xmlNodePtr cur);

  void build_sposet_collection(xmlNodePtr cur);

private:

  ///store the last builder, use if type not provided
  static SPOSetBuilder* last_builder;

  ///reference to the particle pool
  PtclPoolType& ptclPool;
};
}
#endif
