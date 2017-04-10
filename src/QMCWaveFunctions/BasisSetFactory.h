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
#include "QMCWaveFunctions/BasisSetBase.h"

namespace qmcplusplus
{

  ///writes info about contained sposets to stdout
  void write_basis_builders(const std::string& pad="");

  /**returns a named sposet from the global pool
   *  only use in serial portion of execution
   *  ie during initialization prior to threaded code
   */
  SPOSetBase* get_sposet(const std::string& name);



/** derived class from OrbitalBuilderBase
 */
class BasisSetFactory: public OrbitalBuilderBase
{

public:

  ///set of basis set builders resolved by type
  static std::map<std::string,BasisSetBuilder*> basis_builders;

  /** constructor
   * \param els reference to the electrons
   * \param psi reference to the wavefunction
   * \param ions reference to the ions
   */
  BasisSetFactory(ParticleSet& els, TrialWaveFunction& psi, PtclPoolType& psets);

  ~BasisSetFactory();
  bool put(xmlNodePtr cur);

  BasisSetBuilder* createBasisSet(xmlNodePtr cur, xmlNodePtr rootNode=NULL);

  SPOSetBase* createSPOSet(xmlNodePtr cur);

  void build_sposet_collection(xmlNodePtr cur);

private:

  ///store the last builder, use if type not provided
  static BasisSetBuilder* last_builder;

  ///reference to the particle pool
  PtclPoolType& ptclPool;
};
}
#endif
