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

#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"
#include "QMCWaveFunctions/SPOSetBuilder.h"
#include "type_traits/template_types.hpp"

namespace qmcplusplus
{
class SPOSetBuilderFactory : public MPIObjectBase
{
public:
  using PtclPoolType = std::map<std::string, ParticleSet*>;

  /** constructor
   * \param comm communicator
   * \param els reference to the electrons
   * \param ions reference to the ions
   */
  SPOSetBuilderFactory(Communicate* comm, ParticleSet& els, const PtclPoolType& psets);

  ~SPOSetBuilderFactory();

  SPOSetBuilder& createSPOSetBuilder(xmlNodePtr rootNode);

  /** returns a named sposet from the pool
   *  only use in serial portion of execution
   *  ie during initialization prior to threaded code
   */
  SPOSet* getSPOSet(const std::string& name) const;

  SPOSetBuilder& getLastBuilder();

  void buildSPOSetCollection(xmlNodePtr cur);

  bool empty() const { return sposet_builders_.size() == 0; }

private:
  ///writes info about contained sposets to stdout
  void write_sposet_builders_(const std::string& pad = "") const;

  ///set of basis set builders resolved by type
  UPtrVector<SPOSetBuilder> sposet_builders_;

  ///reference to the target particle
  ParticleSet& targetPtcl;

  ///reference to the particle pool
  const PtclPoolType& ptclPool;

  static std::string basisset_tag;
};
} // namespace qmcplusplus
#endif
