//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file WaveFunctionFactory.h
 *@brief Declaration of a WaveFunctionFactory
 */
#ifndef QMCPLUSPLUS_TRIALWAVEFUNCTION_FACTORY_H
#define QMCPLUSPLUS_TRIALWAVEFUNCTION_FACTORY_H

#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"
#include "QMCWaveFunctions/SPOSetBuilderFactory.h"
#include "Message/MPIObjectBase.h"
namespace qmcplusplus
{
/** Factory class to build a many-body wavefunction
 */
class WaveFunctionFactory : public MPIObjectBase
{
public:
  using PtclPoolType = std::map<std::string, ParticleSet*>;

  /** constructor
   * @param psiName name for both the factory and psi
   * @param qp quantum particleset (aka target)
   * @param pset pool of particlesets
   * @param c  communicator
   * @param c  using tasking inside TWF
   */
  WaveFunctionFactory(const std::string& psiName,
                      ParticleSet& qp,
                      const PtclPoolType& pset,
                      Communicate* c,
                      bool tasking = false);

  ///destructor
  ~WaveFunctionFactory();

  ///read from xmlNode
  bool put(xmlNodePtr cur);
  ///get xmlNode
  xmlNodePtr getNode() const { return myNode; }
  ///get targetPsi
  TrialWaveFunction* getTWF() const { return targetPsi.get(); }
  ///get SPOSet
  SPOSet* getSPOSet(const std::string& name) const { return sposet_builder_factory_.getSPOSet(name); }

private:
  /** process xmlNode to populate targetPsi
   */
  bool build(xmlNodePtr cur, bool buildtree = true);

  /** add Fermion wavefunction term */
  bool addFermionTerm(xmlNodePtr cur);

  /** add an OrbitalBuilder and the matching xml node
   * @param b WaveFunctionComponentBuilder*
   * @oaram cur xmlNode for b
   * @return true if successful
   */
  bool addNode(std::unique_ptr<WaveFunctionComponentBuilder> b, xmlNodePtr cur);

  ///many-body wavefunction object
  std::unique_ptr<TrialWaveFunction> targetPsi;
  ///target ParticleSet
  ParticleSet& targetPtcl;
  ///reference to the PtclPoolType
  const PtclPoolType& ptclPool;
  ///input node for a many-body wavefunction
  xmlNodePtr myNode;
  ///builder tree
  UPtrVector<WaveFunctionComponentBuilder> psiBuilder;

  /// factory for all the sposet builders in this WF
  SPOSetBuilderFactory sposet_builder_factory_;
};

} // namespace qmcplusplus
#endif
