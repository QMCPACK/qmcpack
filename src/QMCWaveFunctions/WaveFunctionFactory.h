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
  using PSetMap = std::map<std::string, std::unique_ptr<ParticleSet>>;

  /** constructor
   * @param psiName name for both the factory and psi
   * @param qp quantum particleset (aka target)
   * @param pset pool of particlesets
   * @param c  communicator
   * @param c  using tasking inside TWF
   */
  WaveFunctionFactory(ParticleSet& qp, const PSetMap& pset, Communicate* c);

  ///destructor
  ~WaveFunctionFactory();

  ///read from xmlNode
  std::unique_ptr<TrialWaveFunction> buildTWF(xmlNodePtr cur);

  /// create an empty TrialWaveFunction for testing use.
  std::unique_ptr<TrialWaveFunction> static buildEmptyTWFForTesting(const std::string_view name)
  {
    return std::make_unique<TrialWaveFunction>(name);
  }

private:
  /** add Fermion wavefunction term */
  bool addFermionTerm(TrialWaveFunction& psi, SPOSetBuilderFactory& spo_factory, xmlNodePtr cur);

  ///many-body wavefunction object
  ///target ParticleSet
  ParticleSet& targetPtcl;
  ///reference to the PSetMap
  const PSetMap& ptclPool;
};

} // namespace qmcplusplus
#endif
