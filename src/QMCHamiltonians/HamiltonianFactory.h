//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: John R. Gergely,  University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file HamiltonianFactory.h
 *@brief Declaration of a HamiltonianFactory
 */
#ifndef QMCPLUSPLUS_HAMILTONIAN_FACTORY_H
#define QMCPLUSPLUS_HAMILTONIAN_FACTORY_H

#include "QMCHamiltonians/QMCHamiltonian.h"
namespace qmcplusplus
{
/** Factory class to build a many-body wavefunction
 */
class HamiltonianFactory : public MPIObjectBase
{
public:
  using PSetMap     = std::map<std::string, const std::unique_ptr<ParticleSet>>;
  using PsiPoolType = std::map<std::string, const std::unique_ptr<TrialWaveFunction>>;

  ///constructor
  HamiltonianFactory(const std::string& hName,
                     ParticleSet& qp,
                     const PSetMap& pset,
                     const PsiPoolType& oset,
                     Communicate* c);

  ///read from xmlNode
  bool put(xmlNodePtr cur);

  ///get targetH
  std::unique_ptr<QMCHamiltonian> releaseHamiltonian() { return std::move(targetH); }

private:
  /** process xmlNode to populate targetPsi
   */
  bool build(xmlNodePtr cur);

  void addCoulombPotential(xmlNodePtr cur);
  void addForceHam(xmlNodePtr cur);
  void addPseudoPotential(xmlNodePtr cur);
  void addMPCPotential(xmlNodePtr cur, bool physical = false);

  ///type of the lattice. 0=non-periodic, 1=periodic
  int PBCType;
  ///many-body wavefunction object
  std::unique_ptr<QMCHamiltonian> targetH;
  ///target ParticleSet
  ParticleSet& targetPtcl;
  ///reference to the PSetMap
  const PSetMap& ptclPool;
  ///reference to the TrialWaveFunction Pool
  const PsiPoolType& psiPool;

  ///name of the TrialWaveFunction
  std::string psiName;
};
} // namespace qmcplusplus
#endif
