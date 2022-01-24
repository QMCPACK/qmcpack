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
#include "QMCWaveFunctions/WaveFunctionFactory.h"
namespace qmcplusplus
{
/** Factory class to build a many-body wavefunction
 */
class HamiltonianFactory : public MPIObjectBase
{
public:
  using PtclPoolType = std::map<std::string, ParticleSet*>;
  using PsiPoolType  = std::map<std::string, WaveFunctionFactory*>;

  ///constructor
  HamiltonianFactory(const std::string& hName, ParticleSet& qp, PtclPoolType& pset, PsiPoolType& oset, Communicate* c);

  ///read from xmlNode
  bool put(xmlNodePtr cur);

  /** add a property whose name will be renamed by b
   * @param a target property whose name should be replaced by b
   * @param b new property name
   */
  void renameProperty(const std::string& a, const std::string& b);

  /** renamd a property
   * @param a current name
   *
   * If a is found among the RenamedProperty, a is replaced,
   */
  void renameProperty(std::string& a);

  ///get targetH
  QMCHamiltonian* getH() const { return targetH.get(); }

private:
  /** process xmlNode to populate targetPsi
   */
  bool build(xmlNodePtr cur, bool buildtree);

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
  ///reference to the PtclPoolType
  PtclPoolType& ptclPool;
  ///reference to the WaveFunctionFactory Pool
  PsiPoolType& psiPool;
  ///input node for a many-body wavefunction
  xmlNodePtr myNode;

  ///name of the TrialWaveFunction
  std::string psiName;

  ///list of the old to new name
  std::map<std::string, std::string> RenamedProperty;
};
} // namespace qmcplusplus
#endif
