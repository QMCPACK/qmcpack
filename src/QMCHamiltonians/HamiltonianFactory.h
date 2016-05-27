/////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file HamiltonianFactory.h
 *@brief Declaration of a HamiltonianFactory
 */
#ifndef QMCPLUSPLUS_HAMILTONIAN_FACTORY_H
#define QMCPLUSPLUS_HAMILTONIAN_FACTORY_H

#include <QMCHamiltonians/QMCHamiltonian.h>
#include <QMCWaveFunctions/WaveFunctionFactory.h>
namespace qmcplusplus
{

/** Factory class to build a many-body wavefunction
 */
class HamiltonianFactory: public MPIObjectBase
{
  public:
  typedef std::map<std::string,ParticleSet*> PtclPoolType;
  typedef std::map<std::string,WaveFunctionFactory*> OrbitalPoolType;

  ///type of the lattice. 0=non-periodic, 1=periodic
  int PBCType;
  ///target ParticleSet
  ParticleSet* targetPtcl;
  ///many-body wavefunction object
  QMCHamiltonian* targetH;
  ///reference to the PtclPoolType
  PtclPoolType&  ptclPool;
  ///reference to the WaveFunctionFactory Pool
  OrbitalPoolType&  psiPool;
  ///input node for a many-body wavefunction
  xmlNodePtr myNode;

  ///name of the TrialWaveFunction
  std::string psiName;

  ///list of the old to new name
  std::map<std::string,std::string> RenamedProperty;

  ///constructor
  HamiltonianFactory(ParticleSet* qp, PtclPoolType& pset, OrbitalPoolType& oset,
                     Communicate* c);

  ///destructor
  ~HamiltonianFactory();

  ///read from xmlNode
  bool put(xmlNodePtr cur);

  ///reset member data
  void reset();

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
  void renameProperty( std::string& a);

  void setCloneSize(int np);

  HamiltonianFactory* clone(ParticleSet* qp, TrialWaveFunction* psi,
                            int ip, const std::string& aname);

  std::vector<HamiltonianFactory*> myClones;

  private:
  /** process xmlNode to populate targetPsi
   */
  bool build(xmlNodePtr cur, bool buildtree);

  void addCoulombPotential(xmlNodePtr cur);
  void addForceHam(xmlNodePtr cur);
  void addPseudoPotential(xmlNodePtr cur);
  void addCorePolPotential(xmlNodePtr cur);
  void addModInsKE(xmlNodePtr cur);
  void addMPCPotential(xmlNodePtr cur, bool physical=false);
  void addVHXCPotential(xmlNodePtr cur);

};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
