//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


/**@file ACForce.h
 *@brief Declaration of ACForce, Assaraf-Caffarel ZVZB style force estimation.
 */
#ifndef QMCPLUSPLUS_ACFORCE_H
#define QMCPLUSPLUS_ACFORCE_H

#include "QMCHamiltonians/OperatorBase.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"

namespace qmcplusplus
{
struct ACForce : public OperatorBase
{
  typedef ParticleSet::ParticlePos_t Force_t;
  /** Constructor **/
  ACForce(ParticleSet& source, ParticleSet& target, TrialWaveFunction& psi, QMCHamiltonian& H);
  /** Destructor **/
  ~ACForce(){};
  /** Copy constructor **/
  //ACForce(const ACForce& ac)  {};

  /** I/O Routines */
  bool put(xmlNodePtr cur) { return true; };
  bool get(std::ostream& os) const { return true; };

  /** Cloning **/
  //We don't actually use this makeClone method.  We just put an APP_ABORT here
  OperatorBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);
  //Not derived from base class.  But we need it to properly set the Hamiltonian reference.
  OperatorBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& H);

  /** Initialization/assignment **/
  void resetTargetParticleSet(ParticleSet& P){};
  void addObservables(PropertySetType& plist, BufferType& collectables);
  void setObservables(PropertySetType& plist);
  void setParticlePropertyList(PropertySetType& plist, int offset);

  /** Since we store a reference to QMCHamiltonian, the baseclass method add2Hamiltonian 
 *  isn't sufficient.  We override it here. **/
  void add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& targetH);
  /** Evaluate **/
  Return_t evaluate(ParticleSet& P);

  //** Internal variables **/
  //  I'm assuming that psi, ions, elns, and the hamiltonian are bound to this
  //  instantiation.  Making sure no crosstalk happens is the job of whatever clones this.
  TrialWaveFunction& psi;
  ParticleSet& ions;
  ParticleSet& elns;
  QMCHamiltonian& ham;

  //For indexing observables
  IndexType FirstForceIndex;
  IndexType Nions;

  //Temporary Nion x 3 dimensional arrays for force storage.
  Force_t hf_force;
  Force_t pulay_force;
  Force_t wf_grad;

  //Class info.
  std::string prefix;
  //We also set the following from the OperatorBase class.
  //std::string myName;
};

} // namespace qmcplusplus
#endif
