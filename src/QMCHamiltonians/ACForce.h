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

#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"

namespace qmcplusplus
{
struct ACForce: public QMCHamiltonianBase
{
  typedef ParticleSet::ParticlePos_t Force_t;
  /** Constructor **/
  ACForce(ParticleSet& source, ParticleSet& target, TrialWaveFunction& psi, QMCHamiltonian& H);
  /** Destructor **/
  ~ACForce(){};
  /** Copy constructor **/
  //ACForce(const ACForce& ac)  {};
 
  /** I/O Routines */
  bool put(xmlNodePtr cur){return true;};
  bool get(std::ostream& os) const {return true;};

  /** Cloning **/
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);
  
  /** Initialization/assignment **/
  void resetTargetParticleSet(ParticleSet& P){};
  void addObservables(PropertySetType& plist, BufferType& collectables);
  void setObservables(PropertySetType& plist);
  void setParticlePropertyList(PropertySetType& plist, int offset);

  /** Evaluate **/
  Return_t evaluate(ParticleSet& P);  

  //** Internal variables **/
  TrialWaveFunction psi;
  ParticleSet ions;
  ParticleSet elns;
  QMCHamiltonian ham;

  //For indexing observables
  IndexType FirstForceIndex;
  IndexType Nions;
  
  //Temporary Nion x 3 dimensional arrays for force storage.
  Force_t hf_force;
  Force_t pulay_force;
  Force_t wf_grad;

  //Class info.
  std::string prefix;

};

}
#endif

