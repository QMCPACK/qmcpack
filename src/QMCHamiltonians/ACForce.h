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
#include "QMCHamiltonians/SpaceWarpTransformation.h"

namespace qmcplusplus
{
class ACForce : public OperatorBase
{
public:
  using Forces = ParticleSet::ParticlePos;

  /** Constructor **/
  ACForce(ParticleSet& source, ParticleSet& target, TrialWaveFunction& psi, QMCHamiltonian& H);

  /** Destructor, "final" triggers a clang warning **/
  ~ACForce() override = default;

  /** I/O Routines */
  bool put(xmlNodePtr cur) final;

  bool get(std::ostream& os) const final;

  /** Cloning **/
  //We don't actually use this makeClone method.  We just put an APP_ABORT here
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;

  //Not derived from base class.  But we need it to properly set the Hamiltonian reference.
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& H);

  /** Initialization/assignment **/
  void resetTargetParticleSet(ParticleSet& P) final;

  void addObservables(PropertySetType& plist, BufferType& collectables) final;

  void setObservables(PropertySetType& plist) final;

  void setParticlePropertyList(PropertySetType& plist, int offset) final;

  /** Since we store a reference to QMCHamiltonian, the baseclass method add2Hamiltonian 
 *  isn't sufficient.  We override it here. **/
  void add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& targetH) final;

  /** Evaluate **/
  Return_t evaluate(ParticleSet& P) final;

private:
  ///Finite difference timestep
  RealType delta_;

  //** Internal variables **/
  //  I'm assuming that psi, ions, elns, and the hamiltonian are bound to this
  //  instantiation.  Making sure no crosstalk happens is the job of whatever clones this.
  ParticleSet& ions_;
  ParticleSet& elns_;
  TrialWaveFunction& psi_;
  QMCHamiltonian& ham_;

  ///For indexing observables
  IndexType first_force_index_;

  ///Temporary Nion x 3 dimensional arrays for force storage.
  Forces hf_force_;
  Forces pulay_force_;
  Forces wf_grad_;
  Forces sw_pulay_;
  Forces sw_grad_;

  bool useSpaceWarp_;

  ///The space warp transformation class.
  SpaceWarpTransformation swt_;
};

} // namespace qmcplusplus
#endif
