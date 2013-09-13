//////////////////////////////////////////////////////////////////
// (c) Copyright 2006- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_LOCALECPPOTENTIAL_H
#define QMCPLUSPLUS_LOCALECPPOTENTIAL_H
#include "Particle/ParticleSet.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimLinearSpline.h"
#include "Numerics/OneDimCubicSpline.h"

namespace qmcplusplus
{

/** @ingroup hamiltonian
 * \brief Evaluate the local potentials (either pseudo or full core) around each ion.
 */

struct LocalECPotential: public QMCHamiltonianBase
{

  typedef OneDimGridBase<RealType> GridType;
  typedef OneDimCubicSpline<RealType> RadialPotentialType;

  ///reference to the ionic configuration
  const ParticleSet& IonConfig;
  ///the number of ioncs
  int NumIons;
  ///distance table index
  int myTableIndex;
  ///temporary energy per particle for pbyp move
  RealType PPtmp;
  ///unique set of local ECP to cleanup
  vector<RadialPotentialType*> PPset;
  ///PP[iat] is the local potential for the iat-th particle
  vector<RadialPotentialType*> PP;
  ///effective charge per ion
  vector<RealType> Zeff;
  ///effective charge per species
  vector<RealType> gZeff;
  ///energy per particle
  Vector<RealType> PPart;
  ///single particle trace samples
  Array<TraceReal,1>* Ve_sample;
  Array<TraceReal,1>* Vi_sample;
  const ParticleSet& Peln;
  const ParticleSet& Pion;

  LocalECPotential(const ParticleSet& ions, ParticleSet& els);

  ~LocalECPotential();

  void resetTargetParticleSet(ParticleSet& P);

  virtual void checkout_particle_arrays(TraceManager& tm);
  virtual void delete_particle_arrays();

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  Return_t spevaluate(ParticleSet& P);
  Return_t evaluate_orig(ParticleSet& P);

  Return_t registerData(ParticleSet& P, BufferType& buffer);
  Return_t updateBuffer(ParticleSet& P, BufferType& buffer);
  void copyFromBuffer(ParticleSet& P, BufferType& buf);
  void copyToBuffer(ParticleSet& P, BufferType& buf);
  Return_t evaluatePbyP(ParticleSet& P, int iat);
  void acceptMove(int iat);

  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "LocalECPotential: " << IonConfig.getName();
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  /** Add a RadialPotentialType of a species
   * @param groupID index of the ion species
   * @param ppot local pseudopotential
   * @param z effective charge of groupID particle
   */
  void add(int groupID, RadialPotentialType* ppot, RealType z);
  Return_t evaluateForPbyP(ParticleSet& P);
};
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

