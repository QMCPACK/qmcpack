//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim and Simone Chiesa
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
#ifndef QMCPLUSPLUS_NONLOCAL_ECPOTENTIAL_H
#define QMCPLUSPLUS_NONLOCAL_ECPOTENTIAL_H
#include "QMCHamiltonians/NonLocalECPComponent.h"
#include "QMCHamiltonians/ForceBase.h"

namespace qmcplusplus
{

/** @ingroup hamiltonian
 * \brief Evaluate the semi local potentials
 */
class NonLocalECPotential: public QMCHamiltonianBase, public ForceBase
{
  public:
  ///number of ions
  int NumIons;
  ///index of distance table for the ion-el pair
  int myTableIndex;
  ///the set of local-potentials (one for each ion)
  vector<NonLocalECPComponent*> PP;
  ///unique NonLocalECPComponent to remove
  vector<NonLocalECPComponent*> PPset;
  ///reference to the center ion
  ParticleSet& IonConfig;
  ///target TrialWaveFunction
  TrialWaveFunction& Psi;
  ///true if we should compute forces
  bool ComputeForces;
  ParticleSet::ParticlePos_t PulayTerm;
  ///single particle trace samples
  Array<TraceReal,1>* Ve_sample;
  Array<TraceReal,1>* Vi_sample;
  ParticleSet& Peln;
  ParticleSet& Pion;

  NonLocalECPotential(ParticleSet& ions, ParticleSet& els,
                      TrialWaveFunction& psi, bool computeForces=false);

  ~NonLocalECPotential();

  void resetTargetParticleSet(ParticleSet& P);

  virtual void checkout_particle_arrays(TraceManager& tm);
  virtual void delete_particle_arrays();

  Return_t evaluate(ParticleSet& P);

  Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy);

  Return_t evaluateValueAndDerivatives(ParticleSet& P,
      const opt_variables_type& optvars,
      const vector<RealType>& dlogpsi,
      vector<RealType>& dhpsioverpsi);

  /** Do nothing */
  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "NonLocalECPotential: " << IonConfig.getName();
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  void add(int groupID, NonLocalECPComponent* pp);

  void setRandomGenerator(RandomGenerator_t* rng);

  void addObservables(PropertySetType& plist, BufferType& collectables);

  void setObservables(PropertySetType& plist);

  void setParticlePropertyList(PropertySetType& plist, int offset);

  void registerObservables(vector<observable_helper*>& h5list,
                           hid_t gid) const;
};
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

