//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_SPECIESKINETICENERGY_H
#define QMCPLUSPLUS_SPECIESKINETICENERGY_H

#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus
{

/** SpeciesKineticEnergy evaluate the kinetic energy of each species in the target
 * particle set separately instead of sum over every particle in the set such as BareKinetic
 *  <estimator type="specieskinetic" name="skinetic"/>
 * By default skinetic_u, skinetic_d, etc. columns will be added to scalar.dat, turn off
 *  using scalar="no". If hdf5="yes", then data will be added to stat.h5 as well.
 * The sum of every column that starts with skinetic should be equivalent to the Kinetic column.
 */
class SpeciesKineticEnergy: public QMCHamiltonianBase
{
public:

  SpeciesKineticEnergy(ParticleSet& P);

  bool put(xmlNodePtr cur);         // read input xml node, required
  bool get(std::ostream& os) const; // class description, required

  Return_t evaluate(ParticleSet& P);
  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  { // delegate responsity inline for speed
    return evaluate(P);
  }

  // pure virtual functions require overrider
  void resetTargetParticleSet(ParticleSet& P) { }                         // required
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi); // required

  // allocate multiple columns in scalar.dat
  void addObservables(PropertySetType& plist, BufferType& collectables);
  // fill multiple columns in scalar.dat
  void setObservables(PropertySetType& plist);

  // allow h5 output
  void registerCollectables(std::vector<observable_helper*>& h5desc, hid_t gid) const;

private:
  ParticleSet& tpset; // reference to target particle set
  int num_species;
  std::vector<std::string> species_names;
  std::vector<RealType> species_kinetic,vec_minus_over_2m;
  bool hdf5_out;
  int h5_index; // index of this estimator in the collectables carried by target pset
  //  myIndex: the index of this estimator in the property list in target pset

}; // SpeciesKineticEnergy

} // namespace qmcplusplus
#endif

/***************************************************************************
 * $RCSfile$   $Author: yyang $
 * $Revision: 7049 $   $Date: 2016-08-03 20:47:23 -0500 (Wed, 3 Aug 2017) $
 * $Id: SpeciesKineticEnergy.h 7049 2017-08-03 20:47:23 yyang $
 ***************************************************************************/

