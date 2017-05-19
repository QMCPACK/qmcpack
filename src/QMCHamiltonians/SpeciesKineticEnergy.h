//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Yubo Yang, paul.young.0414@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Yubo Yang, paul.young.0414@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SPECIESKINETICENERGY_H
#define QMCPLUSPLUS_SPECIESKINETICENERGY_H

#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus
{

/** SpeciesKineticEnergy evaluate the kinetic energy of each species in the target
 * particle set separately instead of sum over every particle in the set such as BareKinetic
 *  <estimator type="specieskinetic" name="skinetic"/>
 * By default skinetic_u, skinetic_d, etc. columns will be added to scalar.dat.
 *  If hdf5="yes", then data will be added to stat.h5 as well.
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
