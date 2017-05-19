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

#ifndef QMCPLUSPLUS_LATTICEDEVIATION_H
#define QMCPLUSPLUS_LATTICEDEVIATION_H

#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "ParticleBase/ParticleAttribOps.h"
#include <Particle/DistanceTable.h>
#include <Particle/DistanceTableData.h>

namespace qmcplusplus
{

/** lattice deviation estimator
 *
 * Compute deviation of species="tgroup" in target particle set from species="sgroup" in source particle set. The motivation is to observe the deviation of protons from their crystal sites in an all electron-proton simulation of hydrogen i.e. two-component system
   One can also create any reference point to measure the deviation of an up electron, for example
      <particleset name="wf_center">
         <group name="origin" size="1">
            <attrib name="position" datatype="posArray" condition="0">
                     0.00000000        0.00000000        0.00000000
            </attrib>
         </group>
      </particleset>
      <estimator type="latticedeviation" name="latdev" hdf5="yes" per_xyz="yes" 
           source="wf_center" sgroup="origin" target="e" tgroup="u"/> 
   This estimator outputs to both scalar.dat and stat.h5. The scalar.dat entries are averaged over all
particles, whereas the stat.h5 entries are particle-resolved. The two sets of outputs can be compared
as a consistency check for the estimator.
 */
class LatticeDeviationEstimator: public QMCHamiltonianBase
{
public:

  LatticeDeviationEstimator(ParticleSet& P, ParticleSet& sP, const std::string& tgroup, const std::string& sgroup);
  ~LatticeDeviationEstimator() { }

  bool put(xmlNodePtr cur);         // read input xml node, required
  bool get(std::ostream& os) const; // class description, required

  Return_t evaluate(ParticleSet& P); // main function that calculates the observable
  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  { // delegate responsity inline for speed
    return evaluate(P);
  }

  // allow multiple scalars to be registered in scalar.dat
  void addObservables(PropertySetType& plist, BufferType& collectables);
  void setObservables(PropertySetType& plist);
  //void setParticlePropertyList(PropertySetType& plist, int offset); // is this method ever used?

  // make room in hdf5 observable registry
  void registerCollectables(std::vector<observable_helper*>& h5desc, hid_t gid) const;
  //void addObservables(PropertySetType& plist, BufferType& collectables); // also used for multiple scalars

  // pure virtual functions require overrider
  void resetTargetParticleSet(ParticleSet& P);                            // required
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi); // required

private:
  SpeciesSet&       tspecies; // species table of target particle set
  SpeciesSet&       sspecies; // species table of source particle set
  std::string  tgroup,sgroup; // name of species to track
  DistanceTableData* d_table; // distance table between target and source particle sets
  ParticleSet&   tpset,spset; // save references to source and target particle sets
  int num_sites; // number of lattice sites (i.e. number of source particles)
  bool hdf5_out; // use .h5 file for data (follow SkEstimator)
  int  h5_index; // track the starting memory location in P.Collectables
  bool per_xyz;  // track deviation in each of x,y,z directions
  std::vector<RealType> xyz2; // temporary storage for deviation in each of x,y,z directions
  xmlNodePtr input_xml; // original xml
}; // LatticeDeviationEstimator

} // namespace qmcplusplus
#endif
