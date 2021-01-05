//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: SpinDensity.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SPINDENSITYNEW_H
#define QMCPLUSPLUS_SPINDENSITYNEW_H

#include "SpinDensityInput.h"

#include <vector>
#include <variant>

#include "Configuration.h"
#include "OperatorEstBase.h"
#include "Containers/OhmmsPETE/TinyVector.h"
#include "Utilities/SpeciesSet.h"

namespace qmcplusplus
{

/** Class that collects density per species of particle
 *
 *  commonly used for spin up and down electrons
 *  
 */
class SpinDensityNew : public OperatorEstBase
{
public:
  using POLT    = PtclOnLatticeTraits;
  using Lattice = POLT::ParticleLayout_t;
  using QMCT = QMCTraits;

  //  typedef std::vector<RealType> dens_t;
  //  typedef std::vector<PosType> pts_t;

  //data members
  SpinDensityInput input_;
  SpeciesSet species_;

  // this is a bit of a mess to get from SpeciesSet
  std::vector<int> species_size_;

  //constructor/destructor
  SpinDensityNew(SpinDensityInput& sdi, const SpeciesSet& species);
  ~SpinDensityNew() {}
  SpinDensityNew(const SpinDensityNew& sdn);
  
  //standard interface
  OperatorEstBase* clone();
  void accumulate(RefVector<MCPWalker>& walkers, RefVector<ParticleSet>& psets);

  /** this allows the EstimatorManagerNew to reduce without needing to know the details
   *  of SpinDensityNew's data.
   */
  void collect(const OperatorEstBase&  oeb);

  /** this gets us into the hdf5 file
   *
   *  Just parroting for now don't fully understand.
   *, needs to be unraveled and simplified the hdf5 output is another 
   *  big state big coupling design.
   */
  void registerOperatorEstimator(std::vector<observable_helper*>& h5desc, hid_t gid) const;

private:
  /** data management
   */
  static Data createLocalData(size_t size, DataLocality data_locality);

  //local functions
  void reset();
  void report(const std::string& pad);
};

} // namespace qmcplusplus

#endif /* QMCPLUSPLUS_SPINDENSITYNEW_H */
