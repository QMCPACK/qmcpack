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

  //  typedef std::vector<RealType> dens_t;
  //  typedef std::vector<PosType> pts_t;

  //data members
  SpinDensityInput input_;
  SpeciesSet species_;

  // this is a bit of a mess to get from SpeciesSet
  std::vector<int> species_size;

  
  //constructor/destructor
  SpinDensityNew(SpinDensityInput&& sdi, const SpeciesSet& species);
  ~SpinDensityNew() {}

  //standard interface
  OperatorEstBase* makeClone(ParticleSet& P, TrialWaveFunction& psi);
  QMCT::FullPrecRealType evaluate(ParticleSet& P);

  //obsolete?
  bool get(std::ostream& os) const { return false; }

  //local functions
  void reset();
  void report(const std::string& pad);
};

} // namespace qmcplusplus

#endif /* QMCPLUSPLUS_SPINDENSITYNEW_H */
