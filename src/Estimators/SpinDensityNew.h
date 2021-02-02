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
  using QMCT    = QMCTraits;

  //data members
  SpinDensityInput input_;
  SpeciesSet species_;

  /** @ingroup SpinDensity mutable parameters
   *
   *  they should be limited to values that can be changed from input
   *  or are not present explicitly in the SpinDensityInput
   *  @{
   */

  /// Lattice is always local since it is either in the input or a constructor argument.
  Lattice lattice_;
  SpinDensityInput::DerivedParameters derived_parameters_;
  /**}@*/
  
  // this is a bit of a mess to get from SpeciesSet
  std::vector<int> species_size_;

  /** Constructor for SpinDensityInput that contains an explicitly defined cell
   */
  SpinDensityNew(SpinDensityInput&& sdi, const SpeciesSet& species, DataLocality dl = DataLocality::crowd);
  /** Constructor for SpinDensityInput without explicitly defined cell
   *
   *  the crystal lattice should come from the same particle set as the species set.
   *  in case you are tempted to just pass the ParticleSet don't. It clouds the data dependence of
   *  constructing the estimator and creates a strong coupling between the classes.
   *
   *  Ideally when validating input is built up enough there would be only one constructor with
   *  signature
   *
   *      SpinDensityNew(SpinDensityInput&& sdi, SpinDensityInput::DerivedParameters&& dev_par, SpeciesSet species, DataLocality dl);
   */
  SpinDensityNew(SpinDensityInput&& sdi, const Lattice&, const SpeciesSet& species, const DataLocality dl = DataLocality::crowd);
  SpinDensityNew(const SpinDensityNew& sdn);

  void startBlock(int steps) override;
  //standard interface
  OperatorEstBase* clone() override;
  void accumulate(RefVector<MCPWalker>& walkers, RefVector<ParticleSet>& psets) override;

  /** These absolutely must be of this derived type
   */
  void collect(const RefVector<OperatorEstBase>& operator_estimators) override;

  /** this allows the EstimatorManagerNew to reduce without needing to know the details
   *  of SpinDensityNew's data.
   *
   *  can use base class default until crowd level SpinDensity estimators don't have a copy of the density grid.
   */
  //void collect(const OperatorEstBase&  oeb);

  /** this gets us into the hdf5 file
   *
   *  Just parroting for now don't fully understand.
   *, needs to be unraveled and simplified the hdf5 output is another 
   *  big state big coupling design.
   */
  void registerOperatorEstimator(std::vector<observable_helper*>& h5desc, hid_t gid) const override;

private:
  static std::vector<int> getSpeciesSize(SpeciesSet& species);
  /** derived_parameters_ must be valid i.e. initialized with call to input_.calculateDerivedParameters
   */
  size_t getFullDataSize();
  void accumulateToData(size_t point, QMCT::RealType weight);
  void reset();
  void report(const std::string& pad);
};

} // namespace qmcplusplus

#endif /* QMCPLUSPLUS_SPINDENSITYNEW_H */
