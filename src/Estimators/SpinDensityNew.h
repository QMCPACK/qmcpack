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

#include <optional>
#include <vector>

#include "Configuration.h"
#include "OperatorEstBase.h"
#include "Containers/OhmmsPETE/TinyVector.h"

namespace qmcplusplus
{
class SpeciesSet;
namespace testing
{
class SpinDensityNewTests;
}
/** Class that collects density per species of particle
 *
 *  commonly used for spin up and down electrons
 *
 */
class SpinDensityNew : public OperatorEstBase
{
public:
  using QMCT             = QMCTraits;
  using FullPrecRealType = QMCT::FullPrecRealType;

  /** Construct a spin-density estimator.
   *
   * @param sdi Validated estimator input, moved to this
   *   estimator.
   * @param simulation_lattice The simulation lattice used to locate particles and
   *   enumerate periodic images. It must outlive this estimator and may change during
   *   its lifetime; an implicit measurement grid co-moves with those changes.
   * @param species Species associated with the measured particle set. It must outlive
   *   this estimator.
   * @param dl Data locality for estimator accumulation.
   *
   * All sdi input results in a grid of fixed size being constructed.
   * if sdi.hasCustomCell() then then the volume,shape, and origin in absolute space that
   * grid covers is also fixed. otherwise a grid is constructed based on the initial mapping of
   * sdi values into the the natural coordinates of the simulation
   * lattice.
   * If the simulation_lattice does not change you can continue to
   * interpret the grid in terms of the absolute Bohr units and
   * the cartesian positions of the sdi input. If the simulation_lattice
   * changes over the course of accumulation the density grid only
   * make sense wrt the simulation lattice coordinates.
   */
  SpinDensityNew(SpinDensityInput&& sdi,
                 const Lattice& simulation_lattice,
                 const SpeciesSet& species,
                 const DataLocality dl = DataLocality::crowd);

  /** Constructor used when spawing crowd clones
   *  needs to be public so std::make_unique can call it.
   *  Do not use directly unless you've really thought it through.
   */
  SpinDensityNew(const SpinDensityNew& sdn, DataLocality dl);

  /** This allows us to allocate the necessary data for the DataLocality::queue
   */
  void startBlock(int steps) override;

  /** standard interface
   */
  std::unique_ptr<OperatorEstBase> spawnCrowdClone() const override;

  /** accumulate 1 or more walkers of SpinDensity samples
   *
   *  Accumulation is into a grid that may or may not cover the entire
   *  space particles can sample.
   *
   * Currently each particle is accumulated at most once. i.e. if your
   * grid is defined such that covers the "same" location in the
   * periodic simulation cell more than once you will not acumulate
   * the particle twice.
   * For instance if you define a grid that covers the simulation cell twice
   * your density will infact integrate to np * total_weight not np *
   * total_weight * 2.
   */
  void accumulate(const RefVector<MCPWalker>& walkers,
                  const RefVector<ParticleSet>& psets,
                  const RefVector<TrialWaveFunction>& wfns,
                  const RefVector<QMCHamiltonian>& hams,
                  RandomBase<FullPrecRealType>& rng) override;

  /** this allows the EstimatorManagerNew to reduce without needing to know the details
   *  of SpinDensityNew's data.
   *
   *  can use base class default until crowd level SpinDensity
   *  estimators don't have a copy of the density grid.
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
  void registerOperatorEstimator(hdf_archive& file) override;

private:
  SpinDensityNew(const SpinDensityNew& sdn) = default;

  static std::vector<int> getSpeciesSize(const SpeciesSet& species);
  /** derived_parameters_ must be valid i.e. initialized with call to input_.calculateDerivedParameters
   */
  size_t getFullDataSize() const override;
  void accumulateToData(size_t point, QMCT::RealType weight);
  /// point must initially be the species offset; on success it is the corresponding grid point.
  struct PeriodicFiniteCellBounds
  {
    QMCT::PosType lo;
    QMCT::PosType hi;
  };

  bool getFiniteCellPoint(const QMCT::PosType& position, size_t& point) const;
  PeriodicFiniteCellBounds getPeriodicFiniteCellBounds() const;
  bool getPeriodicFiniteCellPoint(const QMCT::PosType& position,
                                  const PeriodicFiniteCellBounds& bounds,
                                  size_t& point) const;
  void reset();
  void report(const std::string& pad);

  //data members
  const SpinDensityInput input_;
  const SpeciesSet& species_;
  // this is a bit of a mess to get from SpeciesSet
  const std::vector<int> species_size_;

  /** @ingroup SpinDensity mutable parameters
   *
   *  they should be limited to values that can be changed from input
   *  or are not present explicitly in the SpinDensityInput
   *  @{
   */

  /// The simulation lattice is shared so changing cell geometry is observed during accumulation.
  const Lattice& simulation_lattice_;
  /// Present only for an explicit <cell>; it is the immutable measurement lattice.
  const std::optional<Lattice> custom_measurement_lattice_;
  /// Grid dimensions and Cartesian corner are fixed when the estimator is constructed.
  const SpinDensityInput::DerivedParameters derived_parameters_;
  /// Construction-time Cartesian corner expressed in simulation reduced coordinates.
  const QMCT::PosType implicit_corner_u_;
  /**}@*/

  friend class testing::SpinDensityNewTests;
};

} // namespace qmcplusplus

#endif /* QMCPLUSPLUS_SPINDENSITYNEW_H */
