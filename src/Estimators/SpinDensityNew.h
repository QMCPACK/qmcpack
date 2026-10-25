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

#include <limits>
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
   *   its lifetime; an implicit measurement grid co-moves with those changes. A fully
   *   open simulation cell requires an explicit measurement cell in @p sdi. When
   *   @p sdi enables folding, the initial simulation cell must be fully periodic and
   *   commensurate with the explicit measurement cell.
   * @param species Species associated with the measured particle set. It must outlive
   *   this estimator.
   * @param dl Data locality for estimator accumulation.
   *
   * All sdi input results in a grid of fixed size being constructed.
   * If sdi.hasCustomCell() and folding is disabled, the volume, shape, and origin in
   * absolute space that the grid covers are fixed. With folding enabled, the custom
   * geometry is stored in the initial simulation-cell reduced coordinates and co-moves
   * with later simulation-lattice changes. Otherwise a grid is constructed based on the
   * initial mapping of sdi values into the natural coordinates of the simulation
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
  size_t getFullDataSize() const override;
  void accumulateToData(size_t point, QMCT::RealType weight);

  /** @ingroup point resolution functions
   *  All of these functions return whether a point is in the measurement cell the grid covers
   *  Note! point is incremented to the flat index of the grid point corresponding to the position.
   * @{
   */

  /** Get the point when simulation cell has open boundary conditions.
   *  @param position[in]    particle position in native QMCPACK
   *  @param point[in/out]   the index of the grid point the position is binned into.
   */
  bool getCustomMeasurementCellPointForOpenSimulation(const QMCT::PosType& position, size_t& point) const;
  /** Get the point when simulation is periodic and the estimator is folding to a subcell
   *  @param position[in]    particle position in native QMCPACK
   *  @param point[in/out]   the index of the grid point the position is binned into.
   */
  bool getFoldedCustomMeasurementCellPoint(const QMCT::PosType& position, size_t& point) const;

  /** The conservative bounding box for the custom cell in the simulation cell reduced
   *  only used for the non folded custom cell case when simulation cell is periodic
   */
  struct CustomMeasurementCellBounds
  {
    QMCT::PosType lo;
    QMCT::PosType hi;
  };

  /** Get the point when the simulation cell is periodic and the measurement cell is custom.
   *  @param position[in]    particle position in native QMCPACK coordinates.
   *  @param bounds[in]      the conservative bounding box for custom cell
   *  @param point[in/out]   the index of the grid point the position is binned into.
   */
  bool getCustomMeasurementCellPointForPeriodicSimulation(const QMCT::PosType& position,
                                                          const CustomMeasurementCellBounds& bounds,
                                                          size_t& point) const;

  /** This calculates the CustomMeasurementCellBounds
   *  We need to do this for each accumulate call because there is no guarantee the simulation cell hasn't changed.
   */
  CustomMeasurementCellBounds getCustomMeasurementCellBounds() const;

  /**}@*/

  void reset();
  void report(const std::string& pad);

  /// input object
  const SpinDensityInput input_;
  /// The species of particles being accumulated over, determines number of grids
  const SpeciesSet& species_;
  /// this is a bit of a mess to get from SpeciesSet so is pulled out
  const std::vector<int> species_size_;

  static constexpr FullPrecRealType folding_coordinate_tolerance = 64 * std::numeric_limits<QMCT::RealType>::epsilon();
  static constexpr FullPrecRealType folding_commensurability_absolute_tolerance = folding_coordinate_tolerance;
  static constexpr FullPrecRealType folding_commensurability_relative_tolerance =
      64 * std::numeric_limits<FullPrecRealType>::epsilon();

  /** @ingroup Constant Object Parameters
   *  All the following types, members, and functions handle the const
   *  parmeters of a particular SpinDensity estimator. There are serveral different cases the
   *  SpinDensity estimator needs to cover depending on the type of simulation cell, whether
   *  an explicit cell is set by input, and whether it is intended as a unit cell to fold a
   *  simulation supercell down to.
   * @{
   */

  /// the explicit unit cell to fold down to in reduced supercell coordinates
  struct FoldedMeasurementCell
  {
    QMCT::PosType corner_u;
    Lattice lattice_u;
  };

  /// Values derived once from the input and initial simulation lattice.
  struct InstanceParameters
  {
    /// Present only for an explicit <cell>; it is the immutable measurement lattice.
    std::optional<Lattice> custom_measurement_lattice;
    /// Grid dimensions and Cartesian corner are fixed when the estimator is constructed.
    SpinDensityInput::DerivedParameters derived_parameters;
    /// Construction-time Cartesian corner expressed in simulation reduced coordinates.
    QMCT::PosType implicit_corner_u;
    /// Present only with folding; the custom geometry in initial simulation reduced coordinates.
    std::optional<FoldedMeasurementCell> folded_measurement_cell;
  };

  /// Construct the input-dependent parameters and validate the initial measurement geometry.
  static InstanceParameters makeInstanceParameters(const SpinDensityInput& input, const Lattice& simulation_lattice);

  /** Validate the initial custom cell is a possible unit cell
   *  @param ...
   */
  static std::optional<FoldedMeasurementCell> makeFoldedMeasurementCell(
      const SpinDensityInput& input,
      const Lattice& simulation_lattice,
      const std::optional<Lattice>& custom_measurement_lattice,
      const QMCT::PosType& custom_corner);

  /** The simulation lattice is shared so changing cell geometry may is considered possible
   *  If this ever becomes possible think carefully about what the measurement be accumulated now
   *  means.
   */
  const Lattice& simulation_lattice_;
  /// struct containing theb const parameters of this object fixed at construction
  const InstanceParameters instance_parameters_;

  friend class testing::SpinDensityNewTests;
};

} // namespace qmcplusplus

#endif /* QMCPLUSPLUS_SPINDENSITYNEW_H */
