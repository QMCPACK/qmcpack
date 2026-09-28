//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: SpinDensity.cpp
//////////////////////////////////////////////////////////////////////////////////////


#include "SpinDensityNew.h"

#include "hdf5.h"

#include "Message/UniformCommunicateError.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <SpeciesSet.h>

namespace qmcplusplus
{
SpinDensityNew::SpinDensityNew(SpinDensityInput&& input,
                               const Lattice& lattice,
                               const SpeciesSet& species,
                               const DataLocality dl)
    : OperatorEstBase(dl, input.get_name(), std::string{SpinDensityInput::type_tag}),
      input_(std::move(input)),
      species_(species),
      species_size_(getSpeciesSize(species)),
      simulation_lattice_(lattice),
      instance_parameters_(makeInstanceParameters(input_, simulation_lattice_))
{
  data_locality_ = dl;
  data_.resize(getFullDataSize());
  if (input_.get_write_report())
    report("  ");
}

SpinDensityNew::SpinDensityNew(const SpinDensityNew& sdn, DataLocality dl) : SpinDensityNew(sdn)
{
  data_locality_ = dl;
}

SpinDensityNew::InstanceParameters SpinDensityNew::makeInstanceParameters(const SpinDensityInput& input,
                                                                          const Lattice& simulation_lattice)
{
  std::optional<Lattice> custom_measurement_lattice;
  if (input.hasCustomCell())
    custom_measurement_lattice = input.get_cell();

  if (!custom_measurement_lattice && simulation_lattice.SuperCellEnum == SUPERCELL_OPEN)
    throw UniformCommunicateError("SpinDensity input: an explicit cell is required for a fully open simulation cell");

  const Lattice& measurement_lattice = custom_measurement_lattice ? *custom_measurement_lattice : simulation_lattice;
  SpinDensityInput::DerivedParameters derived_parameters = input.calculateDerivedParameters(measurement_lattice);
  QMCT::PosType implicit_corner_u                        = simulation_lattice.toUnit(derived_parameters.corner);
  std::optional<FoldedMeasurementCell> folded_measurement_cell =
      makeFoldedMeasurementCell(input, simulation_lattice, custom_measurement_lattice, derived_parameters.corner);

  return {std::move(custom_measurement_lattice), std::move(derived_parameters), implicit_corner_u,
          std::move(folded_measurement_cell)};
}

std::vector<int> SpinDensityNew::getSpeciesSize(const SpeciesSet& species)
{
  std::vector<int> species_size;
  int index = species.findAttribute("membersize");
  if (index < 0)
    throw std::runtime_error("SpinDensity(P) Species set does not have the required attribute 'membersize'");
  for (int s = 0; s < species.size(); ++s)
    species_size.push_back(species(index, s));
  return species_size;
}

size_t SpinDensityNew::getFullDataSize() const
{
  return species_.size() * instance_parameters_.derived_parameters.npoints;
}

std::optional<SpinDensityNew::FoldedMeasurementCell> SpinDensityNew::makeFoldedMeasurementCell(
    const SpinDensityInput& input,
    const Lattice& simulation_lattice,
    const std::optional<Lattice>& custom_measurement_lattice,
    const QMCT::PosType& custom_corner)
{
  if (!input.hasFolding())
    return std::nullopt;

  assert(custom_measurement_lattice);
  if (simulation_lattice.SuperCellEnum != SUPERCELL_BULK)
    throw UniformCommunicateError("SpinDensity input: folding requires a fully periodic simulation cell");

  Tensor<FullPrecRealType, QMCT::DIM> rounded_transform;
  for (int i = 0; i < QMCT::DIM; ++i)
  {
    const QMCT::PosType custom_u = custom_measurement_lattice->toUnit(simulation_lattice.Rv[i]);
    for (int j = 0; j < QMCT::DIM; ++j)
    {
      const FullPrecRealType nearest_integer = std::round(custom_u[j]);
      const FullPrecRealType tolerance       = folding_commensurability_absolute_tolerance +
          folding_commensurability_relative_tolerance *
              std::max(std::abs(static_cast<FullPrecRealType>(custom_u[j])), std::abs(nearest_integer));
      if (!std::isfinite(custom_u[j]) || std::abs(custom_u[j] - nearest_integer) > tolerance)
        throw UniformCommunicateError(
            "SpinDensity input: folding requires simulation and measurement cells to be commensurate");
      rounded_transform(i, j) = nearest_integer;
    }
  }

  const FullPrecRealType transform_determinant = det(rounded_transform);
  if (!std::isfinite(transform_determinant) || transform_determinant == 0.0)
    throw UniformCommunicateError("SpinDensity input: folding requires a nonsingular cell transformation");

  Lattice reduced_custom_lattice;
  for (int i = 0; i < QMCT::DIM; ++i)
  {
    const QMCT::PosType axis_u = simulation_lattice.toUnit(custom_measurement_lattice->Rv[i]);
    for (int j = 0; j < QMCT::DIM; ++j)
      reduced_custom_lattice.R(i, j) = axis_u[j];
  }
  reduced_custom_lattice.reset();
  return FoldedMeasurementCell{simulation_lattice.toUnit(custom_corner), std::move(reduced_custom_lattice)};
}

std::unique_ptr<OperatorEstBase> SpinDensityNew::spawnCrowdClone() const
{
  std::size_t data_size    = data_.size();
  auto spawn_data_locality = data_locality_;
  if (data_locality_ == DataLocality::rank)
  {
    spawn_data_locality = DataLocality::queue;
    // at construction we don't know what the data requirement is going to be
    // since its steps per block  dependent. so start with 10 steps worth.
    int num_particles = std::accumulate(species_size_.begin(), species_size_.end(), 0);
    data_size         = num_particles * 20;
  }
  UPtr<SpinDensityNew> spawn(std::make_unique<SpinDensityNew>(*this, spawn_data_locality));
  spawn->get_data().resize(data_size);
  return spawn;
}

void SpinDensityNew::startBlock(int steps)
{
  if (data_locality_ == DataLocality::queue)
  {
    int num_particles = std::accumulate(species_size_.begin(), species_size_.end(), 0);
    size_t data_size  = num_particles * steps * 2;
    data_.reserve(data_size);
    data_.resize(0);
  }
}

void SpinDensityNew::accumulate(const RefVector<MCPWalker>& walkers,
                                const RefVector<ParticleSet>& psets,
                                const RefVector<TrialWaveFunction>& wfns,
                                const RefVector<QMCHamiltonian>& hams,
                                RandomBase<FullPrecRealType>& rng)
{
  const auto& dp_ = instance_parameters_.derived_parameters;
  std::optional<CustomMeasurementCellBounds> custom_measurement_cell_bounds;
  if (input_.hasCustomCell() && !input_.hasFolding() && simulation_lattice_.SuperCellEnum != SUPERCELL_OPEN)
    custom_measurement_cell_bounds = getCustomMeasurementCellBounds();

  for (int iw = 0; iw < walkers.size(); ++iw)
  {
    MCPWalker& walker     = walkers[iw];
    ParticleSet& pset     = psets[iw];
    QMCT::RealType weight = walker.Weight;
    assert(weight >= 0);
    // for testing
    walkers_weight_ += weight;
    int p         = 0;
    size_t offset = 0;
    // important notice the offset increment
    for (int s = 0; s < species_.size(); ++s, offset += dp_.npoints)
      for (int ps = 0; ps < species_size_[s]; ++ps, ++p)
      {
        size_t point = offset;
        // This is the simple path, cell is implicit
        if (!input_.hasCustomCell())
        {
          const QMCT::PosType u = simulation_lattice_.toUnit(pset.R[p]) - instance_parameters_.implicit_corner_u;
          for (int d = 0; d < QMCT::DIM; ++d)
            point += dp_.gdims[d] * static_cast<int>(dp_.grid[d] * (u[d] - std::floor(u[d])));
          accumulateToData(point, weight);
        }
        else if (input_.hasFolding())
        {
          if (getFoldedCustomMeasurementCellPoint(pset.R[p], point))
            accumulateToData(point, weight);
        }
        else if (simulation_lattice_.SuperCellEnum == SUPERCELL_OPEN)
        {
          if (getCustomMeasurementCellPointForOpenSimulation(pset.R[p], point))
            accumulateToData(point, weight);
        }
        else if (getCustomMeasurementCellPointForPeriodicSimulation(pset.R[p], *custom_measurement_cell_bounds, point))
          accumulateToData(point, weight);
      }
  }
}

bool SpinDensityNew::getCustomMeasurementCellPointForOpenSimulation(const QMCT::PosType& position, size_t& point) const
{
  const InstanceParameters& parameters = instance_parameters_;
  const QMCT::PosType u =
      parameters.custom_measurement_lattice->toUnit(position - parameters.derived_parameters.corner);
  size_t candidate_point = point;
  for (int d = 0; d < QMCT::DIM; ++d)
  {
    if (u[d] < 0.0 || u[d] >= 1.0)
      return false;
    candidate_point += parameters.derived_parameters.gdims[d] *
        std::min(static_cast<int>(parameters.derived_parameters.grid[d] * u[d]),
                 parameters.derived_parameters.grid[d] - 1);
  }
  point = candidate_point;
  return true;
}

bool SpinDensityNew::getFoldedCustomMeasurementCellPoint(const QMCT::PosType& position, size_t& point) const
{
  assert(instance_parameters_.folded_measurement_cell);
  const FoldedMeasurementCell& folded_measurement_cell = *instance_parameters_.folded_measurement_cell;
  const QMCT::PosType custom_u =
      folded_measurement_cell.lattice_u.toUnit(simulation_lattice_.toUnit(position) - folded_measurement_cell.corner_u);
  size_t candidate_point = point;
  for (int d = 0; d < QMCT::DIM; ++d)
  {
    FullPrecRealType folded_u = custom_u[d] - std::floor(custom_u[d]);
    if (folded_u >= 1.0 - folding_coordinate_tolerance)
      folded_u = 0.0;

    // Rounding can make a coordinate just below one index as the upper boundary.
    const int bin = std::clamp(static_cast<int>(instance_parameters_.derived_parameters.grid[d] * folded_u), 0,
                               instance_parameters_.derived_parameters.grid[d] - 1);
    candidate_point += instance_parameters_.derived_parameters.gdims[d] * bin;
  }
  point = candidate_point;
  return true;
}

SpinDensityNew::CustomMeasurementCellBounds SpinDensityNew::getCustomMeasurementCellBounds() const
{
  CustomMeasurementCellBounds bounds;
  const InstanceParameters& parameters = instance_parameters_;
  bounds.lo                            = simulation_lattice_.toUnit(parameters.derived_parameters.corner);
  bounds.hi                            = bounds.lo;
  for (int j = 0; j < QMCT::DIM; ++j)
  {
    const QMCT::PosType axis_u = simulation_lattice_.toUnit(parameters.custom_measurement_lattice->Rv[j]);
    for (int d = 0; d < QMCT::DIM; ++d)
      if (axis_u[d] < 0.0)
        bounds.lo[d] += axis_u[d];
      else
        bounds.hi[d] += axis_u[d];
  }
  return bounds;
}

bool SpinDensityNew::getCustomMeasurementCellPointForPeriodicSimulation(const QMCT::PosType& position,
                                                                        const CustomMeasurementCellBounds& bounds,
                                                                        size_t& point) const
{
  static_assert(QMCT::DIM == 3, "SpinDensity supports three-dimensional cells only");
  QMCT::PosType simulation_u = simulation_lattice_.toUnit(position);

  std::array<int, 3> nlo{};
  std::array<int, 3> nhi{};
  for (int d = 0; d < QMCT::DIM; ++d)
    if (simulation_lattice_.BoxBConds[d])
    {
      simulation_u[d] -= std::floor(simulation_u[d]);
      nlo[d] = static_cast<int>(std::ceil(bounds.lo[d] - simulation_u[d]));
      nhi[d] = static_cast<int>(std::floor(bounds.hi[d] - simulation_u[d]));
    }

  const QMCT::PosType primary_image = simulation_lattice_.toCart(simulation_u);
  for (int nx = nlo[0]; nx <= nhi[0]; ++nx)
    for (int ny = nlo[1]; ny <= nhi[1]; ++ny)
      for (int nz = nlo[2]; nz <= nhi[2]; ++nz)
      {
        const QMCT::PosType image = primary_image + nx * simulation_lattice_.Rv[0] + ny * simulation_lattice_.Rv[1] +
            nz * simulation_lattice_.Rv[2];
        if (getCustomMeasurementCellPointForOpenSimulation(image, point))
          return true;
      }
  return false;
}

void SpinDensityNew::accumulateToData(size_t point, QMCT::RealType weight)
{
  if (data_locality_ == DataLocality::crowd)
  {
    data_[point] += weight;
  }
  else if (data_locality_ == DataLocality::queue)
  {
    data_.push_back(point);
    data_.push_back(weight);
  }
  else
  {
    throw std::runtime_error("You cannot accumulate to a SpinDensityNew with datalocality of this type");
  }
}

void SpinDensityNew::collect(const RefVector<OperatorEstBase>& type_erased_operator_estimators)
{
  if (data_locality_ == DataLocality::rank)
  {
    for (OperatorEstBase& crowd_oeb : type_erased_operator_estimators)
    {
      // This will throw a std::bad_cast in debug if the calling code hands the
      // wrong type erased operator_estimator type into here.
      // In release we don't want that overhead.
#ifndef NDEBUG
      auto& oeb = dynamic_cast<SpinDensityNew&>(crowd_oeb);
#else
      auto& oeb = static_cast<SpinDensityNew&>(crowd_oeb);
#endif
      auto& data = oeb.get_data();
      for (int id = 0; id < data.size(); id += 2)
      {
        // This is a smell
        size_t point{static_cast<size_t>(data[id])};
        const QMCT::RealType weight{data[id + 1]};
        data_[point] += weight;
        walkers_weight_ += weight;
      }
      oeb.zero();
    }
  }
  else if (data_locality_ == DataLocality::crowd)
  {
    OperatorEstBase::collect(type_erased_operator_estimators);
  }
  else
  {
    throw std::runtime_error("You cannot call collect on a SpinDensityNew with this DataLocality");
  }
}

void SpinDensityNew::report(const std::string& pad)
{
  const InstanceParameters& parameters = instance_parameters_;
  const auto& dp_                      = parameters.derived_parameters;
  app_log() << pad << "SpinDensity report" << std::endl;
  app_log() << pad << "  dim     = " << QMCT::DIM << std::endl;
  app_log() << pad << "  npoints = " << dp_.npoints << std::endl;
  const Lattice& measurement_lattice =
      parameters.custom_measurement_lattice ? *parameters.custom_measurement_lattice : simulation_lattice_;
  const QMCT::PosType measurement_corner =
      parameters.custom_measurement_lattice ? dp_.corner : simulation_lattice_.toCart(parameters.implicit_corner_u);
  app_log() << pad << "  grid    = " << dp_.grid << std::endl;
  app_log() << pad << "  gdims   = " << dp_.gdims << std::endl;
  app_log() << pad << "  corner  = " << measurement_corner << std::endl;
  app_log() << pad << "  center  = " << measurement_corner + measurement_lattice.Center << std::endl;
  app_log() << pad << "  cell " << std::endl;
  for (int d = 0; d < QMCT::DIM; ++d)
    app_log() << pad << "    " << d << " " << measurement_lattice.Rv[d] << std::endl;
  app_log() << pad << "  end cell " << std::endl;
  app_log() << pad << "  nspecies = " << species_.size() << std::endl;
  for (int s = 0; s < species_.size(); ++s)
    app_log() << pad << "    species[" << s << "]"
              << " = " << species_.speciesName[s] << " " << species_size_[s] << std::endl;
  app_log() << pad << "end SpinDensity report" << std::endl;
}

void SpinDensityNew::registerOperatorEstimator(hdf_archive& file)
{
  std::vector<size_t> my_indexes;

  std::vector<int> ng(1, instance_parameters_.derived_parameters.npoints);

  hdf_path hdf_name{my_name_};
  for (int s = 0; s < species_.size(); ++s)
  {
    h5desc_.emplace_back(hdf_name / species_.speciesName[s]);
    auto& oh = h5desc_.back();
    oh.set_dimensions(ng, s * instance_parameters_.derived_parameters.npoints);
  }
}


} // namespace qmcplusplus
