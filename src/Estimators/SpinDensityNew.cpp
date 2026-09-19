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
      lattice_(input_.has_cell() ? input_.get_cell() : lattice),
      simulation_lattice_(lattice)
{
  data_locality_      = dl;
  derived_parameters_ = input_.calculateDerivedParameters(lattice_);
  if (input_.has_cell())
    initializeFiniteCellBounds();
  data_.resize(getFullDataSize());
  if (input_.get_write_report())
    report("  ");
}

SpinDensityNew::SpinDensityNew(const SpinDensityNew& sdn, DataLocality dl) : SpinDensityNew(sdn)
{
  data_locality_ = dl;
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

size_t SpinDensityNew::getFullDataSize() const { return species_.size() * derived_parameters_.npoints; }

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

/** Gets called every step and writes to thread local data.
 *
 *  I tried for readable and not doing the optimizers job.
 *  The offsets into bare data are already bad enough.
 */
void SpinDensityNew::accumulate(const RefVector<MCPWalker>& walkers,
                                const RefVector<ParticleSet>& psets,
                                const RefVector<TrialWaveFunction>& wfns,
                                const RefVector<QMCHamiltonian>& hams,
                                RandomBase<FullPrecRealType>& rng)
{
  auto& dp_ = derived_parameters_;
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
    for (int s = 0; s < species_.size(); ++s, offset += dp_.npoints)
      for (int ps = 0; ps < species_size_[s]; ++ps, ++p)
      {
        size_t point = offset;
        if (!input_.has_cell())
        {
          const QMCT::PosType u = lattice_.toUnit(pset.R[p] - dp_.corner);
          for (int d = 0; d < QMCT::DIM; ++d)
            point += dp_.gdims[d] * static_cast<int>(dp_.grid[d] * (u[d] - std::floor(u[d])));
          accumulateToData(point, weight);
        }
        else if (simulation_lattice_.SuperCellEnum == SUPERCELL_OPEN)
        {
          if (getFiniteCellPoint(pset.R[p], point))
            accumulateToData(point, weight);
        }
        else if (getPeriodicFiniteCellPoint(pset.R[p], point))
          accumulateToData(point, weight);
      }
  }
}

bool SpinDensityNew::getFiniteCellPoint(const QMCT::PosType& position, size_t& point) const
{
  const QMCT::PosType u  = lattice_.toUnit(position - derived_parameters_.corner);
  size_t candidate_point = point;
  for (int d = 0; d < QMCT::DIM; ++d)
  {
    if (!std::isfinite(u[d]) || u[d] < 0.0 || u[d] >= 1.0)
      return false;
    candidate_point += derived_parameters_.gdims[d] *
        std::min(static_cast<int>(derived_parameters_.grid[d] * u[d]), derived_parameters_.grid[d] - 1);
  }
  point = candidate_point;
  return true;
}

void SpinDensityNew::initializeFiniteCellBounds()
{
  finite_cell_lo_ = simulation_lattice_.toUnit(derived_parameters_.corner);
  finite_cell_hi_ = finite_cell_lo_;
  for (int j = 0; j < QMCT::DIM; ++j)
  {
    const QMCT::PosType axis_u = simulation_lattice_.toUnit(lattice_.Rv[j]);
    for (int d = 0; d < QMCT::DIM; ++d)
      if (axis_u[d] < 0.0)
        finite_cell_lo_[d] += axis_u[d];
      else
        finite_cell_hi_[d] += axis_u[d];
  }
}

bool SpinDensityNew::getPeriodicFiniteCellPoint(const QMCT::PosType& position, size_t& point) const
{
  static_assert(QMCT::DIM == 3, "SpinDensity supports three-dimensional cells only");
  QMCT::PosType simulation_u = simulation_lattice_.toUnit(position);

  std::array<int, 3> nlo{};
  std::array<int, 3> nhi{};
  for (int d = 0; d < QMCT::DIM; ++d)
    if (simulation_lattice_.BoxBConds[d])
    {
      simulation_u[d] -= std::floor(simulation_u[d]);
      nlo[d] = static_cast<int>(std::ceil(finite_cell_lo_[d] - simulation_u[d]));
      nhi[d] = static_cast<int>(std::floor(finite_cell_hi_[d] - simulation_u[d]));
    }

  const QMCT::PosType primary_image = simulation_lattice_.toCart(simulation_u);
  for (int nx = nlo[0]; nx <= nhi[0]; ++nx)
    for (int ny = nlo[1]; ny <= nhi[1]; ++ny)
      for (int nz = nlo[2]; nz <= nhi[2]; ++nz)
      {
        const QMCT::PosType image = primary_image + nx * simulation_lattice_.Rv[0] + ny * simulation_lattice_.Rv[1] +
            nz * simulation_lattice_.Rv[2];
        if (getFiniteCellPoint(image, point))
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
  auto& dp_ = derived_parameters_;
  app_log() << pad << "SpinDensity report" << std::endl;
  app_log() << pad << "  dim     = " << QMCT::DIM << std::endl;
  app_log() << pad << "  npoints = " << dp_.npoints << std::endl;
  app_log() << pad << "  grid    = " << dp_.grid << std::endl;
  app_log() << pad << "  gdims   = " << dp_.gdims << std::endl;
  app_log() << pad << "  corner  = " << dp_.corner << std::endl;
  app_log() << pad << "  center  = " << dp_.corner + lattice_.Center << std::endl;
  app_log() << pad << "  cell " << std::endl;
  for (int d = 0; d < QMCT::DIM; ++d)
    app_log() << pad << "    " << d << " " << lattice_.Rv[d] << std::endl;
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

  std::vector<int> ng(1, derived_parameters_.npoints);

  hdf_path hdf_name{my_name_};
  for (int s = 0; s < species_.size(); ++s)
  {
    h5desc_.emplace_back(hdf_name / species_.speciesName[s]);
    auto& oh = h5desc_.back();
    oh.set_dimensions(ng, s * derived_parameters_.npoints);
  }
}


} // namespace qmcplusplus
