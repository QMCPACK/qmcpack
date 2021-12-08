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

#include <iostream>
#include <numeric>
#include <SpeciesSet.h>

namespace qmcplusplus
{
SpinDensityNew::SpinDensityNew(SpinDensityInput&& input, const SpeciesSet& species, DataLocality dl)
    : OperatorEstBase(dl), input_(std::move(input)), species_(species), species_size_(getSpeciesSize(species))
{
  my_name_ = "SpinDensity";

  if (input_.get_cell().explicitly_defined == true)
    lattice_ = input_.get_cell();
  else
    throw std::runtime_error("If SpinDensityInput does not contain a cell definition you must call the constructor "
                             "with an explicit lattice defined");

  derived_parameters_ = input_.calculateDerivedParameters(lattice_);

  data_.resize(getFullDataSize(), 0.0);

  if (input_.get_write_report())
    report("  ");
}

SpinDensityNew::SpinDensityNew(SpinDensityInput&& input,
                               const Lattice& lattice,
                               const SpeciesSet& species,
                               const DataLocality dl)
    : OperatorEstBase(dl),
      input_(std::move(input)),
      species_(species),
      species_size_(getSpeciesSize(species)),
      lattice_(lattice)
{
  my_name_ = "SpinDensity";
  std::cout << "SpinDensity constructor called\n";
  data_locality_ = dl;
  if (input_.get_cell().explicitly_defined == true)
    throw std::runtime_error(
        "SpinDensityNew should not be constructed with both a cell in its input and an lattice input arguement.");
  else if (lattice_.explicitly_defined == false)
    throw std::runtime_error("SpinDensityNew cannot be constructed from a lattice that is not explicitly defined");

  derived_parameters_ = input_.calculateDerivedParameters(lattice_);
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

size_t SpinDensityNew::getFullDataSize() { return species_.size() * derived_parameters_.npoints; }

std::unique_ptr<OperatorEstBase> SpinDensityNew::spawnCrowdClone() const {
  std::size_t data_size = data_.size();
  auto spawn_data_locality = data_locality_;
  if (data_locality_ == DataLocality::rank)
  {
    spawn_data_locality = DataLocality::queue;
    // at construction we don't know what the data requirement is going to be
    // since its steps per block  dependent. so start with 10 steps worth.
    int num_particles = std::accumulate(species_size_.begin(), species_size_.end(), 0);
    data_size  = num_particles * 20;
  }
  UPtr<SpinDensityNew> spawn(std::make_unique<SpinDensityNew>(*this, spawn_data_locality));
  spawn->get_data().resize(data_size);
  return spawn;
}

void SpinDensityNew::startBlock(int steps)
{
  if (data_locality_ == DataLocality::rank)
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
                                RandomGenerator_t& rng)
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
    int p                             = 0;
    size_t offset                     = 0;
    for (int s = 0; s < species_.size(); ++s, offset += dp_.npoints)
      for (int ps = 0; ps < species_size_[s]; ++ps, ++p)
      {
        QMCT::PosType u = lattice_.toUnit(pset.R[p] - dp_.corner);
        size_t point    = offset;
        for (int d = 0; d < QMCT::DIM; ++d)
          point += dp_.gdims[d] * ((int)(dp_.grid[d] * (u[d] - std::floor(u[d])))); //periodic only
        accumulateToData(point, weight);
      }
  }
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

void SpinDensityNew::registerOperatorEstimator(hid_t gid)
{
  std::vector<size_t> my_indexes;
  hid_t sgid = H5Gcreate(gid, my_name_.c_str(), 0);

  //vector<int> ng(DIM);
  //for(int d=0;d<DIM;++d)
  //  ng[d] = grid[d];

  std::vector<int> ng(1);
  ng[0] = derived_parameters_.npoints;

  for (int s = 0; s < species_.size(); ++s)
  {
    h5desc_.emplace_back(std::make_unique<ObservableHelper>(species_.speciesName[s]));
    auto& oh = h5desc_.back();
    oh->set_dimensions(ng, 0);
    oh->open(sgid);
    // bad smell
  }
}


} // namespace qmcplusplus
