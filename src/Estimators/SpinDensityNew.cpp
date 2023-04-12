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

#include <iostream>
#include <numeric>
#include <SpeciesSet.h>

namespace qmcplusplus
{
SpinDensityNew::SpinDensityNew(SpinDensityInput&& input, const SpeciesSet& species, DataLocality dl)
    : OperatorEstBase(dl), input_(std::move(input)), species_(species), species_size_(getSpeciesSize(species))
{
  my_name_ = "SpinDensity";

  data_locality_ = DataLocality::crowd;
  if (input_.get_save_memory())
    dl = DataLocality::rank;

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
  my_name_       = "SpinDensity";
  data_locality_ = dl;
  if (input_.get_cell().explicitly_defined == true)
    lattice_ = input_.get_cell();
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
                                RandomGenerator& rng)
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
    QMCT::PosType u;
    for (int s = 0; s < species_.size(); ++s, offset += dp_.npoints)
      for (int ps = 0; ps < species_size_[s]; ++ps, ++p)
      {
        switch (pset.R.getUnit())
        {
        case PosUnit::Cartesian:
          u = lattice_.toUnit(pset.R[p] - dp_.corner);
          break;
        case PosUnit::Lattice:
#ifndef NDEBUG
          if (!(p < pset.getTotalNum()))
            throw std::runtime_error("p < pset.getTotalNum(): " + std::to_string(p) + " " +
                                     std::to_string(pset.getTotalNum()));
#endif
          u = pset.R[p] - lattice_.toUnit(dp_.corner);
          break;
        }
        size_t point = offset;
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

void SpinDensityNew::report(const std::string& pad) const { report(pad, app_log()); }

void SpinDensityNew::report(const std::string& pad, std::ostream& out) const
{
  auto& dp_ = derived_parameters_;
  out << pad << "SpinDensity report" << std::endl;
  out << pad << "  dim     = " << QMCT::DIM << std::endl;
  out << pad << "  npoints = " << dp_.npoints << std::endl;
  out << pad << "  grid    = " << dp_.grid << std::endl;
  out << pad << "  gdims   = " << dp_.gdims << std::endl;
  out << pad << "  corner  = " << dp_.corner << std::endl;
  out << pad << "  center  = " << dp_.corner + lattice_.Center << std::endl;
  out << pad << "  cell " << std::endl;
  for (int d = 0; d < QMCT::DIM; ++d)
    out << pad << "    " << d << " " << lattice_.Rv[d] << std::endl;
  out << pad << "  end cell " << std::endl;
  out << pad << "  nspecies = " << species_.size() << std::endl;
  for (int s = 0; s < species_.size(); ++s)
    out << pad << "    species[" << s << "]"
        << " = " << species_.speciesName[s] << " " << species_size_[s] << std::endl;
  out << pad << "end SpinDensity report" << std::endl;
}


void SpinDensityNew::registerOperatorEstimator(hdf_archive& file) {}

void SpinDensityNew::write(hdf_archive& file)
{
  auto writeFullPrecData = [&](auto& file, auto& fp_data) {
    std::vector<size_t> my_indexes;
    std::array<hsize_t, 2> ng{1, derived_parameters_.npoints};
    hdf_path hdf_name{my_name_};
    file.push(hdf_name);
    Vector<QMCT::FullPrecRealType> data;
    std::size_t offset = 0;
    for (int s = 0; s < species_.size(); ++s)
    {
      data.attachReference(fp_data.data() + offset, derived_parameters_.npoints);
      file.push(species_.speciesName[s]);
      //file.write(data, "value");
      file.writeSlabReshaped(data, ng, "value");
      offset += derived_parameters_.npoints;
      file.pop();
    }
    file.pop();
  };

#ifdef MIXED_PRECISION
  std::vector<QMCT::FullPrecRealType> expanded_data(data_.size(), 0.0);
  std::copy_n(data_.begin(), data_.size(), expanded_data.begin());
  assert(!data_.empty());
  writeFullPrecData(file, expanded_data);
  // auto total = std::accumulate(data_->begin(), data_->end(), 0.0);
  // std::cout << "data size: " << data_->size() << " : " << total << '\n';

#else
  writeFullPrecData(file, data_);
#endif
}

} // namespace qmcplusplus
