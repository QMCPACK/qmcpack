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

#include <iostream>
#include "SpinDensityNew.h"

namespace qmcplusplus
{
SpinDensityNew::SpinDensityNew(SpinDensityInput&& input, const SpeciesSet& species) : input_(std::move(input)), species_(species)
{
  myName = "SpinDensity";
  std::cout << "SpinDensity constructor called\n";
  // This code is quite suspect.
  // I think it is checking the membersize is either the last attribute or adding it.  If its already there but not
  // not last it fails. Not sure why we care yet.
  int index = species_.findAttribute("membersize");
  // That there will be a number of particles of a particular species is an invariant
  // but the SpeciesSet fails to say that so we have this
  if (index < 0)
    throw std::runtime_error("SpinDensity(P) Species set does not have the required attribute 'membersize'");
  for (int s = 0; s < species_.size(); ++s)
    species_size_.push_back(species_(index, s));
  size_t size_data = species_.size() * input_.get_npoints();
  data_            = createLocalData(size_data, data_locality_);
}

OperatorEstBase* SpinDensityNew::clone()
{
  std::cout << "SpinDensity clone called\n";
  return new SpinDensityNew(*this);
}

SpinDensityNew::SpinDensityNew(const SpinDensityNew& sdn)
    : input_(sdn.input_), species_(sdn.species_), species_size_(sdn.species_size_)
{
  data_locality_ = sdn.data_locality_;

  std::cout << "SpinDensity copy constructor called\n";
  size_t data_size = sdn.data_->size();
  data_            = createLocalData(data_size, data_locality_);
}


/** Gets called every step and writes to likely thread local data.
 *
 *  I tried for readable and not doing the optimizers job.
 *  The offsets into bare data are already bad enough.
 */
void SpinDensityNew::accumulate(RefVector<MCPWalker>& walkers, RefVector<ParticleSet>& psets)
{
  auto accum_per_walker = [this](MCPWalker& walker, ParticleSet& pset) {
    QMCT::RealType weight = walker.Weight;
    // for testing
    walkers_weight_ += weight;
    int p                             = 0;
    std::vector<QMCT::RealType>& data = *data_;
    for (int s = 0; s < species_.size(); ++s)
      for (int ps = 0; ps < species_size_[s]; ++ps, ++p)
      {
        QMCT::PosType u = input_.get_cell().toUnit(pset.R[p] - input_.get_corner());
        size_t point    = input_.get_npoints() * s;
        for (int d = 0; d < QMCT::DIM; ++d)
          point += input_.get_gdims()[d] * ((int)(input_.get_grid()[d] * (u[d] - std::floor(u[d])))); //periodic only
        data[point] += weight;
      }
  };

  for (int iw = 0; iw < walkers.size(); ++iw)
  {
    accum_per_walker(walkers[iw], psets[iw]);
  }
};

void SpinDensityNew::report(const std::string& pad)
{
  app_log() << pad << "SpinDensity report" << std::endl;
  app_log() << pad << "  dim     = " << QMCT::DIM << std::endl;
  app_log() << pad << "  npoints = " << input_.get_npoints() << std::endl;
  app_log() << pad << "  grid    = " << input_.get_grid() << std::endl;
  app_log() << pad << "  gdims   = " << input_.get_gdims() << std::endl;
  app_log() << pad << "  corner  = " << input_.get_corner() << std::endl;
  app_log() << pad << "  center  = " << input_.get_corner() + input_.get_cell().Center << std::endl;
  app_log() << pad << "  cell " << std::endl;
  for (int d = 0; d < QMCT::DIM; ++d)
    app_log() << pad << "    " << d << " " << input_.get_cell().Rv[d] << std::endl;
  app_log() << pad << "  end cell " << std::endl;
  app_log() << pad << "  nspecies = " << species_.size() << std::endl;
  for (int s = 0; s < species_.size(); ++s)
    app_log() << pad << "    species[" << s << "]"
              << " = " << species_.speciesName[s] << " " << species_.attribName.size() << std::endl;
  app_log() << pad << "end SpinDensity report" << std::endl;
}

void SpinDensityNew::registerOperatorEstimator(std::vector<observable_helper*>& h5desc, hid_t gid) const
{
  hid_t sgid = H5Gcreate(gid, myName.c_str(), 0);

  //vector<int> ng(DIM);
  //for(int d=0;d<DIM;++d)
  //  ng[d] = grid[d];

  std::vector<int> ng(1);
  ng[0] = input_.get_npoints();

  for (int s = 0; s < species_.size(); ++s)
  {
    observable_helper* oh = new observable_helper(species_.speciesName[s]);
    oh->set_dimensions(ng, 0);
    oh->open(sgid);
    h5desc.push_back(oh);
  }
}

} // namespace qmcplusplus
