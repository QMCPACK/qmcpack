//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////
#include "Estimators/MagnetizationDensity.h"
#include "Estimators/MagnetizationDensityInput.h"

namespace qmcplusplus
{

MagnetizationDensity::MagnetizationDensity(MagnetizationDensityInput&& minput,
					   const Lattice & lat):
OperatorEstBase(DataLocality::crowd), input_(minput), lattice_(lat)
{
  //Pull consistent corner, grids, etc., from already inititalized input.  
  MagnetizationDensityInput::DerivedParameters derived = minput.calculateDerivedParameters(lat);
  npoints_ = derived.npoints;
  grid_= derived.grid;
  gdims_ = derived.gdims;

  rcorner_ = derived.corner;
  center_  = rcorner_ + lattice_.Center;

  nsamples_ = input_.get_nsamples();
  integrator_ = input_.get_integrator();
  switch(integrator_)
  {
    case Integrator::SIMPSONS:
      app_log()<<"SIMPSONS\n";
      break;
    case Integrator::MONTECARLO:
      app_log()<<"MONTECARLO\n";
      break;
  }

  //Resize the data arrays.
  data_.resize(getFullDataSize());

}

MagnetizationDensity::MagnetizationDensity(const MagnetizationDensity& magdens,DataLocality dl): MagnetizationDensity(magdens)
{
 data_locality_ = dl;
}
void MagnetizationDensity::startBlock(int steps)
{
};
   
size_t MagnetizationDensity::getFullDataSize() { return npoints_ * DIM; }

void MagnetizationDensity::accumulate(const RefVector<MCPWalker>& walkers,
                  const RefVector<ParticleSet>& psets,
                  const RefVector<TrialWaveFunction>& wfns,
                  RandomGenerator& rng) 
{
  for (int iw = 0; iw < walkers.size(); ++iw)
  {
    MCPWalker& walker     = walkers[iw];
    ParticleSet& pset     = psets[iw];
    QMCT::RealType weight = walker.Weight;

    const int np = pset.getTotalNum();
    assert(weight >= 0);
    for (int p = 0; p < np; ++p)
    {
      QMCT::PosType u = lattice_.toUnit(pset.R[p] - rcorner_);
      size_t point    = 0;
      for (int d = 0; d < QMCT::DIM; ++d)
        point += gdims_[d] * ((int)(grid_[d] * (u[d] - std::floor(u[d])))); //periodic only
     // accumulateToData(point, weight);
    }
  }

};

std::unique_ptr<OperatorEstBase> MagnetizationDensity::spawnCrowdClone() const 
{ 
  std::size_t data_size    = data_.size();
  auto spawn_data_locality = data_locality_;

  //Everyone else has this attempt to set up a non-implemented memory saving optimization.
  //We won't rock the boat.  
  if (data_locality_ == DataLocality::rank)
  {
    // This is just a stub until a memory saving optimization is deemed necessary
    spawn_data_locality = DataLocality::queue;
    data_size           = 0;
    throw std::runtime_error("There is no memory savings implementation for MagnetizationDensity");
  }

  UPtr<MagnetizationDensity> spawn(std::make_unique<MagnetizationDensity>(*this, spawn_data_locality));
  spawn->get_data().resize(data_size);
  return spawn;
};

void MagnetizationDensity::registerOperatorEstimator(hdf_archive& file)
{
  std::vector<size_t> my_indexes;

  std::vector<int> ng(DIM*npoints_);

  hdf_path hdf_name{my_name_};
//  for (int s = 0; s < species_.size(); ++s)
//  {
  h5desc_.emplace_back(hdf_name);
  auto& oh = h5desc_.back();
  oh.set_dimensions(ng, 0);
//  }
}

} //namespace qmcplusplus
