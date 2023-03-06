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

MagnetizationDensity::MagnetizationDensity(DataLocality dl):OperatorEstBase(dl)
{
};

MagnetizationDensity::MagnetizationDensity(MagnetizationDensityInput&& minput,
					   const Lattice & lat):
OperatorEstBase(DataLocality::crowd), input_(minput), lattice_(lat)
{

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
  if (input_.get_corner_defined())
  {
    rcorner_ = input_.get_corner();
    center_  = rcorner_ + lattice_.Center;
  }
  else
  {
    if (input_.get_center_defined())
      center_ = input_.get_center();
    else
      center_ = lattice_.Center;
    rcorner_ = center_ -  lattice_.Center;
  }

}
void MagnetizationDensity::startBlock(int steps)
{
};

void MagnetizationDensity::accumulate(const RefVector<MCPWalker>& walkers,
                  const RefVector<ParticleSet>& psets,
                  const RefVector<TrialWaveFunction>& wfns,
                  RandomGenerator& rng) 
{

};

std::unique_ptr<OperatorEstBase> MagnetizationDensity::spawnCrowdClone() const 
{ 
  return 0;
};

} //namespace qmcplusplus
