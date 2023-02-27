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
