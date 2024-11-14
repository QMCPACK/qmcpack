//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include "DMCRefEnergy.h"

namespace qmcplusplus
{

using FullPrecReal = DMCRefEnergy::FullPrecReal;

DMCRefEnergy::DMCRefEnergy(DMCRefEnergyScheme scheme, size_t history_limit)
    : scheme_(scheme), energy_and_variance_(history_limit)
{}

std::tuple<FullPrecReal, FullPrecReal> DMCRefEnergy::getEnergyVariance() const
{
  if (scheme_ == DMCRefEnergyScheme::LIMITED_HISTORY)
  {
    auto avg = energy_and_variance_.weighted_avg();
    return {avg[ENERGY], avg[VARIANCE]};
  }
  else
    return {energy_hist_.mean(), variance_hist_.mean()};
}

void DMCRefEnergy::pushWeightEnergyVariance(FullPrecReal weight, FullPrecReal ene, FullPrecReal var)
{
  if (scheme_ == DMCRefEnergyScheme::LIMITED_HISTORY)
    energy_and_variance_.push({weight, {ene, var}});
  else
  {
    energy_hist_(ene);
    variance_hist_(var);
  }
}

size_t DMCRefEnergy::count() const
{
  if (scheme_ == DMCRefEnergyScheme::LIMITED_HISTORY)
    return energy_and_variance_.size();
  else
  {
    assert(energy_hist_.count() == variance_hist_.count());
    return energy_hist_.count();
  }
}

} // namespace qmcplusplus
