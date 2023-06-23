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
// -*- C++ -*-
#ifndef QMCPLUSPLUS_DMCREFENERGY_H
#define QMCPLUSPLUS_DMCREFENERGY_H

#include <tuple>
#include <Configuration.h>
#include <Estimators/SizeLimitedDataQueue.hpp>
#include <Estimators/accumulators.h>
#include "DMCRefEnergyScheme.h"

namespace qmcplusplus
{
/** Handle updating Eref used for calculating the trial energy.
 */
class DMCRefEnergy
{
public:
  using FullPrecReal = QMCTraits::FullPrecRealType;

  enum DataLayout
  {
    ENERGY = 0,
    VARIANCE,
    DATA_SIZE
  };

private:
  /// scheme
  DMCRefEnergyScheme scheme_;

  // legacy scheme data
  ///a simple accumulator for energy
  accumulator_set<FullPrecReal> energy_hist_;
  ///a simple accumulator for variance
  accumulator_set<FullPrecReal> variance_hist_;

  // limited memory scheme data
  SizeLimitedDataQueue<FullPrecReal, DataLayout::DATA_SIZE> energy_and_variance_;

public:
  DMCRefEnergy(DMCRefEnergyScheme scheme, size_t history_limit);

  /// return energy and variance
  std::tuple<FullPrecReal, FullPrecReal> getEnergyVariance() const;

  /// record weight, energy and variance.
  void pushWeightEnergyVariance(FullPrecReal weight, FullPrecReal ene, FullPrecReal var);

  /// return record count.
  size_t count() const;
};

} // namespace qmcplusplus
#endif
