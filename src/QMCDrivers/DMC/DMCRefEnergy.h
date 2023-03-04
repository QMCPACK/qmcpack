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

namespace qmcplusplus
{
/** Handle updating Eref used for calculating the trial energy.
 */
class DMCRefEnergy
{
public:
  using FullPrecRealType = QMCTraits::FullPrecRealType;
  enum Scheme
  {
    LEGACY,
    LIMITED_HISTORY
  } scheme_;

  enum DataLayout
  {
    ENERGY = 0,
    VARIANCE,
    DATA_SIZE
  };

private:
  // legacy scheme data
  ///a simple accumulator for energy
  accumulator_set<FullPrecRealType> energy_hist_;
  ///a simple accumulator for variance
  accumulator_set<FullPrecRealType> variance_hist_;

  // limited memory scheme data
  SizeLimitedDataQueue<FullPrecRealType, DataLayout::DATA_SIZE> energy_and_variance_;

public:
  DMCRefEnergy(Scheme scheme, size_t history_limit);

  /// return energy and variance
  std::tuple<FullPrecRealType, FullPrecRealType> getEnergyVariance() const;

  /// record weight, energy and variance.
  void pushWeightEnergyVariance(FullPrecRealType weight, FullPrecRealType ene, FullPrecRealType var);

  /// return record count.
  size_t count() const;
};

} // namespace qmcplusplus
#endif
