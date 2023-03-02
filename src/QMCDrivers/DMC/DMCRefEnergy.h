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
  WEIGHT = 0,
  ENERGY,
  VARIANCE,
  DATA_SIZE
};

private:
  // legacy scheme
  ///a simple accumulator for energy
  accumulator_set<FullPrecRealType> energy_hist_;
  ///a simple accumulator for variance
  accumulator_set<FullPrecRealType> variance_hist_;

  // limited memory scheme
  SizeLimitedDataQueue<FullPrecRealType, DataLayout::DATA_SIZE> energy_and_variance_;

  public:
  DMCRefEnergy(Scheme scheme, size_t history_limit);

  std::tuple<FullPrecRealType, FullPrecRealType> getEnergyVariance() const;

  void pushWeightEnergyVariance(FullPrecRealType weight, FullPrecRealType ene, FullPrecRealType var);

  size_t count() const;

};

} // namespace qmcplusplus
#endif
