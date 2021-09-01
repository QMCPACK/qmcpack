//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: MomentumEstimator.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "MomentumDistribution.h"

#include <iostream>
#include <numeric>

namespace qmcplusplus
{
MomentumDistribution::MomentumDistribution(MomentumDistributionInput&& mdi, DataLocality dl)
    : OperatorEstBase(dl), input_(std::move(mdi))
{
  myName         = "MomentumDistribution";
  data_locality_ = dl;
}


MomentumDistribution* MomentumDistribution::clone()
{
  std::cout << "MomentumDistribution clone called\n";
  return new MomentumDistribution(std::move(input_),data_locality_);
}

MomentumDistribution::MomentumDistribution(const MomentumDistribution& md)
    : OperatorEstBase(md),
      input_(std::move(md.input_))
{
  if (data_locality_ == DataLocality::crowd)
  {
    size_t data_size = md.data_->size();
    data_            = createLocalData(data_size, data_locality_);
  }
  else if (data_locality_ == DataLocality::rank)
  {
    assert(md.data_locality_ == DataLocality::rank);
    size_t data_size  = 10; // jtk fix
    data_locality_    = DataLocality::queue;
    data_             = createLocalData(data_size, data_locality_);
  }
}

void MomentumDistribution::startBlock(int steps)
{
  if (data_locality_ == DataLocality::rank)
  {
    size_t data_size  = 10; // jtk fix
    data_->reserve(data_size);
    data_->resize(0);
  }
}

/** Gets called every step and writes to thread local data.
 *
 */
void MomentumDistribution::accumulate(const RefVector<MCPWalker>& walkers, const RefVector<ParticleSet>& psets)
{
};


void MomentumDistribution::collect(const RefVector<OperatorEstBase>& type_erased_operator_estimators)
{
  if (data_locality_ == DataLocality::crowd)
  {
    OperatorEstBase::collect(type_erased_operator_estimators);
  }
  else
  {
    throw std::runtime_error("You cannot call collect on a MomentumDistribution with this DataLocality");
  }
}


void MomentumDistribution::registerOperatorEstimator(hid_t gid)
{
}


} // namespace qmcplusplus
