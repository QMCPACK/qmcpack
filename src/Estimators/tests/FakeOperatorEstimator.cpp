//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "FakeOperatorEstimator.h"
#include "type_traits/DataLocality.h"

namespace qmcplusplus
{

  FakeOperatorEstimator::FakeOperatorEstimator(int num_ranks, DataLocality data_locality) :
    OperatorEstBase(data_locality)
{
  data_locality_ = data_locality;
  data_.resize(num_ranks * 10);
}

FakeOperatorEstimator::FakeOperatorEstimator(const FakeOperatorEstimator& foe)
  : OperatorEstBase(foe)
{
  data_.resize(foe.data_.size());
}

}
