//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include <vector>
#include <type_traits>

#ifndef QMCPLUSPLUS_ESTIMATOR_BUFFER_HPP
#define QMCPLUSPLUS_ESTIMATOR_BUFFER_HPP

#include "PooledData.h"

namespace qmcplusplus
{

template<typename T>
class EstimatorBuffer
{
public:
  EstimatorBuffer() = default;
  void reserve(std::size_t);
private:

  std::vector<T> data_;
};

} // namespace qmcplusplus
#endif
