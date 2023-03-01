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
#ifndef QMCPLUSPLUS_SIZELIMITEDDATAQUEUE_H
#define QMCPLUSPLUS_SIZELIMITEDDATAQUEUE_H

#include <deque>
#include <array>
#include <cassert>

namespace qmcplusplus
{
template<typename T, size_t NUM_DATA_FIELDS>
class SizeLimitedDataQueue
{
private:
  enum
  {
    WEIGHT         = 0,
    FIRST_PROPERTY = 1,
    NUM_FIELDS     = NUM_DATA_FIELDS + 1
  };
  using value_type = std::array<T, NUM_FIELDS>;
  std::deque<value_type> data;
  const size_t size_limit_;

public:
  SizeLimitedDataQueue(size_t size_limit) : size_limit_(size_limit) {}

  void push(const value_type& val)
  {
    if (data.size() == size_limit_)
      data.pop_front();
    assert(data.size() < size_limit_);
    data.push_back(val);
  }

  void push(value_type&& val)
  {
    if (data.size() == size_limit_)
      data.pop_front();
    assert(data.size() < size_limit_);
    data.push_back(val);
  }

  value_type weighted_avg() const
  {
    value_type avg;
    std::fill(avg.begin(), avg.end(), T(0));
    for (auto& element : data)
    {
      T weight = element[WEIGHT];
      avg[WEIGHT] += weight;
      for (size_t i = FIRST_PROPERTY; i < element.size(); i++)
        avg[i] += element[i] * weight;
    }
    T total_weight = avg[WEIGHT];
    avg[WEIGHT] /= data.size();
    for (size_t i = FIRST_PROPERTY; i < avg.size(); i++)
      avg[i] /= total_weight;
    return avg;
  }

  auto size() const { return data.size(); }
};

} // namespace qmcplusplus
#endif
