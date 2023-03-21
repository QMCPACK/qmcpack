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

/** collect data with a history limit.
 * data stored in std::deque<std::array<T, NUM_FIELDS>>
 */
template<typename T, size_t NUM_FIELDS>
class SizeLimitedDataQueue
{
public:
  struct HistoryElement
  {
    T weight;
    std::array<T, NUM_FIELDS> properties;
  };

  using value_type = HistoryElement;

  SizeLimitedDataQueue(size_t size_limit) : size_limit_(size_limit) {}

  /// add a new record
  void push(const value_type& val)
  {
    if (data.size() == size_limit_)
      data.pop_front();
    assert(data.size() < size_limit_);
    data.push_back(val);
  }

  /// add a new record
  void push(value_type&& val)
  {
    if (data.size() == size_limit_)
      data.pop_front();
    assert(data.size() < size_limit_);
    data.push_back(val);
  }

  /// return weighted average
  auto weighted_avg() const
  {
    std::array<T, NUM_FIELDS> avg;
    std::fill(avg.begin(), avg.end(), T(0));
    T weight_sum = 0;
    for (auto& element : data)
    {
      weight_sum += element.weight;
      for (size_t i = 0; i < NUM_FIELDS; i++)
        avg[i] += element.properties[i] * element.weight;
    }
    for (size_t i = 0; i < NUM_FIELDS; i++)
      avg[i] /= weight_sum;
    return avg;
  }

  /// return the number of records
  auto size() const { return data.size(); }

private:
  std::deque<value_type> data;
  const size_t size_limit_;
};

} // namespace qmcplusplus
#endif
