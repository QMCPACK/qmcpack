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

namespace qmcplusplus
{

template<typename T>
class EstimatorBuffer
{
public:
  EstimatorBuffer() = default;

  template<class... ARGS>
  std::size_t copyIn(T weight, ARGS&... args)
  {
    return totalSize(std::forward<ARGS>(args)...);
  }


private:
  
  std::size_t totalSize(const std::vector<T>& vec) { return vec.size(); }
  template<typename... ARGS>
  std::size_t totalSize(const std::vector<T>& vec, ARGS... vecs)
  {
    return vec.size() + totalSize(std::forward<ARGS>(vecs)...);
  }


  std::vector<T> data_;
};

} // namespace qmcplusplus
#endif
