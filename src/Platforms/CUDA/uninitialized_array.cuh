//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Youngjun Lee, leey@anl.gov, Argonne National Laboratory
//
// File created by: Youngjun Lee, leey@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_UNINITIALIZED_ARRAY_CUH
#define QMCPLUSPLUS_UNINITIALIZED_ARRAY_CUH

#include <cstddef>

namespace qmcplusplus
{
namespace device
{
template<class T, std::size_t N>
struct uninitialized_array
{
  using value_type                  = T;
  static constexpr std::size_t size = N;
  alignas(T) unsigned char data_[N * sizeof(T)];

  __device__ T* data() { return reinterpret_cast<T*>(data_); }

  __device__ const T* data() const { return reinterpret_cast<const T*>(data_); }

  __device__ T& operator[](unsigned int idx) { return data()[idx]; }

  __device__ T const& operator[](unsigned int idx) const { return data()[idx]; }
};
} // namespace device
} // namespace qmcplusplus
#endif
