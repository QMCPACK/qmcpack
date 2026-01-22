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

#ifndef QMCPLUSPLUS_UNINITIALIZED_ARRAY_HPP
#define QMCPLUSPLUS_UNINITIALIZED_ARRAY_HPP

#include <cstddef>

#if defined(__CUDACC__) || defined(__HIPCC__)
#  define QMCPACK_HOST_DEVICE __host__ __device__
#else
#  define QMCPACK_HOST_DEVICE
#endif

namespace qmcplusplus
{
namespace device
{
template <class T, std::size_t N>
struct uninitialized_array
{
  using value_type = T;
  static constexpr std::size_t size = N;
  alignas(T) unsigned char data_[N * sizeof(T)];

  QMCPACK_HOST_DEVICE T* data()
  {
    return reinterpret_cast<T*>(data_);
  }

  QMCPACK_HOST_DEVICE const T* data() const
  {
    return reinterpret_cast<const T*>(data_);
  }

  QMCPACK_HOST_DEVICE T& operator[](unsigned int idx)
  {
    return data()[idx];
  }

  QMCPACK_HOST_DEVICE T const& operator[](unsigned int idx) const
  {
    return data()[idx];
  }
};
} // namespace device
} // namespace qmcplusplus
#undef QMCPACK_HOST_DEVICE
#endif
