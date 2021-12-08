//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_ACCESS_TRAITS_H
#define QMCPLUSPLUS_ACCESS_TRAITS_H

#include <algorithm>

namespace qmcplusplus
{
/** template class analogous to std::allocator_traits.
 *  * defines the is_host_accessible and is_dual_space traits
 *  * abstracts the data movement on,off and from place to place on device.
 *  * abstracts the fill function for the allocator.
 */
template<class Allocator>
struct qmc_allocator_traits
{
  using value_type = typename Allocator::value_type;

  static const bool is_host_accessible = true;
  static const bool is_dual_space      = false;

  static void fill_n(value_type* ptr, size_t n, const value_type& value) { std::fill_n(ptr, n, value); }

  // So we can write generic tests that work with all QMCPACK allocators
  static void attachReference(Allocator& from, Allocator& to, value_type* from_data, value_type* ref) {}
  // These abstract synchronous transfers, async semantics are vender specific
  static void updateTo(Allocator& a, value_type* host_ptr, size_t n) {}
  static void updateFrom(Allocator& a, value_type* host_ptr, size_t n) {}
  static void deviceSideCopyN(Allocator& a, size_t to, size_t n, size_t from) {}
};

template<class Allocator>
using IsHostSafe = std::enable_if_t<qmc_allocator_traits<Allocator>::is_host_accessible>;

template<class Allocator>
using IsNotHostSafe = std::enable_if_t<!qmc_allocator_traits<Allocator>::is_host_accessible>;

template<class Allocator>
using IsDualSpace = std::enable_if_t<qmc_allocator_traits<Allocator>::is_dual_space>;
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_ACCESS_TRAITS_H
