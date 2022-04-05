//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_TYPE_MAPPING_H
#define QMCPLUSPLUS_TYPE_MAPPING_H

#include <type_traits>

namespace qmcplusplus
{
// This namespace is to protect against the clash between this disjuntion and the std one
// when c++17 is finally enabled for QMCPACK.  When 14 support is dropped confirm this can be dropped and move
// the the std one.

namespace c14disjunction
{
// std::disjuntion is in std::c++17 so remove this when c++14 support is dropped
template<class...>
struct disjunction : std::false_type
{};
template<class B1>
struct disjunction<B1> : B1
{};
template<class B1, class... Bn>
struct disjunction<B1, Bn...> : std::conditional_t<bool(B1::value), B1, disjunction<Bn...>>
{};
} // namespace c14disjunction

template<typename V1, typename V2, typename T>
struct OnTypesEqual : std::integral_constant<bool, std::is_same<V1, V2>::value>
{
  using type = T;
};

template<typename T>
struct default_type : std::true_type
{
  using type = T;
};
} // namespace qmcplusplus

#endif
