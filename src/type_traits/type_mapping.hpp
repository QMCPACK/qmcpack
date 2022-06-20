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

template<typename V1, typename V2, typename T>
struct OnTypesEqual : std::bool_constant<std::is_same<V1, V2>::value>
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
