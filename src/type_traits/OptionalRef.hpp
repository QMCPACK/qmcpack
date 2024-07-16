//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_OPTIONALREF_HPP
#define QMCPLUSPLUS_OPTIONALREF_HPP

#include <functional>
#include <optional>

namespace qmcplusplus
{
template<typename T>
using OptionalRef = std::optional<std::reference_wrapper<T>>;

template<typename T>
constexpr auto makeOptionalRef(T& obj)
{
  return std::make_optional<std::reference_wrapper<T>>(obj);
}
} // namespace qmcplusplus
#endif
