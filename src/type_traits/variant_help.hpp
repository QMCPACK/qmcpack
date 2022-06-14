//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_VARIANT_HELP_HPP
#define QMCPLUSPLUS_VARIANT_HELP_HPP

/** \file
 *  A collection of functions to help with the use of variants
 */

namespace qmcplusplus
{

template<class T, class R>
constexpr bool has(const R& this_one)
{
  return std::holds_alternative<T>(this_one);
}

} // namespace qmcplusplus

#endif
