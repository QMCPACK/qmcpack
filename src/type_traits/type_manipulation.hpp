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

#ifndef QMCPLUSPLUS_TYPE_MANIPULATIONS_H
#define QMCPLUSPLUS_TYPE_MANIPULATIONS_H

#include <type_traits>

/** Type with the bottom const removed in types of X (const) * const *
 */
template<typename CT>
struct BottomConstRemoved {
  using type = typename std::add_pointer<typename std::remove_const<typename std::remove_pointer<CT>::type>::type>::type;
};

#endif
