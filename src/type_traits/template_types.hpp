//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_TEMPLATE_TYPES_HPP
#define QMCPLUSPLUS_TEMPLATE_TYPES_HPP

#include <vector>
#include <functional>

namespace qmcplusplus
{

template<typename T>
using RefVector = std::vector<std::reference_wrapper<T>>;

}
#endif
