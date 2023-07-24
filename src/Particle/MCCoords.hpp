//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//                    Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MCCOORDS_HPP
#define QMCPLUSPLUS_MCCOORDS_HPP

#include "Configuration.h"
#include "MCCoordsT.hpp"

namespace qmcplusplus
{
template<CoordsType MCT>
using MCCoords = MCCoordsT<QMCTraits::ValueType, MCT>;

} // namespace qmcplusplus

#endif
