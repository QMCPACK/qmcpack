//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "MCCoords.hpp"

namespace qmcplusplus
{
template<MCCoordsTypes MCT>
void MCCoords<MCT>::resize(const std::size_t size)
{
  rs.resize(size);
}

// template<>
// void MCCoords<MCCoordsTypes::RSSPINS>::resize(const std::size_t size)
// {
//   rs.resize(size);
//   spins.resize(size);
// }

template struct MCCoords<MCCoordsTypes::RS>;
template struct MCCoords<MCCoordsTypes::RSSPINS>;
}
