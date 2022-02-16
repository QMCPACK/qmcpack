//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "MCCoords.hpp"

namespace qmcplusplus
{
void MCCoords<CoordsType::POS>::resize(const std::size_t size) { positions.resize(size); }

void MCCoords<CoordsType::POS_SPIN>::resize(const std::size_t size)
{
  positions.resize(size);
  spins.resize(size);
}
template struct MCCoords<CoordsType::POS>;
template struct MCCoords<CoordsType::POS_SPIN>;
} // namespace qmcplusplus
