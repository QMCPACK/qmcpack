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
void MCCoords<CoordsType::POS>::getSubset(const std::size_t offset,
                                          const std::size_t size,
                                          MCCoords<CoordsType::POS>& out) const
{
  std::copy_n(positions.begin() + offset, size, out.positions.begin());
}

void MCCoords<CoordsType::POS_SPIN>::getSubset(const std::size_t offset,
                                               const std::size_t size,
                                               MCCoords<CoordsType::POS_SPIN>& out) const
{
  std::copy_n(positions.begin() + offset, size, out.positions.begin());
  std::copy_n(spins.begin() + offset, size, out.spins.begin());
}

template struct MCCoords<CoordsType::POS>;
template struct MCCoords<CoordsType::POS_SPIN>;
} // namespace qmcplusplus
