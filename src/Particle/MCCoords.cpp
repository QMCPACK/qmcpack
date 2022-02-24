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
MCCoords<CoordsType::POS> MCCoords<CoordsType::POS>::getSubset(const std::size_t offset, const std::size_t size)
{
  MCCoords<CoordsType::POS> out(size);
  auto start = positions.begin() + offset * size;
  auto end   = start + size;
  std::copy(start, end, out.positions.begin());
  return out;
}

MCCoords<CoordsType::POS_SPIN> MCCoords<CoordsType::POS_SPIN>::getSubset(const std::size_t offset,
                                                                         const std::size_t size)
{
  MCCoords<CoordsType::POS_SPIN> out(size);
  auto pos_start  = positions.begin() + offset * size;
  auto pos_end    = pos_start + size;
  auto spin_start = spins.begin() + offset * size;
  auto spin_end   = spin_start + size;
  std::copy(pos_start, pos_end, out.positions.begin());
  std::copy(spin_start, spin_end, out.spins.begin());
  return out;
}

template struct MCCoords<CoordsType::POS>;
template struct MCCoords<CoordsType::POS_SPIN>;
} // namespace qmcplusplus
