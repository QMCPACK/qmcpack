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

void MCCoords<CoordsType::POS>::getSubset(const std::size_t offset,
                                          const std::size_t size,
                                          MCCoords<CoordsType::POS>& out)
{
  assert(out.positions.size() == size);
  auto start = positions.begin() + offset * size;
  auto end   = start + size;
  std::copy(start, end, out.positions.begin());
}

void MCCoords<CoordsType::POS_SPIN>::resize(const std::size_t size)
{
  positions.resize(size);
  spins.resize(size);
}

void MCCoords<CoordsType::POS_SPIN>::getSubset(const std::size_t offset,
                                               const std::size_t size,
                                               MCCoords<CoordsType::POS_SPIN>& out)
{
  assert(out.positions.size() == size);
  assert(out.spins.size() == size);
  auto pos_start  = positions.begin() + offset * size;
  auto pos_end    = pos_start + size;
  auto spin_start = spins.begin() + offset * size;
  auto spin_end   = spin_start + size;
  std::copy(pos_start, pos_end, out.positions.begin());
  std::copy(spin_start, spin_end, out.spins.begin());
}

template struct MCCoords<CoordsType::POS>;
template struct MCCoords<CoordsType::POS_SPIN>;
} // namespace qmcplusplus
