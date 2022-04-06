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

MCCoords<CoordsType::POS>& MCCoords<CoordsType::POS>::operator+=(const MCCoords<CoordsType::POS>& rhs)
{
  assert(positions.size() == rhs.positions.size());
  std::transform(positions.begin(), positions.end(), rhs.positions.begin(), positions.begin(),
                 [](const QMCTraits::PosType& x, const QMCTraits::PosType& y) { return x + y; });
  return *this;
}

void MCCoords<CoordsType::POS_SPIN>::getSubset(const std::size_t offset,
                                               const std::size_t size,
                                               MCCoords<CoordsType::POS_SPIN>& out) const
{
  std::copy_n(positions.begin() + offset, size, out.positions.begin());
  std::copy_n(spins.begin() + offset, size, out.spins.begin());
}

MCCoords<CoordsType::POS_SPIN>& MCCoords<CoordsType::POS_SPIN>::operator+=(const MCCoords<CoordsType::POS_SPIN>& rhs)
{
  assert(positions.size() == rhs.positions.size());
  std::transform(positions.begin(), positions.end(), rhs.positions.begin(), positions.begin(),
                 [](const QMCTraits::PosType& x, const QMCTraits::PosType& y) { return x + y; });
  std::transform(spins.begin(), spins.end(), rhs.spins.begin(), spins.begin(),
                 [](const QMCTraits::FullPrecRealType& x, const QMCTraits::FullPrecRealType& y) { return x + y; });
  return *this;
}

template struct MCCoords<CoordsType::POS>;
template struct MCCoords<CoordsType::POS_SPIN>;
} // namespace qmcplusplus
