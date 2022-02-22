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
#include "type_traits/complex_help.hpp"
#include <algorithm>

#include <vector>

namespace qmcplusplus
{
enum class CoordsType
{
  POS,
  POS_SPIN
};

template<CoordsType MCT>
struct MCCoords;

template<>
struct MCCoords<CoordsType::POS>
{
  void resize(const std::size_t size);
  /** get subset of MCCoords
   * [param,out] out
   */
  void getSubset(const std::size_t offset, const std::size_t size, MCCoords<CoordsType::POS>& out);

  std::vector<QMCTraits::PosType> positions;
};

template<>
struct MCCoords<CoordsType::POS_SPIN>
{
  void resize(const std::size_t size);
  /** get subset of MCCoords
   * [param,out] out
   */
  void getSubset(const std::size_t offset, const std::size_t size, MCCoords<CoordsType::POS_SPIN>& out);

  std::vector<QMCTraits::PosType> positions;
  std::vector<QMCTraits::FullPrecRealType> spins;
};

extern template struct MCCoords<CoordsType::POS>;
extern template struct MCCoords<CoordsType::POS_SPIN>;

template<CoordsType CT>
MCCoords<CT> operator+(const MCCoords<CT>& lhs, const MCCoords<CT>& rhs)
{
  MCCoords<CT> out;
  assert(lhs.positions.size() == rhs.positions.size());
  out.resize(lhs.positions.size());
  std::transform(lhs.positions.begin(), lhs.positions.end(), rhs.positions.begin(), out.positions.begin(),
                 [](const QMCTraits::PosType& x, const QMCTraits::PosType& y) { return x + y; });
  if constexpr (CT == CoordsType::POS_SPIN)
    std::transform(lhs.spins.begin(), lhs.spins.end(), rhs.spins.begin(), out.spins.begin(),
                   [](const QMCTraits::FullPrecRealType& x, const QMCTraits::FullPrecRealType& y) { return x + y; });
  return out;
}

template<CoordsType CT>
MCCoords<CT> operator-(const MCCoords<CT>& lhs, const MCCoords<CT>& rhs)
{
  MCCoords<CT> out;
  assert(lhs.positions.size() == rhs.positions.size());
  out.resize(lhs.positions.size());
  std::transform(lhs.positions.begin(), lhs.positions.end(), rhs.positions.begin(), out.positions.begin(),
                 [](const QMCTraits::PosType& x, const QMCTraits::PosType& y) { return x - y; });
  if constexpr (CT == CoordsType::POS_SPIN)
    std::transform(lhs.spins.begin(), lhs.spins.end(), rhs.spins.begin(), out.spins.begin(),
                   [](const QMCTraits::FullPrecRealType& x, const QMCTraits::FullPrecRealType& y) { return x - y; });
  return out;
}

} // namespace qmcplusplus

#endif
