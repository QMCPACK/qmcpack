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
  MCCoords(const std::size_t size) : positions(size) {}

  MCCoords& operator+=(const MCCoords& rhs);

  /** get subset of MCCoords
   * [param,out] out
   */
  void getSubset(const std::size_t offset, const std::size_t size, MCCoords<CoordsType::POS>& out) const;

  std::vector<QMCTraits::PosType> positions;
};

template<>
struct MCCoords<CoordsType::POS_SPIN>
{
  MCCoords(const std::size_t size) : positions(size), spins(size) {}

  MCCoords& operator+=(const MCCoords& rhs);

  /** get subset of MCCoords
   * [param,out] out
   */
  void getSubset(const std::size_t offset, const std::size_t size, MCCoords<CoordsType::POS_SPIN>& out) const;

  std::vector<QMCTraits::PosType> positions;
  std::vector<QMCTraits::FullPrecRealType> spins;
};

extern template struct MCCoords<CoordsType::POS>;
extern template struct MCCoords<CoordsType::POS_SPIN>;
} // namespace qmcplusplus

#endif
