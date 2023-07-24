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

#ifndef QMCPLUSPLUS_MCCOORDST_HPP
#define QMCPLUSPLUS_MCCOORDST_HPP

#include "ParticleSetTraits.h"
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

template<typename T, CoordsType MCT>
struct MCCoordsT;

template<typename T>
struct MCCoordsT<T, CoordsType::POS>
{
  using PosType = typename ParticleSetTraits<T>::PosType;

  MCCoordsT(const std::size_t size) : positions(size) {}

  MCCoordsT& operator+=(const MCCoordsT& rhs);

  /** get subset of MCCoordsT
     * [param,out] out
     */
  void getSubset(const std::size_t offset, const std::size_t size, MCCoordsT<T, CoordsType::POS>& out) const;

  std::vector<PosType> positions;
};

template<typename T>
struct MCCoordsT<T, CoordsType::POS_SPIN>
{
  using PosType          = typename ParticleSetTraits<T>::PosType;
  using FullPrecRealType = typename ParticleSetTraits<T>::FullPrecRealType;

  MCCoordsT(const std::size_t size) : positions(size), spins(size) {}

  MCCoordsT& operator+=(const MCCoordsT& rhs);

  /** get subset of MCCoordsT
     * [param,out] out
     */
  void getSubset(const std::size_t offset, const std::size_t size, MCCoordsT<T, CoordsType::POS_SPIN>& out) const;

  std::vector<PosType> positions;
  std::vector<FullPrecRealType> spins;
};
} // namespace qmcplusplus

#endif
