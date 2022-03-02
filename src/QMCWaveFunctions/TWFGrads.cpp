//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "TWFGrads.hpp"
#include <algorithm>

namespace qmcplusplus
{
TWFGrads<CoordsType::POS>::TWFGrads(const std::size_t size) : grads_positions(size) {}

TWFGrads<CoordsType::POS>& TWFGrads<CoordsType::POS>::operator+=(const TWFGrads<CoordsType::POS>& rhs)
{
  assert(grads_positions.size() == rhs.grads_positions.size());
  std::transform(grads_positions.begin(), grads_positions.end(), rhs.grads_positions.begin(), grads_positions.begin(),
                 [](const QMCTraits::GradType& x, const QMCTraits::GradType& y) { return x + y; });
  return *this;
}

TWFGrads<CoordsType::POS_SPIN>::TWFGrads(const std::size_t size) : grads_positions(size), grads_spins(size) {}

TWFGrads<CoordsType::POS_SPIN>& TWFGrads<CoordsType::POS_SPIN>::operator+=(const TWFGrads<CoordsType::POS_SPIN>& rhs)
{
  assert(grads_positions.size() == rhs.grads_positions.size());
  std::transform(grads_positions.begin(), grads_positions.end(), rhs.grads_positions.begin(), grads_positions.begin(),
                 [](const QMCTraits::GradType& x, const QMCTraits::GradType& y) { return x + y; });
  std::transform(grads_spins.begin(), grads_spins.end(), rhs.grads_spins.begin(), grads_spins.begin(),
                 [](const QMCTraits::ComplexType& x, const QMCTraits::ComplexType& y) { return x + y; });
  return *this;
}

template struct TWFGrads<CoordsType::POS>;
template struct TWFGrads<CoordsType::POS_SPIN>;
} // namespace qmcplusplus
