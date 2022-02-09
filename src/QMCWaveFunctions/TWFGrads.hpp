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


#ifndef QMCPLUSPLUS_TWFGRADS_HPP
#define QMCPLUSPLUS_TWFGRADS_HPP

#include "Configuration.h"
#include "Particle/MCCoords.hpp"

namespace qmcplusplus
{
template<CoordsType CT = CoordsType::POS>
struct TWFGrads
{
  void resize(const std::size_t size) { grads_positions.resize(size); }
  std::vector<QMCTraits::GradType> grads_positions;
};

template<>
struct TWFGrads<CoordsType::POS_SPIN>
{
  void resize(const std::size_t size)
  {
    grads_positions.resize(size);
    grads_spins.resize(size);
  }
  std::vector<QMCTraits::GradType> grads_positions;
  std::vector<QMCTraits::ComplexType> grads_spins;
};

extern template struct TWFGrads<CoordsType::POS>;
extern template struct TWFGrads<CoordsType::POS_SPIN>;
} // namespace qmcplusplus

#endif
