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

namespace qmcplusplus
{

template<CoordsType CT>
void TWFGrads<CT>::resize(const std::size_t size)
{
  grads_positions.resize(size);
}

template struct TWFGrads<CoordsType::POS>;
template struct TWFGrads<CoordsType::POS_SPIN>;
}
