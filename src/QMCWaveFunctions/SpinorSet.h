//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers
//
// File developed by: Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//                    Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by:  Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SPINORSET_H
#define QMCPLUSPLUS_SPINORSET_H

#include "Configuration.h"
#include "QMCWaveFunctions/SpinorSetT.h"

namespace qmcplusplus
{
using SpinorSet = SpinorSetT<QMCTraits::ValueType>;

} // namespace qmcplusplus
#endif
