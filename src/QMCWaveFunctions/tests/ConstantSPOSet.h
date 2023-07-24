//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 Raymond Clay and QMCPACK developers.
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_CONSTANTSPOSET_H
#define QMCPLUSPLUS_CONSTANTSPOSET_H

#include "Configuration.h"
#include "QMCWaveFunctions/tests/ConstantSPOSetT.h"

namespace qmcplusplus
{
using ConstantSPOSet = ConstantSPOSetT<QMCTraits::ValueType>;

} // namespace qmcplusplus
#endif
