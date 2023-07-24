//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SOA_LCAO_SPINOR_BUILDER_H
#define QMCPLUSPLUS_SOA_LCAO_SPINOR_BUILDER_H

#include "Configuration.h"
#include "QMCWaveFunctions/LCAO/LCAOSpinorBuilderT.h"

namespace qmcplusplus
{
using LCAOSpinorBuilder = LCAOSpinorBuilderT<QMCTraits::ValueType>;

} // namespace qmcplusplus
#endif
