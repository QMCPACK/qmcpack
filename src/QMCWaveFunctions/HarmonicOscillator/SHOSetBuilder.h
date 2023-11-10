//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SHO_BASIS_BUILDER_H
#define QMCPLUSPLUS_SHO_BASIS_BUILDER_H

#include "Configuration.h"
#include "QMCWaveFunctions/HarmonicOscillator/SHOSetBuilderT.h"

namespace qmcplusplus
{
using SHOSetBuilder = SHOSetBuilderT<QMCTraits::ValueType>;
} // namespace qmcplusplus

#endif
