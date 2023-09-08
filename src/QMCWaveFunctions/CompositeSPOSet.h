//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_COMPOSITE_SPOSET_H
#define QMCPLUSPLUS_COMPOSITE_SPOSET_H

#include "Configuration.h"
#include "QMCWaveFunctions/CompositeSPOSetT.h"

namespace qmcplusplus
{
using CompositeSPOSet        = CompositeSPOSetT<QMCTraits::ValueType>;
using CompositeSPOSetBuilder = CompositeSPOSetBuilderT<QMCTraits::ValueType>;
} // namespace qmcplusplus

#endif
