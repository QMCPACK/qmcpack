//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_DISTANCETABLEDATAIMPL_H
#define QMCPLUSPLUS_DISTANCETABLEDATAIMPL_H

#include "Configuration.h"
#include "Particle/DistanceTableT.h"

namespace qmcplusplus
{
using DistanceTable = DistanceTableT<QMCTraits::ValueType>;
using DistanceTableAA = DistanceTableAAT<QMCTraits::ValueType>;
using DistanceTableAB = DistanceTableABT<QMCTraits::ValueType>;
} // namespace qmcplusplus
#endif
