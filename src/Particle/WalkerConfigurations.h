//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file WalkerConfigurations.h
 * @brief Declaration of a WalkerConfigurations
 */
#ifndef QMCPLUSPLUS_WALKERCONFIGURATIONS_H
#define QMCPLUSPLUS_WALKERCONFIGURATIONS_H
#include "Configuration.h"
#include "Particle/WalkerConfigurationsT.h"

namespace qmcplusplus
{
using WalkerConfigurations = WalkerConfigurationsT<QMCTraits::ValueType>;
} // namespace qmcplusplus
#endif
