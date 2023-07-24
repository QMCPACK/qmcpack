//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file SPOSetBuilder.h
 * @brief Declaration of a base class of SPOSet Builders
 */
#ifndef QMCPLUSPLUS_SPOSET_BUILDER_H
#define QMCPLUSPLUS_SPOSET_BUILDER_H

#include "Configuration.h"
#include "QMCWaveFunctions/SPOSetBuilderT.h"

namespace qmcplusplus
{
using SPOSetBuilder = SPOSetBuilderT<QMCTraits::ValueType>;

} // namespace qmcplusplus
#endif
