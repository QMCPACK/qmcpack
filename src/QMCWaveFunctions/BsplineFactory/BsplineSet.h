//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file BsplineSet.h
 *
 * BsplineSet is a SPOSet derived class and serves as a base class for B-spline SPO C2C/C2R/R2R implementation
 */
#ifndef QMCPLUSPLUS_BSPLINESET_H
#define QMCPLUSPLUS_BSPLINESET_H

#include "Configuration.h"
#include "QMCWaveFunctions/BsplineFactory/BsplineSetT.h"

namespace qmcplusplus
{
using BsplineSet = BsplineSetT<QMCTraits::ValueType>;

} // namespace qmcplusplus
#endif
