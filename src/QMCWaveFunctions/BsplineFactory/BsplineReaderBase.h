//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file BsplineReaderBase.h
 *
 * base class to read data and manage spline tables
 */
#ifndef QMCPLUSPLUS_BSPLINE_READER_BASE_H
#define QMCPLUSPLUS_BSPLINE_READER_BASE_H

#include "Configuration.h"
#include "QMCWaveFunctions/BsplineFactory/BsplineReaderBaseT.h"

namespace qmcplusplus
{
using BsplineReaderBase = BsplineReaderBaseT<QMCTraits::ValueType>;

} // namespace qmcplusplus
#endif
