//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@inte.com, Intel Corp.
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file
 *
 * The most general reader class for the following classes using the full single grid for the supercell
 * - SplineR2R
 * - SplineC2C
 * - SplineC2R
 * Each band is initialized with UBspline_3d_d and both real and imaginary parts are passed to the objects
 * which will convert the data type to their internal precision. 
 */
#ifndef QMCPLUSPLUS_SPLINESET_READER_H
#define QMCPLUSPLUS_SPLINESET_READER_H

#include "Configuration.h"
#include "QMCWaveFunctions/BsplineFactory/SplineSetReaderT.h"

namespace qmcplusplus
{
template<typename SPLINEBASE>
using SplineSetReader = SplineSetReaderT<SPLINEBASE>;

} // namespace qmcplusplus
#endif
