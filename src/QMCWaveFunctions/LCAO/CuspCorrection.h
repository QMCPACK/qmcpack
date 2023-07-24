//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file CuspCorrection.h
  * @brief Corrections to electron-nucleus cusp for all-electron molecular calculations.
  */

#ifndef QMCPLUSPLUS_CUSPCORRECTION_H
#define QMCPLUSPLUS_CUSPCORRECTION_H

#include "Configuration.h"
#include "QMCWaveFunctions/LCAO/CuspCorrectionT.h"

namespace qmcplusplus
{
using CuspCorrectionParameters = CuspCorrectionParametersT<QMCTraits::ValueType>;

using CuspCorrection = CuspCorrectionT<QMCTraits::ValueType>;

} // namespace qmcplusplus

#endif
