//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_DMCREFENERGYSCHEME_H
#define QMCPLUSPLUS_DMCREFENERGYSCHEME_H

namespace qmcplusplus
{
/** DMCRefEnergy schemes
 */
enum class DMCRefEnergyScheme
{
  UNLIMITED_HISTORY,
  LIMITED_HISTORY
};
} // namespace qmcplusplus
#endif
