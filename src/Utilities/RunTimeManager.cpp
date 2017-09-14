//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file RunTimeManager.cpp
 *  @brief Class for determining elapsed run time enabling simulations to adjust to time limits.

 */
#include <Utilities/RunTimeManager.h>


namespace qmcplusplus
{

RunTimeManagerClass RunTimeManager;

double
LoopTimer::get_time_per_iteration() {
  if (nloop > 0)
  {
    return total_time/nloop;
  }
  return 0.0;
}

}
