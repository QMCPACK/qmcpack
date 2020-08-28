//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "Clock.h"

namespace qmcplusplus
{
double fake_cpu_clock::fake_cpu_clock_value     = 0.0;
double fake_cpu_clock::fake_cpu_clock_increment = 1.0;

} // namespace qmcplusplus
