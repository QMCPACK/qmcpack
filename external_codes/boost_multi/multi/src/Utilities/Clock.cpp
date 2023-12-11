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

FakeChronoClock::time_point FakeChronoClock::fake_chrono_clock_value   = FakeChronoClock::time_point();
FakeChronoClock::duration FakeChronoClock::fake_chrono_clock_increment = std::chrono::seconds(1);

} // namespace qmcplusplus
