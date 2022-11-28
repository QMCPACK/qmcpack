//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_CLOCK_H
#define QMCPLUSPLUS_CLOCK_H

#include <stddef.h>
#include <chrono>

namespace qmcplusplus
{

using ChronoClock = std::chrono::system_clock;

// Implements a std::chrono clock
// See https://github.com/korfuri/fake_clock

class FakeChronoClock
{
public:
  using duration   = std::chrono::nanoseconds;
  using rep        = duration::rep;
  using period     = duration::period;
  using time_point = std::chrono::time_point<FakeChronoClock>;
  static time_point now() noexcept
  {
    fake_chrono_clock_value += fake_chrono_clock_increment;
    return fake_chrono_clock_value;
  }

  static time_point fake_chrono_clock_value;
  static duration fake_chrono_clock_increment;
};

} // namespace qmcplusplus
#endif
