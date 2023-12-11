//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file Timer.h
 * @brief Timer class
 */
#ifndef QMCPLUSPLUS_TIMER_H
#define QMCPLUSPLUS_TIMER_H

#include "Utilities/Clock.h"

namespace qmcplusplus
{
struct Timer
{
  using Clock = std::chrono::system_clock;
  Clock::time_point start_time;
  inline Timer() { start_time = Clock::now(); }
  inline void restart() { start_time = Clock::now(); }
  inline double elapsed() const
  {
    std::chrono::duration<double> elapsed = Clock::now() - start_time;
    return elapsed.count();
  }
};
} // namespace qmcplusplus
#endif
