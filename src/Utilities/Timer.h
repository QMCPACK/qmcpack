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
  double start_time;
  inline Timer() { start_time = CPUClock()(); }
  inline void restart() { start_time = CPUClock()(); }
  inline double elapsed() const { return CPUClock()() - start_time; }
};
} // namespace qmcplusplus
#endif
