//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file Timer.h
 * @brief Timer class using boost::timer
 */
#ifndef QMCPLUSPLUS_TIMER_H
#define QMCPLUSPLUS_TIMER_H

#include <Utilities/Clock.h>

namespace qmcplusplus
{
struct Timer
{
  double start_time;
  inline Timer()
  {
    start_time=cpu_clock();
  }
  inline void restart()
  {
    start_time=cpu_clock();
  }
  inline double elapsed() const
  {
    return cpu_clock()-start_time;
  }
};
}
#endif
