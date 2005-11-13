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
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
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

#if defined(HAVE_LIBBOOST)
#include <boost/timer.hpp>
namespace qmcplusplus {
  typedef boost::timer Timer;
}
#else
#include "Utilities/Clock.h"
namespace qmcplusplus {

  struct Timer: private Pooma::Clock {

    Timer() { }
    inline void restart() { start();}
    inline double elapsed() { 
      stop();
      return cpu_time();
    }
  };
}
#endif
#endif
