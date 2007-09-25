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

#include "Message/OpenMP.h"

#if defined(ENABLE_OPENMP)
namespace qmcplusplus  {
  /** Timer using omp_get_wtime 
   */
  struct Timer
  {
    double start_time;
    inline Timer() { start_time=omp_get_wtime();}
    inline void restart() {start_time=omp_get_wtime();}
    inline double elapsed() const {
      return omp_get_wtime()-start_time;
    }
  };
}
#else /* use boost or pooma */
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
#endif
