//////////////////////////////////////////////////////////////////
// (c) Copyright 2008- by Ken Esler
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Ken Esler
//   National Center for Supercomputing Applications &
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file NewTimer.h
 * @brief NewTimer class various high-resolution timers.
 */
#ifndef QMCPLUSPLUS_NEW_TIMER_H
#define QMCPLUSPLUS_NEW_TIMER_H

#include "Message/OpenMP.h"
#include <vector>
#include <string>
#include <algorithm>

class Communicate;

namespace qmcplusplus  {

#if defined(ENABLE_OPENMP)
#if defined(BGP_BUG)
  class NewTimer
  {
  protected:
    double start_time;
    double total_time;
    long num_calls;
    std::string name;
  public:
    inline void start() 
    { }
    
    inline void stop()  
    { }

    inline double    get_total() const  
    { return total_time;             }
    
    inline long  get_num_calls() const  
    { return num_calls;              }
    
    inline std::string get_name() const 
    { return name;                   }

    inline void reset()           
    { num_calls = 0; total_time=0.0; }
        
    NewTimer(const std::string& myname) : 
      total_time(0.0), num_calls(0), name(myname)
    { }

    void set_name(const std::string& myname)
    {
      name=myname;
    }
  };
#else
  /* Timer using omp_get_wtime  */
  class NewTimer
  {
  protected:
    double start_time;
    double total_time;
    long num_calls;
    std::string name;
  public:
    inline void start() 
    { start_time = omp_get_wtime(); }
    
    inline void stop()  
    { total_time += omp_get_wtime() - start_time;  num_calls++;   }

    inline double    get_total() const  
    { return total_time;             }
    
    inline long  get_num_calls() const  
    { return num_calls;              }
    
    inline std::string get_name() const 
    { return name;                   }

    inline void reset()           
    { num_calls = 0; total_time=0.0; }
        
    NewTimer(const std::string& myname) : 
      total_time(0.0), num_calls(0), name(myname)
    { }

    void set_name(const std::string& myname)
    {
      name=myname;
    }
  };
#endif
#else /* use boost or pooma */
#include <sys/time.h>
  /* Timer using gettimeofday  */
  class NewTimer
  {
  protected:
    suseconds_t start_time;
    double total_time;
    long num_calls;
    struct timeval tv;
    std::string name;
  public:
    inline void start() 
    { 
      gettimeofday(&tv, NULL);
       start_time=(double)tv.tv_sec+(1.e-6)*tv.tv_usec;
    }
    inline void stop()  
    { 
      gettimeofday(&tv, NULL);
      total_time += (double)tv.tv_sec+(1.e-6)*tv.tv_usec-start_time;
      num_calls++;
    }
    inline double    get_total() const 
    { return total_time; }

    inline long  get_num_calls() const 
    { return num_calls;  }

    inline std::string get_name() const
    { return name; }

    inline void reset()          { num_calls = 0; total_time=0.0; }

    NewTimer(const std::string& myname) : 
      total_time(0.0), num_calls(0), name(myname)
    { }

    void set_name(const std::string& myname)
    {
      name=myname;
    }
  };
#endif

  struct TimerComparator
  {
    inline bool operator()(const NewTimer *a, const NewTimer *b)
    {
      return a->get_name() < b->get_name();
    }
  };


  class TimerManagerClass
  {
  protected:
    std::vector<NewTimer*> TimerList;
  public:
    inline void addTimer (NewTimer* t)
    {
#pragma omp critical
      {
        TimerList.push_back(t);
      }
    }

    void reset();
    void print (Communicate* comm);
  };

  extern TimerManagerClass TimerManager;
}

#endif
