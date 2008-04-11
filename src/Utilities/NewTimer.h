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
#ifndef QMCPLUSPLUS_NEW_TIMER_H
#define QMCPLUSPLUS_NEW_TIMER_H

#include "Message/OpenMP.h"
#include <vector>
#include <string>
#include <algorithm>

namespace qmcplusplus  {

#if defined(ENABLE_OPENMP)
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
    
    inline const std::string& get_name() const 
    { return name;                   }

    inline void reset()           
    { num_calls = 0; total_time=0.0; }
        
    NewTimer(std::string myname) : 
      total_time(0.0), num_calls(0.0), name(myname)
    { }
  };
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
      start_time = tv.tv_usec; 
    }
    inline void stop()  
    { 
      gettimeofday(&tv, NULL);
      total_time += 1.0e-6*(double)(tv.tv_usec - start_time);
      num_calls++;
    }
    inline double    get_total() const 
    { return total_time; }

    inline long  get_num_calls() const 
    { return num_calls;  }

    inline std::string get_name() const
    { return name; }

    inline void reset()          { num_calls = 0; total_time=0.0; }

    NewTimer(string myname) : 
      total_time(0.0), num_calls(0.0), name(myname)
    { }
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
    void addTimer (NewTimer* timer)
    {
      TimerList.push_back(timer);
    }

    void print () 
    {
      TimerComparator comp;
      std::sort(TimerList.begin(), TimerList.end(), comp);
      std::vector<std::string> nameList;
      std::vector<double> timeList;
      std::vector<long>   callList;
      std::string lastName = "";
      int numDistinct = 0;
      for (int i=0; i<TimerList.size(); i++) {
	NewTimer &timer = *TimerList[i];
	if (timer.get_name() == lastName && lastName != "") {
	  timeList[numDistinct-1]  += timer.get_total();
	  callList[numDistinct-1] += timer.get_num_calls();
	}
	else {
	  nameList.push_back(timer.get_name());
	  timeList.push_back(timer.get_total());
	  callList.push_back(timer.get_num_calls());
	  lastName == timer.get_name();
	  numDistinct++;
	}
      }
      fprintf (stderr, "Routine name                   Total time    Num Calls    Time per call\n");
      for (int i=0; i<numDistinct; i++) {
	fprintf (stderr, "%30s  %5.3f  %ld  %5.9f\n", nameList[i].c_str(),
		 timeList[i], callList[i], timeList[i]/(double)callList[i]);
      }
    }
  };

  extern TimerManagerClass TimerManager;
}


#endif
