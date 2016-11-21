//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


/** @file NewTimer.cpp
 * @brief Implements TimerManager
 */
#include "Utilities/NewTimer.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include <map>
#include <limits>
#include <cstdio>
namespace qmcplusplus
{
TimerManagerClass TimerManager;

void TimerManagerClass::addTimer(NewTimer* t)
{
  #pragma omp critical
  {
    t->set_manager(this);
    TimerList.push_back(t);
  }
}

void TimerManagerClass::reset()
{
  for (int i=0; i<TimerList.size(); i++)
    TimerList[i]->reset();
}

void TimerManagerClass::collate_flat_profile(Communicate *comm,
                                        std::map<std::string, int> &nameList,
                                        std::vector<double> &timeList,
                                        std::vector<long> &callList)
{
  for(int i=0; i<TimerList.size(); ++i)
  {
    NewTimer &timer = *TimerList[i];
    nameList_t::iterator it(nameList.find(timer.get_name()));
    if(it == nameList.end())
    {
      int ind=nameList.size();
      nameList[timer.get_name()]=ind;
      timeList.push_back(timer.get_total());
      callList.push_back(timer.get_num_calls());
    }
    else
    {
      int ind=(*it).second;
      timeList[ind]+=timer.get_total();
      callList[ind]+=timer.get_num_calls();
    }
  }

  if (comm)
  {
    comm->allreduce(timeList);
    comm->allreduce(callList);
  }
}

void TimerManagerClass::collate_stack_profile(Communicate *comm,
                                        std::map<std::string, int> &nameList,
                                        std::vector<double> &timeList,
                                        std::vector<double> &timeExclList,
                                        std::vector<long> &callList)
{
  for(int i=0; i<TimerList.size(); ++i)
  {
    NewTimer &timer = *TimerList[i];

    std::map<std::string,double>::iterator stack_name_it = timer.get_per_stack_total_time().begin();
    for (; stack_name_it != timer.get_per_stack_total_time().end(); stack_name_it++)
    {
      const std::string &stack_name = stack_name_it->first;
      nameList_t::iterator it(nameList.find(stack_name));

      double totalTime = timer.get_total(stack_name);
      double exclTime = timer.get_exclusive_time(stack_name);
      int ncalls = timer.get_num_calls(stack_name);

      if(it == nameList.end())
      {
        int ind=nameList.size();
        nameList[stack_name]=ind;
        timeList.push_back(totalTime);
        timeExclList.push_back(exclTime);
        callList.push_back(ncalls);
      }
      else
      {
        int ind=(*it).second;
        timeList[ind]+=totalTime;
        timeExclList[ind]+=exclTime;
        callList[ind]+=timer.get_num_calls();
      }
    }
  }

  if (comm)
  {
    comm->allreduce(timeList);
    comm->allreduce(timeExclList);
    comm->allreduce(callList);
  }
}

void
TimerManagerClass::print(Communicate* comm)
{
#if ENABLE_TIMER
  std::map<std::string,int> nameList;
  std::vector<double> timeList;
  std::vector<long>   callList;

  collate_flat_profile(comm, nameList, timeList, callList);

  if(comm->rank() == 0)
  {
    #pragma omp master
    {
      std::map<std::string,int>::iterator it(nameList.begin()), it_end(nameList.end());
      while(it != it_end)
      {
        int i=(*it).second;
        //if(callList[i]) //skip zeros
        fprintf (stderr, "%-40s  %9.4f  %13ld  %16.9f  %12.6f TIMER\n"
        , (*it).first.c_str()
        , timeList[i], callList[i]
        , timeList[i]/(static_cast<double>(callList[i])+std::numeric_limits<double>::epsilon())
        , timeList[i]/static_cast<double>(omp_get_max_threads()*comm->size()));
        ++it;
      }
    }
  }
#endif
}


std::string
NewTimer::get_stack_name()
{
  std::string stack_name = name;
  NewTimer *current = parent;
  while (current)
  {
    stack_name += "/";
    stack_name += current->get_name();
    current = current->parent;
  }
  return stack_name;
}

double
NewTimer::get_exclusive_time(const std::string &stack_name)
{
  double exclusive_total = per_stack_total_time[stack_name];
  for (int i = 0; i < children.size(); i++)
  {
    std::string tmp_stack = children[i]->get_name();
    tmp_stack += "/";
    tmp_stack += stack_name;

    exclusive_total -= children[i]->get_total(tmp_stack);
  }
  return exclusive_total;
}

}
