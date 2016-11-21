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

void TimerManagerClass::get_stack_roots(std::vector<NewTimer *> &roots)
{
  for (int i = 0; i < TimerList.size(); i++)
  {
    if (TimerList[i]->get_parent() == NULL)
    {
      roots.push_back(TimerList[i]);
    }
  }
}

void TimerManagerClass::collate_flat_profile(Communicate *comm, FlatProfileData &p)
{
  for(int i=0; i<TimerList.size(); ++i)
  {
    NewTimer &timer = *TimerList[i];
    nameList_t::iterator it(p.nameList.find(timer.get_name()));
    if(it == p.nameList.end())
    {
      int ind=p.nameList.size();
      p.nameList[timer.get_name()]=ind;
      p.timeList.push_back(timer.get_total());
      p.callList.push_back(timer.get_num_calls());
    }
    else
    {
      int ind=(*it).second;
      p.timeList[ind]+=timer.get_total();
      p.callList[ind]+=timer.get_num_calls();
    }
  }

  if (comm)
  {
    comm->allreduce(p.timeList);
    comm->allreduce(p.callList);
  }
}

void TimerManagerClass::collate_stack_profile(Communicate *comm, StackProfileData &p)
{
  for(int i=0; i<TimerList.size(); ++i)
  {
    NewTimer &timer = *TimerList[i];

    std::map<std::string,double>::iterator stack_name_it = timer.get_per_stack_total_time().begin();
    for (; stack_name_it != timer.get_per_stack_total_time().end(); stack_name_it++)
    {
      const std::string &stack_name = stack_name_it->first;
      nameList_t::iterator it(p.nameList.find(stack_name));

      double totalTime = timer.get_total(stack_name);
      double exclTime = timer.get_exclusive_time(stack_name);
      int ncalls = timer.get_num_calls(stack_name);

      if(it == p.nameList.end())
      {
        int ind=p.nameList.size();
        p.nameList[stack_name]=ind;
        p.timeList.push_back(totalTime);
        p.timeExclList.push_back(exclTime);
        p.callList.push_back(ncalls);
      }
      else
      {
        int ind=(*it).second;
        p.timeList[ind]+=totalTime;
        p.timeExclList[ind]+=exclTime;
        p.callList[ind]+=timer.get_num_calls();
      }
    }
  }

  if (comm)
  {
    comm->allreduce(p.timeList);
    comm->allreduce(p.timeExclList);
    comm->allreduce(p.callList);
  }
}

void
TimerManagerClass::print(Communicate* comm)
{
  printf("Stack timer profile\n");
  print_stack(comm);
  printf("\nFlat profile\n");
  print_flat(comm);

}

void
TimerManagerClass::print_flat(Communicate* comm)
{
#if ENABLE_TIMER
  FlatProfileData p;

  collate_flat_profile(comm, p);

  if(comm->rank() == 0)
  {
    #pragma omp master
    {
      std::map<std::string,int>::iterator it(p.nameList.begin()), it_end(p.nameList.end());
      while(it != it_end)
      {
        int i=(*it).second;
        //if(callList[i]) //skip zeros
        printf ("%-40s  %9.4f  %13ld  %16.9f  %12.6f TIMER\n"
        , (*it).first.c_str()
        , p.timeList[i], p.callList[i]
        , p.timeList[i]/(static_cast<double>(p.callList[i])+std::numeric_limits<double>::epsilon())
        , p.timeList[i]/static_cast<double>(omp_get_max_threads()*comm->size()));
        ++it;
      }
    }
  }
#endif
}

void
TimerManagerClass::print_stack(Communicate* comm)
{
#if ENABLE_TIMER
  StackProfileData p;

  collate_stack_profile(comm, p);
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
