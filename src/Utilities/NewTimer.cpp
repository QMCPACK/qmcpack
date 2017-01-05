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
#include <libxml/xmlwriter.h>
#include "Configuration.h"
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
    if (t->get_name().find(TIMER_STACK_SEPARATOR) != std::string::npos)
    { 
      app_log() << "Warning: Timer name (" << t->get_name()
                << ") should not contain the character "
                << TIMER_STACK_SEPARATOR << std::endl;
    }
    t->set_manager(this);
    t->set_active_by_timer_threshold(timer_threshold);
    TimerList.push_back(t);
  }
}

void TimerManagerClass::reset()
{
  for (int i=0; i<TimerList.size(); i++)
    TimerList[i]->reset();
}

void TimerManagerClass::set_timer_threshold(const timer_levels threshold)
{
  timer_threshold = threshold;
  for (int i=0; i<TimerList.size(); i++)
  {
    TimerList[i]->set_active_by_timer_threshold(timer_threshold);
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

struct ProfileData
{
    double time;
    double calls;

    ProfileData &operator+=(const ProfileData &pd)
    {
        time += pd.time;
        calls += pd.calls;
        return *this;
    }
};

int
get_level(const std::string &stack_name)
{
  int level = 0;
  for (int i = 0; i < stack_name.length(); i++)
  {
    if (stack_name[i] == TIMER_STACK_SEPARATOR)
    {
      level++;
    }
  }
  return level;
}

std::string
get_leaf_name(const std::string &stack_name)
{
  int pos = stack_name.find_last_of(TIMER_STACK_SEPARATOR);
  if (pos == std::string::npos)
  {
    return stack_name;
  }

  return stack_name.substr(pos+1, stack_name.length()-pos);
}

void TimerManagerClass::collate_stack_profile(Communicate *comm, StackProfileData &p)
{
#ifdef USE_STACK_TIMERS
  // Put stacks from all timers into one data structure
  // By naming the timer stacks as 'timer1/timer2', etc, the ordering done by the
  // map's keys will also place the stacks in depth-first order.
  // The order in which sibling timers are encountered in the code is not
  // preserved. They will be ordered alphabetically instead.
  std::map<std::string, ProfileData> all_stacks;
  for(int i=0; i<TimerList.size(); ++i)
  {
    NewTimer &timer = *TimerList[i];
    std::map<std::string,double>::iterator stack_name_it = timer.get_per_stack_total_time().begin();
    for (; stack_name_it != timer.get_per_stack_total_time().end(); stack_name_it++)
    {
        ProfileData pd;
        const std::string &stack_name = stack_name_it->first;
        pd.time = timer.get_total(stack_name);
        pd.calls = timer.get_num_calls(stack_name);

        all_stacks[stack_name] += pd;
    }
  }

  // Fill in the output data structure (but don't compute exclusive time yet)
  std::map<std::string, ProfileData>::iterator si = all_stacks.begin();
  int idx = 0;
  for(;si != all_stacks.end(); ++si)
  {
    std::string stack_name = si->first;
    p.nameList[stack_name] = idx;
    p.names.push_back(stack_name);
    p.timeList.push_back(si->second.time);
    p.timeExclList.push_back(si->second.time);
    p.callList.push_back(si->second.calls);
    idx++;
  }

  // Subtract times of immediate children to get exclusive time
  for (idx = 0; idx < p.timeList.size(); idx++)
  {
    int start_level = get_level(p.names[idx]);
    for (int i = idx+1; i < p.timeList.size(); i++)
    {
      int level = get_level(p.names[i]);
      if (level == start_level+1)
      {
        p.timeExclList[idx] -= p.timeExclList[i];
      }
      if (level == start_level)
      {
        break;
      }
    }
  }
#endif
}

void
TimerManagerClass::print(Communicate* comm)
{
#if ENABLE_TIMERS
#ifdef USE_STACK_TIMERS
  if(comm == NULL || comm->rank() == 0)
  {
    printf("Stack timer profile\n");
  }
  print_stack(comm);
#else
  if(comm == NULL || comm->rank() == 0)
  {
    printf("\nFlat profile\n");
  }
  print_flat(comm);
#endif
#endif
}

void
TimerManagerClass::print_flat(Communicate* comm)
{
#if ENABLE_TIMERS
  FlatProfileData p;

  collate_flat_profile(comm, p);

  if(comm == NULL || comm->rank() == 0)
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
#if ENABLE_TIMERS
  StackProfileData p;

  collate_stack_profile(comm, p);

  if(comm == NULL || comm->rank() == 0)
  {
    printf("%-38s  %-9s  %-9s  %-10s  %-13s\n","Timer","Inclusive_time","Exclusive_time","Calls","Time_per_call");
    for (int i = 0; i < p.names.size(); i++)
    {
      std::string stack_name = p.names[i];
      int level = get_level(stack_name);
      std::string name = get_leaf_name(stack_name);
      std::string indent_str(2*level, ' ');
      std::string indented_str = indent_str + name;
      printf ("%-40s  %9.4f  %9.4f  %13ld  %16.9f\n"
              , indented_str.c_str()
              , p.timeList[i]
              , p.timeExclList[i]
              , p.callList[i]
              , p.timeList[i]/(static_cast<double>(p.callList[i])+std::numeric_limits<double>::epsilon()));
    }
  }
#endif
}

void
TimerManagerClass::output_timing(Communicate *comm, Libxml2Document &doc, xmlNodePtr root)
{
#if ENABLE_TIMERS
#ifdef USE_STACK_TIMERS
  StackProfileData p;

  collate_stack_profile(comm, p);

  if(comm == NULL || comm->rank() == 0)
  {
    xmlNodePtr timing_root = doc.addChild(root, "timing");
    std::vector<xmlNodePtr> node_stack;
    node_stack.push_back(timing_root);
    xmlNodePtr current_root = timing_root;

    for (int i = 0; i < p.names.size(); i++)
    {
      std::string stack_name = p.names[i];
      int level = get_level(stack_name);
      std::string name = get_leaf_name(stack_name);

      std::string indent_str(2*level, ' ');

      xmlNodePtr timer = doc.addChild(current_root, "timer");
      doc.addChild(timer, "name", name);
      doc.addChild(timer, "time_incl", p.timeList[i]);
      doc.addChild(timer, "time_excl", p.timeExclList[i]);
      doc.addChild(timer, "calls", p.callList[i]);

      int next_level = level;
      if (i+1 < p.names.size())
      {
        next_level = get_level(p.names[i+1]);
      }

      if (next_level > level)
      {
        xmlNodePtr next_node = doc.addChild(timer, "includes");
        node_stack.push_back(next_node);
        current_root = next_node;
      }
      if (next_level < level)
      {
        for (int j = 0; j < level-next_level; j++)
        {
          node_stack.pop_back();
          current_root = node_stack.back();
        }
      }
    }
  }

#endif
#endif
}

void
NewTimer::set_active_by_timer_threshold(const timer_levels threshold)
{
  if (timer_level <= threshold)
  {
    active = true;
  }
  else
  {
    active = false;
  }

}

}
