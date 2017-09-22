//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file RunTimeManager.h
 * @brief Class for determining elapsed run time enabling simulations to adjust to time limits.
 */
#ifndef QMCPLUSPLUS_RUNTIME_MANAGER_H
#define QMCPLUSPLUS_RUNTIME_MANAGER_H

#include <Utilities/Clock.h>


namespace qmcplusplus
{

class RunTimeManagerClass
{
public:
  void start()  { start_time = cpu_clock(); }
  double elapsed() { return cpu_clock() - start_time; }
  // Initialize the start time at static class initialization time
  RunTimeManagerClass() { start(); }
private:
  double start_time;
};

extern RunTimeManagerClass RunTimeManager;

class LoopTimer
{
public:
  void start() { start_time = cpu_clock(); }
  void stop() { nloop++; total_time += cpu_clock() - start_time; }
  double get_time_per_iteration();

  LoopTimer() : nloop(0), start_time(0.0), total_time(0.0) {}
private:
  int nloop;
  double start_time;
  double total_time;
};

}
#endif
