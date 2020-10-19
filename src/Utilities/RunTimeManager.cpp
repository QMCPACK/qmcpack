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


/** @file RunTimeManager.cpp
 *  @brief Class for determining elapsed run time enabling simulations to adjust to time limits.

 */
#include "RunTimeManager.h"
#include <sstream>


namespace qmcplusplus
{
RunTimeManager<CPUClock> run_time_manager;

template class RunTimeManager<CPUClock>;
template class RunTimeManager<FakeCPUClock>;

template<class CLOCK>
double LoopTimer<CLOCK>::get_time_per_iteration()
{
  if (nloop > 0)
  {
    return total_time / nloop;
  }
  return 0.0;
}

template class LoopTimer<CPUClock>;
template class LoopTimer<FakeCPUClock>;

template<class CLOCK>
bool RunTimeControl<CLOCK>::enough_time_for_next_iteration(LoopTimer<CLOCK>& loop_timer)
{
  m_loop_time      = loop_timer.get_time_per_iteration();
  m_elapsed        = runtimeManager.elapsed();
  m_remaining      = MaxCPUSecs - m_elapsed;
  bool enough_time = true;
  if ((m_loop_margin * m_loop_time + m_runtime_safety_padding) > m_remaining)
    enough_time = false;

  return enough_time;
}

template<class CLOCK>
std::string RunTimeControl<CLOCK>::time_limit_message(const std::string& driverName, int block)
{
  std::stringstream log;
  log << "Time limit reached for " << driverName << ", stopping after block " << block - 1 << std::endl;
  log << "  Iteration time per " << driverName << " block (seconds) = " << m_loop_time << std::endl;
  log << "  Elapsed time (seconds)       = " << m_elapsed << std::endl;
  log << "  Remaining time (seconds)      = " << m_remaining << std::endl;
  return log.str();
}

template class RunTimeControl<CPUClock>;
template class RunTimeControl<FakeCPUClock>;

} // namespace qmcplusplus
