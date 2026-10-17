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
#include <fstream>
#include <iomanip>
#include <cstdio>

namespace qmcplusplus
{
RunTimeManager<ChronoClock> run_time_manager;

template class RunTimeManager<ChronoClock>;
template class RunTimeManager<FakeChronoClock>;

template<class CLOCK>
LoopTimer<CLOCK>::LoopTimer() : nloop(0), ticking(false), total_time(0.0)
{}

template<class CLOCK>
void LoopTimer<CLOCK>::start()
{
  if (ticking)
    throw std::runtime_error("LoopTimer started already!");
  start_time = CLOCK::now();
  ticking    = true;
}

template<class CLOCK>
void LoopTimer<CLOCK>::stop()
{
  if (!ticking)
    throw std::runtime_error("LoopTimer didn't start but called stop!");
  nloop++;
  std::chrono::duration<double> elapsed = CLOCK::now() - start_time;
  total_time += elapsed.count();

  ticking = false;
}

template<class CLOCK>
double LoopTimer<CLOCK>::get_time_per_iteration() const
{
  if (nloop > 0)
    return total_time / nloop;
  return 0.0;
}

template class LoopTimer<ChronoClock>;
template class LoopTimer<FakeChronoClock>;

template<class CLOCK>
RunTimeControl<CLOCK>::RunTimeControl(RunTimeManager<CLOCK>& rm,
                                      int maxCPUSecs,
                                      const std::string& stop_file_prefix,
                                      bool cleanup)
    : MaxCPUSecs(maxCPUSecs),
      runtimeManager(rm),
      stop_filename_(stop_file_prefix + ".STOP"),
      stop_status_(StopStatus::CONTINUE)
{
  if (stop_file_prefix.empty())
    throw std::runtime_error("Stop control filename prefix must not be empty!");

  if (cleanup)
  {
    std::remove(stop_filename_.c_str());
    if (std::ifstream(stop_filename_.c_str()))
      throw std::runtime_error("Failed to delete the existing stop control file \"" + stop_filename_ +
                               "\", cannot continue!");
  }

  m_runtime_safety_padding = 30.0; // generous 30 seconds to allow for shut down?
  m_loop_margin            = 1.1;  // 10% margin on average loop time?
}

template<class CLOCK>
typename RunTimeControl<CLOCK>::StopStatus RunTimeControl<CLOCK>::enough_time_for_next_iteration() const
{
  if (m_elapsed >= MaxCPUSecs)
    return StopStatus::MAX_SECONDS_PASSED;

  if ((m_loop_margin * m_loop_time + m_runtime_safety_padding) > m_remaining)
    return StopStatus::NOT_ENOUGH_TIME;

  return StopStatus::CONTINUE;
}

template<class CLOCK>
typename RunTimeControl<CLOCK>::StopStatus RunTimeControl<CLOCK>::stop_file_requested() const
{
  if (std::ifstream(stop_filename_.c_str()))
    return StopStatus::STOP_FILE;
  else
    return StopStatus::CONTINUE;
}

template<class CLOCK>
bool RunTimeControl<CLOCK>::checkStop(LoopTimer<CLOCK>& loop_timer)
{
  m_loop_time = loop_timer.get_time_per_iteration();
  m_elapsed   = runtimeManager.elapsed();
  m_remaining = MaxCPUSecs - m_elapsed;

  stop_status_ = enough_time_for_next_iteration();
  bool need_to_stop = stop_status_ != StopStatus::CONTINUE;

  if (!need_to_stop)
  {
    stop_status_ = stop_file_requested();
    need_to_stop = stop_status_ != StopStatus::CONTINUE;
  }

  return need_to_stop;
}

template<class CLOCK>
std::string RunTimeControl<CLOCK>::generateProgressMessage(const LoopTimer<CLOCK>& loop_timer,
                                                           const std::string& driverName,
                                                           int block,
                                                           int num_blocks) const
{
  std::stringstream log;
  if (block == 0 || block + 1 == num_blocks / 4 || block + 1 == num_blocks / 2 || block + 1 == (num_blocks * 3) / 4 ||
      block + 1 == num_blocks)
  {
    log << "Completed block " << std::setw(4) << block + 1 << " of " << num_blocks << " average "
        << std::setprecision(4) << loop_timer.get_time_per_iteration() << " secs/block after " << std::setprecision(4)
        << runtimeManager.elapsed() << " secs" << std::endl;
  }
  return log.str();
}

template<class CLOCK>
std::string RunTimeControl<CLOCK>::generateStopMessage(const std::string& driverName,
                                                       int completed_iterations,
                                                       const std::string& iteration_name) const
{
  std::stringstream log;
  log << "RunTimeControl takes action in " << driverName;
  if (iteration_name == "block")
    log << " driver";
  log << "." << std::endl;
  if (stop_status_ == StopStatus::MAX_SECONDS_PASSED)
    log << "Time limit reached. Stopping after " << iteration_name << " " << completed_iterations << std::endl
        << "Hard limit (seconds) " << MaxCPUSecs << ", elapsed (seconds) " << m_elapsed << std::endl;
  else if (stop_status_ == StopStatus::NOT_ENOUGH_TIME)
  {
    log << "Insufficient time for next " << iteration_name << ". Stopping after " << iteration_name << " "
        << completed_iterations << std::endl;
    log << "  Iteration time per " << driverName << " " << iteration_name << " (seconds) = " << m_loop_time
        << std::endl;
    log << "  Elapsed   time (seconds) = " << m_elapsed << std::endl;
    log << "  Remaining time (seconds) = " << m_remaining << std::endl;
  }
  else if (stop_status_ == StopStatus::STOP_FILE)
    log << "Stop requested from the control file \"" + stop_filename_ + "\", stopping after " << iteration_name << " "
        << completed_iterations << std::endl;
  else
    throw std::runtime_error("Unidentified stop status!");

  return log.str();
}

template class RunTimeControl<ChronoClock>;
template class RunTimeControl<FakeChronoClock>;

} // namespace qmcplusplus
