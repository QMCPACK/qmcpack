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
#include "ProjectData.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "OutputManager.h"
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
RunTimeControl<CLOCK>::RunTimeControl(RunTimeManager<CLOCK>& rm,
                                      const ProjectData& project_data,
                                      const std::string& scope_name,
                                      Communicate& comm)
    : project_data_(project_data), scope_name_(scope_name), runtimeManager(rm), comm_(comm), m_elapsed(0.0)
{
  if (project_data_.getTitle().empty())
    throw std::runtime_error("Stop control filename prefix must not be empty!");
}

template<class CLOCK>
std::string RunTimeControl<CLOCK>::generateStopMessage(StopStatus stop_status) const
{
  const int max_cpu_secs          = project_data_.getMaxCPUSeconds();
  const std::string stop_filename = project_data_.getTitle() + ".STOP";

  std::stringstream log;
  log << "RunTimeControl takes action in " << scope_name_ << "." << std::endl;
  if (stop_status == StopStatus::MAX_SECONDS_PASSED)
    log << "Time limit reached." << std::endl
        << "Hard limit (seconds) " << max_cpu_secs << ", elapsed (seconds) " << m_elapsed << std::endl;
  else if (stop_status == StopStatus::STOP_FILE)
    log << "Stop requested from the control file \"" << stop_filename << "\"" << std::endl;

  return log.str();
}

template<class CLOCK>
bool RunTimeControl<CLOCK>::checkStop()
{
  if (runtimeManager.isStopNeeded())
  {
    stop_message_.clear();
    return true;
  }

  bool stop_requested = false;
  if (comm_.rank() == 0)
  {
    m_elapsed                       = runtimeManager.elapsed();
    const int max_cpu_secs          = project_data_.getMaxCPUSeconds();
    const std::string stop_filename = project_data_.getTitle() + ".STOP";

    StopStatus stop_status = StopStatus::CONTINUE;
    if (m_elapsed >= max_cpu_secs)
      stop_status = StopStatus::MAX_SECONDS_PASSED;
    else if (std::ifstream(stop_filename.c_str()))
    {
      std::remove(stop_filename.c_str());
      stop_status = StopStatus::STOP_FILE;
    }

    if (stop_status != StopStatus::CONTINUE)
    {
      stop_message_  = generateStopMessage(stop_status);
      stop_requested = true;
    }
  }

  comm_.bcast(stop_requested);

  if (stop_requested)
  {
    app_log() << stop_message_;
    runtimeManager.markStop();
    return true;
  }
  else
  {
    stop_message_.clear();
    return false;
  }
}

template class RunTimeControl<ChronoClock>;
template class RunTimeControl<FakeChronoClock>;

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
LoopRunTimeControl<CLOCK>::LoopRunTimeControl(RunTimeManager<CLOCK>& rm,
                                              const ProjectData& project_data,
                                              const std::string& scope_name,
                                              const std::string& iteration_name,
                                              Communicate& comm)
    : project_data_(project_data),
      scope_name_(scope_name),
      iteration_name_(iteration_name),
      runtimeManager(rm),
      comm_(comm)
{
  if (project_data_.getTitle().empty())
    throw std::runtime_error("Stop control filename prefix must not be empty!");

  m_runtime_safety_padding = 30.0; // generous 30 seconds to allow for shut down?
  m_loop_margin            = 1.1;  // 10% margin on average loop time?
}

template<class CLOCK>
StopStatus LoopRunTimeControl<CLOCK>::enough_time_for_next_iteration() const
{
  if (m_elapsed >= project_data_.getMaxCPUSeconds())
    return StopStatus::MAX_SECONDS_PASSED;

  if ((m_loop_margin * m_loop_time + m_runtime_safety_padding) > m_remaining)
    return StopStatus::NOT_ENOUGH_TIME;

  return StopStatus::CONTINUE;
}

template<class CLOCK>
StopStatus LoopRunTimeControl<CLOCK>::stop_file_requested() const
{
  const std::string stop_filename = project_data_.getTitle() + ".STOP";
  if (std::ifstream(stop_filename.c_str()))
  {
    std::remove(stop_filename.c_str());
    return StopStatus::STOP_FILE;
  }
  else
    return StopStatus::CONTINUE;
}

template<class CLOCK>
bool LoopRunTimeControl<CLOCK>::checkStop(LoopTimer<CLOCK>& loop_timer, int completed_iterations)
{
  if (runtimeManager.isStopNeeded())
  {
    stop_message_.clear();
    return true;
  }

  bool stop_requested = false;
  if (comm_.rank() == 0)
  {
    m_loop_time = loop_timer.get_time_per_iteration();
    m_elapsed   = runtimeManager.elapsed();
    m_remaining = project_data_.getMaxCPUSeconds() - m_elapsed;

    StopStatus stop_status = enough_time_for_next_iteration();

    if (stop_status == StopStatus::CONTINUE)
      stop_status = stop_file_requested();

    if (stop_status != StopStatus::CONTINUE)
    {
      stop_message_  = generateStopMessage(stop_status, completed_iterations);
      stop_requested = true;
    }
  }

  comm_.bcast(stop_requested);

  if (stop_requested)
  {
    app_log() << stop_message_;
    runtimeManager.markStop();
    return true;
  }
  else
  {
    stop_message_.clear();
    return false;
  }
}

template<class CLOCK>
std::string LoopRunTimeControl<CLOCK>::generateProgressMessage(const LoopTimer<CLOCK>& loop_timer,
                                                               int block,
                                                               int num_blocks) const
{
  std::stringstream log;
  if (block == 0 || block + 1 == num_blocks / 4 || block + 1 == num_blocks / 2 || block + 1 == (num_blocks * 3) / 4 ||
      block + 1 == num_blocks)
  {
    log << "Completed " << iteration_name_ << " " << std::setw(4) << block + 1 << " of " << num_blocks << " average "
        << std::setprecision(4) << loop_timer.get_time_per_iteration() << " secs/" << iteration_name_ << " after "
        << std::setprecision(4) << runtimeManager.elapsed() << " secs" << std::endl;
  }
  return log.str();
}

template<class CLOCK>
std::string LoopRunTimeControl<CLOCK>::generateStopMessage(StopStatus stop_status, int completed_iterations) const
{
  std::stringstream log;
  log << "LoopRunTimeControl takes action in " << scope_name_ << "." << std::endl;
  if (stop_status == StopStatus::MAX_SECONDS_PASSED)
    log << "Time limit reached. Stopping before " << iteration_name_ << " " << completed_iterations << std::endl
        << "Hard limit (seconds) " << project_data_.getMaxCPUSeconds() << ", elapsed (seconds) " << m_elapsed << std::endl;
  else if (stop_status == StopStatus::NOT_ENOUGH_TIME)
  {
    log << "Insufficient time for the current " << iteration_name_ << ". Stopping before " << iteration_name_ << " "
        << completed_iterations << std::endl;
    log << "  Iteration time per " << scope_name_ << " " << iteration_name_ << " (seconds) = " << m_loop_time
        << std::endl;
    log << "  Elapsed   time (seconds) = " << m_elapsed << std::endl;
    log << "  Remaining time (seconds) = " << m_remaining << std::endl;
  }
  else if (stop_status == StopStatus::STOP_FILE)
    log << "Stop requested from the control file \"" + project_data_.getTitle() + ".STOP\", stopping before " << iteration_name_ << " "
        << completed_iterations << std::endl;
  else
    throw std::runtime_error("Unidentified stop status!");

  return log.str();
}

template class LoopRunTimeControl<ChronoClock>;
template class LoopRunTimeControl<FakeChronoClock>;

} // namespace qmcplusplus
