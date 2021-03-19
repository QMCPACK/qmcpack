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

#include "Utilities/Clock.h"
#include <string>


namespace qmcplusplus
{
template<class CLOCK = CPUClock>
class RunTimeManager
{
public:
  inline double elapsed() { return CLOCK()() - start_time; }
  // Initialize the start time at static class initialization time
  RunTimeManager() : start_time(CLOCK()()) {}

  bool isStopNeeded() const { return need_to_stop_; }
  void markStop() { need_to_stop_ = true; }

private:
  const double start_time;
  bool need_to_stop_;
};

extern RunTimeManager<CPUClock> run_time_manager;

extern template class RunTimeManager<CPUClock>;
extern template class RunTimeManager<FakeCPUClock>;

template<class CLOCK = CPUClock>
class LoopTimer
{
public:
  LoopTimer();
  void start();
  void stop();
  double get_time_per_iteration() const;

private:
  int nloop;
  bool ticking;
  double start_time;
  double total_time;
};

extern template class LoopTimer<CPUClock>;
extern template class LoopTimer<FakeCPUClock>;

template<class CLOCK = CPUClock>
class RunTimeControl
{
  const int MaxCPUSecs;
  double m_runtime_safety_padding;
  double m_loop_margin;
  double m_loop_time;
  double m_elapsed;
  double m_remaining;
  RunTimeManager<CLOCK>& runtimeManager;
  /// the prefix of the stop file (stop_file_prefix + ".STOP")
  const std::string stop_filename_;

  enum class StopStatus
  {
    CONTINUE,           // keep running
    MAX_SECONDS_PASSED, // all already passed max_seconds
    NOT_ENOUGH_TIME,    // not enough time for next iteration
    STOP_FILE,          // reqsuted stop from a file
  } stop_status_;

  bool enough_time_for_next_iteration(LoopTimer<CLOCK>& loop_timer);
  bool stop_file_reqeusted();

public:
  /** constructor
   * @param rm the RunTimeManager attached to
   * @param maxCPUSecs maxmimal allowed seconds from the beginning of the simulation
   * @param stop_file_prefix the prefix of the stop file
   * @param cleanup if true, clean up stop files left from previous simulations. Fine control needed for ensemble runs.
   */
  RunTimeControl(RunTimeManager<CLOCK>& rm, int maxCPUSecs, const std::string& stop_file_prefix, bool cleanup = true);

  /** check if the run needs to stop because of walltime or stop control file.
   * it should be called at the end of each block in a driver.
   */
  bool checkStop(LoopTimer<CLOCK>& loop_timer);

  /// generate stop message explaining why
  std::string generateStopMessage(const std::string& driverName, int block) const;

  // for testing
  void runtime_padding(int runtime_padding) { m_runtime_safety_padding = runtime_padding; }
  void loop_margin(int loopMargin) { m_loop_margin = loopMargin; }
};

extern template class RunTimeControl<CPUClock>;
extern template class RunTimeControl<FakeCPUClock>;

} // namespace qmcplusplus
#endif
