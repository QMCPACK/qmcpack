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
template<class CLOCK = ChronoClock>
class RunTimeManager
{
public:
  inline double elapsed()
  {
    std::chrono::duration<double> elapsed = CLOCK::now() - start_time;
    return elapsed.count();
  }
  // Initialize the start time at static class initialization time
  RunTimeManager() : start_time(CLOCK::now()) {}

  bool isStopNeeded() const { return need_to_stop_; }
  void markStop() { need_to_stop_ = true; }

private:
  const typename CLOCK::time_point start_time;
  bool need_to_stop_;
};

extern RunTimeManager<ChronoClock> run_time_manager;

extern template class RunTimeManager<ChronoClock>;
extern template class RunTimeManager<FakeChronoClock>;

template<class CLOCK = ChronoClock>
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
  typename CLOCK::time_point start_time;
  double total_time;
};

extern template class LoopTimer<ChronoClock>;
extern template class LoopTimer<FakeChronoClock>;

template<class CLOCK = ChronoClock>
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
  bool stop_file_requested();

public:
  /** constructor
   * @param rm the RunTimeManager attached to
   * @param maxCPUSecs maxmimal allowed seconds from the beginning of the simulation
   * @param stop_file_prefix the prefix of the stop file
   * @param cleanup if true, clean up stop files left from previous simulations. Rank 0 handle this.
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

extern template class RunTimeControl<ChronoClock>;
extern template class RunTimeControl<FakeChronoClock>;

} // namespace qmcplusplus
#endif
