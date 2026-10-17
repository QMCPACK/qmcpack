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
/** @brief Coordinates application-wide execution timing and early termination signals.
 *
 * RunTimeManager tracks the total wallclock time since construction (typically at static
 * application startup) and maintains a shared stop signal to halt downstream calculations
 * cleanly when runtime limits or stop files are encountered.
 *
 * @tparam CLOCK Clock type providing now() and time_point (defaults to ChronoClock; FakeChronoClock for testing).
 */
template<class CLOCK = ChronoClock>
class RunTimeManager
{
public:
  /** @brief Get the elapsed wallclock time in seconds since this manager was initialized.
   *  @return Elapsed time in seconds.
   */
  inline double elapsed() const
  {
    std::chrono::duration<double> elapsed = CLOCK::now() - start_time;
    return elapsed.count();
  }
  /// Default constructor, records start timestamp.
  RunTimeManager() : start_time(CLOCK::now()), need_to_stop_(false) {}

  /** @brief Check whether an early simulation stop has been flagged.
   *  @return true if an early stop was requested by any driver or controller.
   */
  bool isStopNeeded() const { return need_to_stop_; }

  /// Flag that the simulation must stop early, preventing subsequent sections from running.
  void markStop() { need_to_stop_ = true; }

private:
  /// The reference timestamp taken at instantiation.
  const typename CLOCK::time_point start_time;
  /// Global early-termination flag.
  bool need_to_stop_{false};
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

  /** @brief Check whether there is sufficient remaining walltime for another block or iteration.
   *
   * Checks if total elapsed time has exceeded MaxCPUSecs, or if projected time for the next iteration
   * (average iteration time scaled by loop margin + runtime safety padding) exceeds remaining time.
   *
   * @return MAX_SECONDS_PASSED if hard time limit is reached,
   *         NOT_ENOUGH_TIME if remaining time cannot safely accommodate the next iteration,
   *         or CONTINUE if sufficient time remains.
   */
  StopStatus enough_time_for_next_iteration() const;

  /** @brief Check if an explicit stop file (.STOP) exists on disk.
   *  @return STOP_FILE if the stop file exists, or CONTINUE otherwise.
   */
  StopStatus stop_file_requested() const;

public:
  /** constructor
   * @param rm the RunTimeManager attached to
   * @param maxCPUSecs maxmimal allowed seconds from the beginning of the simulation
   * @param stop_file_prefix the prefix of the stop file
   * @param cleanup if true, clean up stop files left from previous simulations. Rank 0 handle this.
   */
  RunTimeControl(RunTimeManager<CLOCK>& rm, int maxCPUSecs, const std::string& stop_file_prefix, bool cleanup = true);

  /** check if the run needs to stop because of walltime or stop control file.
   * Call before each driver block or before additional loop iterations.
   */
  bool checkStop(LoopTimer<CLOCK>& loop_timer);

  /// generate terse progress messages
  std::string generateProgressMessage(const LoopTimer<CLOCK>& loop_timer,
                                      const std::string& driverName,
                                      int block,
                                      int num_blocks) const;
  /// Generate a stop message with the completed block or iteration count.
  std::string generateStopMessage(const std::string& driverName,
                                  int completed_iterations,
                                  const std::string& iteration_name = "block") const;

  // for testing
  void runtime_padding(int runtime_padding) { m_runtime_safety_padding = runtime_padding; }
  void loop_margin(int loopMargin) { m_loop_margin = loopMargin; }
};

extern template class RunTimeControl<ChronoClock>;
extern template class RunTimeControl<FakeChronoClock>;

} // namespace qmcplusplus
#endif
