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
#include <string_view>

class Communicate;

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

enum class StopStatus
{
  CONTINUE,           // keep running
  MAX_SECONDS_PASSED, // already passed max_seconds
  NOT_ENOUGH_TIME,    // not enough time for next iteration
  STOP_FILE,          // requested stop from a file
};

class ProjectData;

/** @brief Scope-level runtime and stop-condition controller.
 *
 * RunTimeControl handles early simulation termination at discrete checkpoints or scope boundaries
 * (e.g. before launching a QMC section in QMCMain::execute() or between distinct QMC subsections
 * inside a loop). It inspects whether the hard wallclock limit (from ProjectData::getMaxCPUSeconds())
 * has been exceeded or whether a user-created stop file (`<title>.STOP`) is present on disk.
 *
 * Rank 0 conducts the walltime and filesystem checks, broadcasts the decision to all ranks in the
 * communicator, prints an informative message to app_log(), and flags the attached RunTimeManager
 * via markStop(). If RunTimeManager already indicates that a stop is needed, checkStop() acts as a
 * lightweight no-op returning true immediately.
 *
 * Intended use:
 * Use RunTimeControl at coarse boundaries outside repetitive loops where iteration timing / projection
 * is not applicable, for example:
 * @code
 * RunTimeControl<> runtime_control(run_time_manager, project_data, "QMCMain", comm);
 * for (auto& action : qmc_actions)
 * {
 *   if (runtime_control.checkStop("executeAction"))
 *     break;
 *   executeAction(action);
 * }
 * @endcode
 *
 * @tparam CLOCK Clock type (defaults to ChronoClock; FakeChronoClock used in testing).
 */
template<class CLOCK = ChronoClock>
class RunTimeControl
{
  const ProjectData& project_data_;
  const std::string scope_name_;
  RunTimeManager<CLOCK>& runtimeManager;
  Communicate& comm_;

  double m_elapsed{0.0};
  std::string stop_message_;

  /// Return the stop control filename for the current project.
  std::string getStopFilename() const;

  /// Generate a stop message.
  std::string generateStopMessage(StopStatus stop_status, std::string_view phase) const;

public:
  /** @brief Constructor.
   *  @param rm the RunTimeManager tracking application-wide elapsed time and stop state
   *  @param project_data ProjectData containing project title (for the `.STOP` filename) and max CPU seconds
   *  @param scope_name descriptive name of the scope for diagnostic logging (e.g. "QMCMain")
   *  @param comm communicator used for rank checks and broadcasting stop signals across processes
   */
  RunTimeControl(RunTimeManager<CLOCK>& rm,
                 const ProjectData& project_data,
                 const std::string& scope_name,
                 Communicate& comm);

  /** @brief Check whether the run needs to stop because of walltime or stop control file.
   *  Rank 0 performs the time/file checks, broadcasts to comm, logs the message,
   *  and marks stop on runtimeManager. If runtimeManager already needs stop, returns true as a no-op.
   *  @param phase execution phase name (e.g. "executeQMC") included in the stop message (must not be empty)
   *  @return true if stop is requested, or false if continuing.
   */
  bool checkStop(std::string_view phase);

  /// Retrieve the generated stop message (primarily for unit testing; empty on non-root ranks).
  const std::string& get_stop_message() const { return stop_message_; }
};

extern template class RunTimeControl<ChronoClock>;
extern template class RunTimeControl<FakeChronoClock>;

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

/** @brief Iteration- and block-level runtime and stop-condition controller.
 *
 * LoopRunTimeControl manages early termination within iterative computational loops (e.g. Monte Carlo
 * blocks in QMC drivers, wave function optimization iterations, or multi-iteration QMC loops). In addition
 * to checking hard wallclock limits and the `<title>.STOP` control file, LoopRunTimeControl uses a LoopTimer
 * to estimate whether sufficient remaining walltime exists to safely execute another iteration or block:
 *
 *   (loop_margin * average_iteration_time + runtime_safety_padding) <= remaining_walltime
 *
 * If insufficient time remains, or if a hard limit / stop file is encountered, rank 0 generates a detailed
 * termination message, broadcasts the stop signal to all ranks, logs the message, removes the consumed
 * `<title>.STOP` file, and flags the attached RunTimeManager via markStop().
 *
 * If RunTimeManager is already flagged for stopping, checkStop() returns true immediately as a silent no-op.
 *
 * Intended use:
 * Instantiate inside driver run() methods or iterative loops, calling checkStop() at the start of each block:
 * @code
 * LoopTimer<> loop_timer;
 * LoopRunTimeControl<> runtime_control(run_time_manager, project_data, "VMC driver", "block", comm);
 * for (int block = 0; block < num_blocks; ++block)
 * {
 *   if (runtime_control.checkStop(loop_timer, block))
 *     break;
 *
 *   {
 *     ScopeGuard<LoopTimer<>> block_timer(loop_timer);
 *     runBlock(block);
 *   }
 *   app_log() << runtime_control.generateProgressMessage(loop_timer, block, num_blocks);
 * }
 * @endcode
 *
 * @tparam CLOCK Clock type (defaults to ChronoClock; FakeChronoClock used in testing).
 */
template<class CLOCK = ChronoClock>
class LoopRunTimeControl
{
  const ProjectData& project_data_;
  const std::string scope_name_;
  const std::string iteration_name_;
  double m_runtime_safety_padding{30.0};
  double m_loop_margin{1.1};
  double m_loop_time{0.0};
  double m_elapsed{0.0};
  double m_remaining{0.0};
  RunTimeManager<CLOCK>& runtimeManager;
  Communicate& comm_;

  std::string stop_message_;

  /// Return the stop control filename for the current project.
  std::string getStopFilename() const;

  /** @brief Check whether there is sufficient remaining walltime for another block or iteration.
   *
   * Checks if total elapsed time has exceeded maximum CPU seconds, or if projected time for the next iteration
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

  /// Generate a stop message with the completed block or iteration count.
  std::string generateStopMessage(StopStatus stop_status, int completed_iterations) const;

public:
  /** @brief Constructor.
   *  @param rm the RunTimeManager tracking application-wide elapsed time and stop state
   *  @param project_data ProjectData containing project title (for the `.STOP` filename) and max CPU seconds
   *  @param scope_name descriptive name of the driver or loop scope for logging (e.g. "DMC driver", "executeLoop")
   *  @param iteration_name name of the loop unit (e.g. "block", "loop iteration", "optimization iteration")
   *  @param comm communicator used for rank checks and broadcasting stop signals across processes
   */
  LoopRunTimeControl(RunTimeManager<CLOCK>& rm,
                     const ProjectData& project_data,
                     const std::string& scope_name,
                     const std::string& iteration_name,
                     Communicate& comm);

  /** @brief Check if the run needs to stop because of walltime or stop control file.
   *  Rank 0 performs the time/file checks, broadcasts to comm, logs the message,
   *  and marks stop on runtimeManager. If runtimeManager already needs stop, returns true as a no-op.
   *  Call before each driver block or before additional loop iterations.
   *  @param loop_timer LoopTimer tracking the duration of completed iterations
   *  @param completed_iterations count of already completed iterations or blocks
   *  @return true if stop is requested, or false if continuing.
   */
  bool checkStop(LoopTimer<CLOCK>& loop_timer, int completed_iterations);

  /** @brief Generate terse progress messages at milestone intervals (e.g., 0%, 25%, 50%, 75%, 100%).
   *  @param loop_timer LoopTimer tracking iteration durations
   *  @param block 0-based index of current block
   *  @param num_blocks total number of blocks
   *  @return formatted progress message string, or empty string if not at a milestone block.
   */
  std::string generateProgressMessage(const LoopTimer<CLOCK>& loop_timer, int block, int num_blocks) const;

  // for testing
  const std::string& get_stop_message() const { return stop_message_; }
  void set_runtime_padding(double runtime_padding) { m_runtime_safety_padding = runtime_padding; }
  void set_loop_margin(double loop_margin) { m_loop_margin = loop_margin; }
};

extern template class LoopRunTimeControl<ChronoClock>;
extern template class LoopRunTimeControl<FakeChronoClock>;

} // namespace qmcplusplus
#endif
