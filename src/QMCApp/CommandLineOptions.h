//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_COMMAND_LINE_OPTIONS_H
#define QMCPLUSPLUS_COMMAND_LINE_OPTIONS_H

#include <string>
#include <vector>
#include <optional>
#include <ostream>
#include "Platforms/Host/OutputManager.h"

namespace qmcplusplus
{

class CommandLineOptions
{
public:
  CommandLineOptions() = default;

  /** parse command line arguments
   * @param argc argument count
   * @param argv argument vector
   */
  void parse(int argc, char** argv);

  /// Return true if debug output is requested
  bool get_debug() const noexcept { return debug_; }

  /// Return true if dryrun is requested
  bool get_dryrun() const noexcept { return dryrun_; }

  /// Return true if noprint is requested
  bool get_noprint() const noexcept { return noprint_; }

  /// Return true if help is requested
  bool get_help() const noexcept { return help_; }

  /// Return true if version is requested
  bool get_version() const noexcept { return version_; }

  /// Return true if execution should normally stop (e.g. help, version)
  bool should_normal_stop() const noexcept { return help_ || version_; }

  /// Return true if timers are enabled
  bool get_enable_timers() const noexcept { return enable_timers_; }

  /// Return timer threshold level
  const std::string& get_timer_level() const noexcept { return timer_level_; }

  /// Return verbosity level if specified
  const std::optional<Verbosity>& get_verbosity() const noexcept { return verbosity_; }

  /// Return virtual node ranks if specified
  const std::optional<int>& get_vnode_ranks() const noexcept { return vnode_ranks_; }

  /// Return parsed input XML filenames
  const std::vector<std::string>& get_inputs() const noexcept { return inputs_; }

  /// Print version and usage information to output stream
  void print_version_and_usage(std::ostream& os) const;

private:
  bool debug_{false};
  bool dryrun_{false};
  bool noprint_{false};
  bool help_{false};
  bool version_{false};
  bool enable_timers_{false};
  std::string timer_level_;
  std::optional<Verbosity> verbosity_;
  std::optional<int> vnode_ranks_;
  std::vector<std::string> inputs_;
};

} // namespace qmcplusplus

#endif
