//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_GLOBAL_OBJECTS_H
#define QMCPLUSPLUS_GLOBAL_OBJECTS_H

#include <ostream>
#include <string>

namespace qmcplusplus
{
/** class to definte global variables to keep track a run
 */
struct QMCState
{
  ///size of memory allocated in byte per MPI
  size_t memory_allocated{0};
  ///init for <qmc/> section
  int qmc_counter{0};
  ///true, if density is used, e.g. MPC
  bool use_density{false};
  ///true, if it is a dryrun
  bool dryrun{false};
  ///true, print out file
  bool io_node{true};

  ///constructor
  QMCState();
  ///initialize options from the command-line
  void initialize(int argc, char** argv);
  ///print command-line options
  void print_options(std::ostream& os);
  /** print memory increase
   * @param who the name of the class/function calling this
   * @param before memory_allocated before calling print
   */
  void print_memory_change(const std::string& who, size_t before);

  /// Print git info (commit hash, etc) if project was build from git repository
  void print_git_info_if_present(std::ostream& os);

  /** Setup the filename infix used for ensemble jobs.
   *  When num_jobs > 1, sets the infix to ".gXXX" based on my_job_id (e.g., ".g001" for my_job_id = 1).
   *  Otherwise, clears the infix to an empty string.
   *  @param my_job_id ID of the current job in the ensemble
   *  @param num_jobs total number of ensemble jobs
   */
  void setupFilenameInfixForEnsembleJobs(int my_job_id, int num_jobs);

  /// Return the filename infix for ensemble jobs (e.g., ".g001" if num_jobs > 1, or empty string otherwise).
  const std::string& get_filename_infix_for_ensemble_jobs() const noexcept { return filename_infix_for_ensemble_jobs_; }

private:
  /// Filename infix for ensemble jobs (e.g., ".g001" for job 1, or empty string if not an ensemble run)
  std::string filename_infix_for_ensemble_jobs_;
};

///a unique QMCState during a run
extern QMCState qmc_common;

} // namespace qmcplusplus

#endif
