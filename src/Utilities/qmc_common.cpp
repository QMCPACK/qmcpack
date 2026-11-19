//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "qmc_common.h"
#include <array>
#include <cassert>
#include <cstdio>
#include <stdexcept>
#include "config.h"
#include "qmcpack_version.h"
#include "Message/Communicate.h"
#include "Platforms/Host/OutputManager.h"

namespace qmcplusplus
{
QMCState::QMCState() = default;

void QMCState::print_options(std::ostream& os)
{
  os << "  Global options " << std::endl;
  if (dryrun)
    os << "  dryrun : All QMC and loop sections will be skipped." << std::endl;
}

void QMCState::print_memory_change(const std::string& who, size_t before)
{
  before = memory_allocated - before;
  app_log() << "MEMORY increase " << (before >> 20) << " MB " << who << std::endl;
}

void QMCState::print_git_info_if_present(std::ostream& os)
{
#ifdef QMCPACK_GIT_BRANCH
  os << std::endl;
  os << "  Git branch: " << QMCPACK_GIT_BRANCH << std::endl;
  os << "  Last git commit: " << QMCPACK_GIT_HASH << std::endl;
  os << "  Last git commit date: " << QMCPACK_GIT_COMMIT_LAST_CHANGED << std::endl;
  os << "  Last git commit subject: " << QMCPACK_GIT_COMMIT_SUBJECT << std::endl;
#endif
}

void QMCState::setupFilenameInfixForEnsembleJobs(int my_job_id, int num_jobs)
{
  assert(num_jobs > 0);
  assert(my_job_id >= 0);

  if (num_jobs > 1)
  {
    std::array<char, 16> infix;
    int len = std::snprintf(infix.data(), infix.size(), ".g%03d", my_job_id);
    if (len < 0)
      throw std::runtime_error("Error generating filename_infix_for_ensemble_jobs");
    filename_infix_for_ensemble_jobs_ = std::string(infix.data(), len);
  }
  else
    filename_infix_for_ensemble_jobs_.clear();
}


QMCState qmc_common;

} // namespace qmcplusplus
