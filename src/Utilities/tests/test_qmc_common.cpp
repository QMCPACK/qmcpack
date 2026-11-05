//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include <catch2/catch_test_macros.hpp>
#include "qmc_common.h"

namespace qmcplusplus
{

TEST_CASE("QMCState::setupFilenameInfixForEnsembleJobs", "[utilities]")
{
  QMCState state;
  CHECK(state.get_filename_infix_for_ensemble_jobs().empty());

  // Single job case (num_jobs == 1) should produce an empty infix
  state.setupFilenameInfixForEnsembleJobs(0, 1);
  CHECK(state.get_filename_infix_for_ensemble_jobs().empty());

  // Ensemble jobs (num_jobs > 1)
  state.setupFilenameInfixForEnsembleJobs(0, 2);
  CHECK(state.get_filename_infix_for_ensemble_jobs() == ".g000");

  state.setupFilenameInfixForEnsembleJobs(1, 2);
  CHECK(state.get_filename_infix_for_ensemble_jobs() == ".g001");

  state.setupFilenameInfixForEnsembleJobs(12, 16);
  CHECK(state.get_filename_infix_for_ensemble_jobs() == ".g012");

  state.setupFilenameInfixForEnsembleJobs(123, 150);
  CHECK(state.get_filename_infix_for_ensemble_jobs() == ".g123");

  // Setting back to a single job should clear the infix
  state.setupFilenameInfixForEnsembleJobs(0, 1);
  CHECK(state.get_filename_infix_for_ensemble_jobs().empty());
}

} // namespace qmcplusplus
