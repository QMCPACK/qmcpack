//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_VALIDQMCINPUTSECTIONS_H
#define QMCPLUSPLUS_VALIDQMCINPUTSECTIONS_H

/** @file
 *  @brief avoids many repeated xml heredoc sections
 */

#include <array>

namespace qmcplusplus
{
namespace testing
{
// clang-format: off
constexpr std::array<const char*, 3> valid_vmc_input_sections{
    R"(
  <qmc method="vmc" move="pbyp">
    <estimator name="LocalEnergy" hdf5="no" />
    <parameter name="walkers">                1 </parameter>
    <parameter name="stepsbetweensamples">    1 </parameter>
    <parameter name="warmupSteps">            5 </parameter>
    <parameter name="substeps">               5 </parameter>
    <parameter name="steps">                  1 </parameter>
    <parameter name="blocks">                 2 </parameter>
    <parameter name="timestep">             1.0 </parameter>
    <parameter name="usedrift">              no </parameter>
  </qmc>
)",
     R"(
  <qmc method="vmc_batch" move="pbyp">
    <parameter name="crowds">                 8 </parameter>
    <estimator name="LocalEnergy" hdf5="no" />
    <parameter name="walkers">                32 </parameter>
    <parameter name="stepsbetweensamples">    1 </parameter>
    <parameter name="warmupSteps">            5 </parameter>
    <parameter name="substeps">               5 </parameter>
    <parameter name="steps">                  1 </parameter>
    <parameter name="blocks">                 2 </parameter>
    <parameter name="timestep">             1.0 </parameter>
    <parameter name="usedrift">              no </parameter>
  </qmc>
)",
     R"(
  <qmc method="vmc_batch" move="pbyp">
    <parameter name="crowds">                 1 </parameter>
    <estimator name="LocalEnergy" hdf5="no" />
    <parameter name="walkers">                1 </parameter>
    <parameter name="stepsbetweensamples">    1 </parameter>
    <parameter name="warmupSteps">            5 </parameter>
    <parameter name="substeps">               5 </parameter>
    <parameter name="steps">                  1 </parameter>
    <parameter name="blocks">                 2 </parameter>
    <parameter name="timestep">             1.0 </parameter>
    <parameter name="usedrift">              no </parameter>
  </qmc>
)"};


// to avoid creating a situation where section test xml is in two places
constexpr int valid_vmc_input_vmc_index       = 0;
constexpr int valid_vmc_input_vmc_batch_index = 1;
constexpr int valid_vmc_input_vmc_tiny_index = 2;

constexpr std::array<const char*, 2> valid_dmc_input_sections{
    R"(
  <qmc method="dmc" move="pbyp" gpu="yes">
    <estimator name="LocalEnergy" hdf5="no" />
    <parameter name="targetwalkers">        256 </parameter>
    <parameter name="warmupSteps">          100 </parameter>
    <parameter name="steps">                 10 </parameter>
    <parameter name="blocks">               100 </parameter>
    <parameter name="timestep">            0.01 </parameter>
    <parameter name="reconfiguration">       no </parameter>
    <parameter name="nonlocalmoves">         no </parameter>
  </qmc>
)",
     R"(
  <qmc method="dmc_batch" move="pbyp">
    <parameter name="crowds">                 8 </parameter>
    <estimator name="LocalEnergy" hdf5="no" />
    <parameter name="walkers">                32 </parameter>
    <parameter name="warmupSteps">            5 </parameter>
    <parameter name="substeps">               5 </parameter>
    <parameter name="steps">                  1 </parameter>
    <parameter name="blocks">                 2 </parameter>
    <parameter name="timestep">             1.0 </parameter>
    <parameter name="usedrift">              no </parameter>
  </qmc>
)"};

// to avoid creating a situation where section test xml is in two places
constexpr int valid_dmc_input_dmc_index       = 0;
constexpr int valid_dmc_input_dmc_batch_index = 1;

// clang-format: on
} // namespace testing
} // namespace qmcplusplus

#endif
