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
#include <string_view>

namespace qmcplusplus
{
namespace testing
{

class VmcLegacyInput
{
public:
  static std::string_view getXml() { return xml; }

private:
  static constexpr std::string_view xml{
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
)"};
};

class VmcInputs
{
public:
  enum class valid
  {
    CROWDS = 0,
    TINY,
    DEP_BATCH
  };
  static std::string_view getXml(valid val) { return xml[static_cast<std::size_t>(val)]; }
  static auto begin() { return xml.begin(); }
  static auto end() { return xml.end(); }

  // clang-format: off
  static constexpr std::array<std::string_view, 3> xml{
      R"(
  <qmc method="vmc" move="pbyp">
    <parameter name="crowds">                 8 </parameter>
    <estimators>
      <estimator name="LocalEnergy" hdf5="no" />
    </estimators>
    <parameter name="total_walkers">          32 </parameter>
    <parameter name="warmupSteps">            5 </parameter>
    <parameter name="substeps">               5 </parameter>
    <parameter name="steps">                  1 </parameter>
    <parameter name="blocks">                 2 </parameter>
    <parameter name="timestep">             1.0 </parameter>
    <parameter name="usedrift">              no </parameter>
  </qmc>
)",
      R"(
  <qmc method="vmc" move="pbyp">
    <parameter name="crowds">                 1 </parameter>
    <estimator name="LocalEnergy" hdf5="no" />
    <parameter name="walkers_per_rank">       1 </parameter>
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
    <parameter name="total_walkers">          32 </parameter>
    <parameter name="warmupSteps">            5 </parameter>
    <parameter name="substeps">               5 </parameter>
    <parameter name="steps">                  1 </parameter>
    <parameter name="blocks">                 2 </parameter>
    <parameter name="timestep">             1.0 </parameter>
    <parameter name="usedrift">              no </parameter>
  </qmc>
)"};
};

class DmcInputs
{
public:
  enum class valid
  {
    CROWDS = 0,
    DEP_BATCH,
    EMPERIOD
  };

  static std::string_view getXml(valid val) { return xml[static_cast<std::size_t>(val)]; }
  static auto begin() { return xml.begin(); }
  static auto end() { return xml.end(); }

private:
  static constexpr std::array<std::string_view, 3> xml{
      R"XML(
  <qmc method="dmc" move="pbyp">
    <parameter name="crowds">                 4 </parameter>
    <estimators>
      <estimator type="LocalEnergy" hdf5="no" />
    </estimators>
    <parameter name="total_walkers">          8 </parameter>
    <parameter name="reserve">             1.25 </parameter>
    <parameter name="warmupSteps">            5 </parameter>
    <parameter name="substeps">               5 </parameter>
    <parameter name="steps">                  1 </parameter>
    <parameter name="blocks">                 2 </parameter>
    <parameter name="timestep">             1.0 </parameter>
    <parameter name="usedrift">              no </parameter>
  </qmc>
)XML",
      R"XML(
  <qmc method="dmc_batch" move="pbyp">
<parameter name="crowds">                 4 </parameter>
    <estimator name="LocalEnergy" hdf5="no" />
    <parameter name="total_walkers">          8 </parameter>
    <parameter name="reserve">             1.25 </parameter>
    <parameter name="warmupSteps">            5 </parameter>
    <parameter name="substeps">               5 </parameter>
    <parameter name="steps">                  1 </parameter>
    <parameter name="blocks">                 2 </parameter>
    <parameter name="timestep">             1.0 </parameter>
    <parameter name="usedrift">              no </parameter>
  </qmc>
)XML",
      R"XML(
  <qmc method="dmc" move="pbyp">
    <estimators>
      <estimator type="PerParticleHamiltonianLogger" name="dmc_vem_test" to_stdout="false" />
    </estimators>
    <parameter name="crowds">                 4 </parameter>
    <estimator name="LocalEnergy" hdf5="no" />
    <parameter name="total_walkers">          16 </parameter>
    <parameter name="reserve">             1.25 </parameter>
    <parameter name="warmupSteps">            10 </parameter>
    <parameter name="substeps">               1 </parameter>
    <parameter name="steps">                  16 </parameter>
    <parameter name="blocks">                 2 </parameter>
    <parameter name="timestep">             1.0 </parameter>
    <parameter name="usedrift">              no </parameter>
    <parameter name="estimator_measurement_period">  4 </parameter>
  </qmc>
)XML"};
};


/** As far as I can tell these are no longer valid */
constexpr std::array<const char*, 2> valid_opt_input_sections{
    R"(
  <qmc method="opt" move="pbyp" gpu="yes">
  </qmc>
)",
    R"(
  <qmc method="opt" move="pbyp">
    <optimize method="test" output_param_file="yes"/>
    <parameter name="opt_num_crowds">                 4 </parameter>
    <parameter name="opt_crowd_size">                 8 </parameter>
  </qmc>
)"};

// clang-format: on
} // namespace testing
} // namespace qmcplusplus

#endif
