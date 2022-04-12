//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_VALIDSPINDENSITYINPUT_H
#define QMCPLUSPLUS_VALIDSPINDENSITYINPUT_H

#include <array>

namespace qmcplusplus
{
namespace testing
{
// clang-format: off
  constexpr std::array<std::string_view, 3> valid_spin_density_input_sections{
    R"XML(
<estimator name="spindensity_new" type="spindensity" report="yes">
  <parameter name="grid">
    10 10 10
  </parameter>
  <parameter name="center">
    0.0 0.0 0.0
  </parameter>
  <parameter name="cell">
    3.37316115        3.37316115        0.00000000
    0.00000000        3.37316115        3.37316115
    3.37316115        0.00000000        3.37316115
  </parameter>
</estimator>
)XML",
    R"XML(
<estimator name="spindensity_new" type="spindensity" report="yes">
  <parameter name="dr">
    .4777 .4777 .4777
  </parameter>
  <parameter name="center">
    0.0 0.0 0.0
  </parameter>
  <parameter name="cell">
    3.37316115        3.37316115        0.00000000
    0.00000000        3.37316115        3.37316115
    3.37316115        0.00000000        3.37316115
  </parameter>
</estimator>
)XML",
    R"XML(
<estimator name="spindensity_new" type="spindensity" report="yes">
  <parameter name="dr">
    .4777 .4777 .4777
  </parameter>
  <parameter name="center">
    0.0 0.0 0.0
  </parameter>
</estimator>
)XML"};

// clang-format: on
constexpr int valid_spindensity_input_grid    = 0;
constexpr int valid_spindensity_input_dr      = 1;
constexpr int valid_spindensity_input_no_cell = 2;

} // namespace testing
} // namespace qmcplusplus

#endif /* QMCPLUSPLUS_VALIDSPINDENSITYINPUT_H */
