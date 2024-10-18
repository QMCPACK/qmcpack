//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_INVALID_MAGNETIZATIONDENSITY_INPUT_H
#define QMCPLUSPLUS_INVALID_MAGNETIZATIONDENSITY_INPUT_H

#include <array>
#include <string_view>

namespace qmcplusplus
{
namespace testing
{
// clang-format: off
constexpr std::array<std::string_view, 3> invalid_mag_density_input_sections{
    R"XML(
<estimator type="MagnetizationDensity" name="magdensity">
  <parameter name="integrator"   >  tacocat       </parameter>
  <parameter name="samples"      >  64            </parameter>
  <parameter name="center"       >  0.0 0.0 0.1   </parameter>
  <parameter name="grid"         >  4 4 4         </parameter>
</estimator>
)XML",
    R"XML(
<estimator type="MagnetizationDensity" name="magdensity">
  <parameter name="center"       >  0.0 1.0 0.0 </parameter>
  <parameter name="corner"       >  0.0 0.1 0.0 </parameter>
  <parameter name="samples"      >  128         </parameter>
  <parameter name="dr"           >  0.1 0.1 0.1 </parameter>
</estimator>
)XML",
    R"XML(
<estimator type="MagnetizationDensity" name="magdensity">
  <parameter name="grid"       >  0.9 1 5   </parameter>
</estimator>
)XML"
    // clang-format: on
};

constexpr int invalid_magnetization_density_integrator   = 0;
constexpr int invalid_magnetization_density_cornercenter = 1;
constexpr int invalid_magnetization_density_badgrid      = 2;
} // namespace testing
} // namespace qmcplusplus

#endif
