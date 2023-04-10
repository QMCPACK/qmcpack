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

#ifndef QMCPLUSPLUS_VALID_MAGNETIZATIONDENSITY_INPUT_H
#define QMCPLUSPLUS_VALID_MAGNETIZATIONDENSITY_INPUT_H

#include <array>
#include <string_view>

namespace qmcplusplus
{
namespace testing
{
namespace magdensity
{
enum Inputs
{
  valid_magdensity_input = 0,
  valid_magdensity_input_dr,
  valid_magdensity_input_grid,
  valid_magdensity_input_unittest
};

// clang-format: off
constexpr std::array<std::string_view, 4> valid_mag_density_input_sections{
    R"XML(
<estimator type="MagnetizationDensity" name="magdensity">
  <parameter name="integrator"   >  simpsons       </parameter>
  <parameter name="samples"      >  64             </parameter>
  <parameter name="center"       >  0.0 0.0 0.1    </parameter>
  <parameter name="grid"         >  4 3 2          </parameter>
</estimator>
)XML",
    R"XML(
<estimator type="MagnetizationDensity" name="magdensity">
  <parameter name="integrator"   >  montecarlo       </parameter>
  <parameter name="samples"      >  128              </parameter>
  <parameter name="dr"           >  0.01 0.02 0.03   </parameter>
</estimator>
)XML",
    R"XML(
<estimator type="MagnetizationDensity" name="magdensity">
  <parameter name="samples"      >  32            </parameter>
  <parameter name="corner"       >  0.0 0.0 0.1   </parameter>
  <parameter name="grid"         >  4 3 2         </parameter>
</estimator>
)XML",
    R"XML(
<estimator type="MagnetizationDensity" name="magdensity">
  <parameter name="integrator"   >  simpsons       </parameter>
  <parameter name="samples"      >  9            </parameter>
  <parameter name="corner"       >  0.0 0.0 0.0   </parameter>
  <parameter name="grid"         >  2 2 2         </parameter>
</estimator>
)XML"
    // clang-format: on
};
} // namespace magdensity

} // namespace testing
} // namespace qmcplusplus

#endif
