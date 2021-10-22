//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_INVALID_OBDM_INPUT_H
#define QMCPLUSPLUS_INVALID_OBDM_INPUT_H

#include <array>

namespace qmcplusplus
{
namespace testing
{
// clang-format: off
constexpr std::array<const char*, 2> invalid_one_body_density_matrices_input_sections{
    R"(
<estimator type="dm1b" name="DensityMatrices">
  <parameter name="basis"        >  spo_u spo_uv  </parameter>
  <parameter name="evaluator"    >  matrix        </parameter>
  <parameter name="integrator"   >  path          </parameter>
  <parameter name="scale"        > -0.2           </parameter>
  <parameter name="samples"      >  64            </parameter>
  <parameter name="timestep"     >  0.5           </parameter>
  <parameter name="use_drift"    >  no            </parameter>
</estimator>
)",
    R"(
<estimator type="dm1b" name="DensityMatrices">
  <parameter name="basis"        >  dm_basis      </parameter>
  <parameter name="evaluator"    >  loop          </parameter>
  <parameter name="integrator"   >  uniform       </parameter>
  <parameter name="samples"      >  128           </parameter>
  <parameter name="scale"        >  1.1           </parameter>
  <parameter name="timestep"     >  0.5           </parameter>
  <parameter name="use_drift"    >  yes           </parameter>
</estimator>
)"
    // clang-format: on
};

constexpr int invalid_obdm_input_bad_integrator = 0;
constexpr int invalid_obdm_input_bad_scale      = 1;
} // namespace testing
} // namespace qmcplusplus

#endif
