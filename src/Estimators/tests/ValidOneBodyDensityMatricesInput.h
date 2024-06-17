//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_VALID_OBDM_INPUT_H
#define QMCPLUSPLUS_VALID_OBDM_INPUT_H

#include <array>
#include <string_view>
#include "ValidInputHelp.h"

namespace qmcplusplus
{
namespace testing
{

class ValidOneBodyDensityMatricesInput
{
public:
  enum class valid : std::size_t
  {
    VANILLA = 0,
    SCALE,
    GRID
  };

  TEST_INPUT_ACCESSORS(valid)

private:
  static constexpr std::array<std::string_view, 3> xml{
      R"XML(
<estimator type="OneBodyDensityMatrices" name="OneBodyDensityMatrices">
  <parameter name="basis"        >  spo_ud spo_dm </parameter>
  <parameter name="evaluator"    >  matrix        </parameter>
  <parameter name="integrator"   >  density       </parameter>
  <parameter name="samples"      >  64            </parameter>
  <parameter name="timestep"     >  0.5           </parameter>
  <parameter name="center"       >  0.0 0.0 0.1   </parameter>
  <parameter name="use_drift"    >  yes           </parameter>
</estimator>
)XML",
      R"XML(
<estimator type="OneBodyDensityMatrices" name="OneBodyDensityMatrices">
  <parameter name="basis"        >  spo_ud spo_dm  </parameter>
  <parameter name="evaluator"    >  matrix         </parameter>
  <parameter name="integrator"   >  uniform       </parameter>
  <parameter name="samples"      >  128           </parameter>
  <parameter name="scale"        >  0.8           </parameter>
  <parameter name="timestep"     >  0.5           </parameter>
  <parameter name="use_drift"    >  yes           </parameter>
</estimator>
)XML",
      R"XML(
<estimator type="OneBodyDensityMatrices" name="OneBodyDensityMatrices">
  <parameter name="basis"        >  spo_ud spo_dm </parameter>
  <parameter name="evaluator"    >  matrix        </parameter>
  <parameter name="integrator"   >  uniform_grid  </parameter>
  <parameter name="points"       >  22            </parameter>
  <parameter name="scale"        >  0.8           </parameter>
  <parameter name="timestep"     >  0.5           </parameter>
  <parameter name="use_drift"    >  no            </parameter>
</estimator>
)XML"};

public:
};

} // namespace testing
} // namespace qmcplusplus

#endif
