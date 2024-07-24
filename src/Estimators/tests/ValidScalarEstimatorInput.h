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

#ifndef QMCPLUSPLUS_VALIDSCALARESTIMATORINPUT_H
#define QMCPLUSPLUS_VALIDSCALARESTIMATORINPUT_H

#include <array>
#include <string_view>
#include "ValidInputHelp.h"

namespace qmcplusplus
{
namespace testing
{

class ScalarEstimatorInputs
{
public:
  enum class valid : std::size_t
  {
    LOCAL_ENERGY = 0,
    CS_LOCAL_ENERGY,
    LOCAL_ENERGY_LEGACY,
    RMC_LOCAL_ENERGY
  };

  TEST_INPUT_ACCESSORS(valid)

private:
  static constexpr std::array<std::string_view, 4> xml{
      R"XML(
   <estimator type="LocalEnergy" hdf5="no"/>
	  )XML",
      R"XML(
   <estimator type="CSLocalEnergy" npsi="2"/>
	  )XML",
      R"XML(
   <estimator name="eLocal" hdf5="yes"/>
	  )XML",
      R"XML(
   <estimator name="RMC" nObs="2"/>
	  )XML"};
};

} // namespace testing
} // namespace qmcplusplus

#endif /* QMCPLUSPLUS_VALIDSPINDENSITYINPUT_H */
