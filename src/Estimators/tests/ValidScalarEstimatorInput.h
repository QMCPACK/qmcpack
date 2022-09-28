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

#ifndef QMCPLUSPLUS_VALIDSCALARESTIMATORINPUT_H
#define QMCPLUSPLUS_VALIDSCALARESTIMATORINPUT_H

#include <array>

namespace qmcplusplus
{
namespace testing
{

constexpr std::array<std::string_view, 5> valid_scalar_estimator_input_sections{
    R"XML(
   <estimator type="LocalEnergy" hdf5="no"/>
	  )XML",
    R"XML(
   <estimator type="CSLocalEnergy" npsi="2"/>
	  )XML",
    R"XML(
   <estimator name="CSLocalEnergy" npsi="2"/>
	  )XML",
    R"XML(
   <estimator name="eLocal" hdf5="yes"/>
	  )XML",
    R"XML(
   <estimator name="RMC" nObs="2"/>
	  )XML"};

constexpr int local_energy_input           = 0;
constexpr int cs_local_energy_input        = 1;
constexpr int cs_local_energy_input_legacy = 2;
constexpr int local_energy_input_legacy    = 3;
constexpr int rmc_local_energy_input       = 4;

} // namespace testing
} // namespace qmcplusplus

#endif /* QMCPLUSPLUS_VALIDSPINDENSITYINPUT_H */
