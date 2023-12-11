//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_VALIDMOMENTUMDISTRIBUTIONINPUT_H
#define QMCPLUSPLUS_VALIDMOMENTUMDISTRIBUTIONINPUT_H

#include <array>

namespace qmcplusplus
{
namespace testing
{
constexpr std::array<std::string_view, 2> valid_momentum_distribution_input_sections{
    R"(
<estimator type="MomentumDistribution" name="nofk" samples="5" kmax="3"/>
)",
    R"(
<estimator type="MomentumDistribution" name="nofk" samples="5" kmax0="3" kmax1="2" kmax2="1.5"/>
)"};

} // namespace testing
} // namespace qmcplusplus

#endif /* QMCPLUSPLUS_VALIDMOMENTUMDISTRIBUTIONINPUT_H */
