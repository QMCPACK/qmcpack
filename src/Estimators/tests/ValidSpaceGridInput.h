//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_VALID_SPACEGRID_INPUT_H
#define QMCPLUSPLUS_VALID_SPACEGRID_INPUT_H

#include <array>
#include <string_view>

namespace qmcplusplus
{
namespace testing
{

struct ValidSpaceGridInput
{
  enum valid
  {
    DEFAULT = 0,
    ORIGIN,
    CYLINDRICAL,
    SPHERICAL,
    WITH_STEP,
    WITHOUT_STEP,
    NUM_CASES
  };
  static constexpr std::array<std::string_view, NUM_CASES> xml{
      R"XML(
  <spacegrid coord="cartesian">
    <axis p1="a1" scale=".5" label="x" grid="-1 (.1) 1"/>
    <axis p1="a2" scale=".5" label="y" grid="-1 (.1) 1"/>
    <axis p1="a3" scale=".5" label="z" grid="-1 (.1) 1"/>
  </spacegrid>
)XML",
      R"XML(
  <spacegrid coord="cartesian">
     <origin p1="zero"/>
     <axis p1="a1" scale=".5" label="x" grid="-1 (.1) 1"/>
     <axis p1="a2" scale=".5" label="y" grid="-1 (.1) 1"/>
     <axis p1="a3" scale=".5" label="z" grid="-1 (.1) 1"/>
  </spacegrid>
)XML",
      R"XML(
  <spacegrid coord="cylindrical">
     <origin p1="zero"/>
     <axis p1="r1" scale=".5" label="r" grid="0 (.1) 1"/>
     <axis p1="r2" scale=".5" label="phi" grid="0 (.1) 1"/>
     <axis p1="r3" scale=".5" label="z" grid="-1 (.1) 1"/>
  </spacegrid>
)XML",
      R"XML(
  <spacegrid coord="spherical">
    <origin p1="ion1"/>
    <axis p1="r1" scale="6.9" label="r"     grid="0 (0.1) 1"/>
    <axis p1="r2" scale="6.9" label="phi"   grid="0 (0.1) 1"/>
    <axis p1="r3" scale="6.9" label="theta" grid="0 (0.1) 1"/>
  </spacegrid>
)XML",
      R"XML(
  <spacegrid coord="spherical">
    <origin p1="ion1"/>
    <axis p1="r1" scale="6.9" label="r" grid="0 (0.1) 1"/>
    <axis p1="r2" scale="6.9" label="phi" grid="0 (0.1) 1"/>
    <axis p1="r3" scale="6.9" label="theta" grid="0 (0.1) 1"/>
  </spacegrid>
)XML",
      R"XML(
  <spacegrid coord="spherical">
    <origin p1="ion2"/>
    <axis p1="r1" scale="6.9" label="r"     grid="0 1"/>
    <axis p1="r2" scale="6.9" label="phi"   grid="0 1"/>
    <axis p1="r3" scale="6.9" label="theta" grid="0 1"/>
  </spacegrid>
)XML"};
};

struct InvalidSpaceGridInput
{
  static constexpr std::array<std::string_view, 5> xml{
      R"XML(
  <spacegrid coord="cartesian">
    <axis p1="a1" scale=".5" label="x" grid="-1 (.1) 1"/>
    <axis p1="a2" scale=".5" label="q" grid="-1 (.1) 1"/>
    <axis p1="a3" scale=".5" label="z" grid="-1 (.1) 1"/>
  </spacegrid>
      )XML",
      R"XML(
  <spacegrid coord="sphericalt">
    <origin p1="ion1p"/>
    <axis p1="r6" scale="6.9" label="r"     grid="0 (0.1) 1"/>
    <axis p1="r2" scale="6.9" label="phi"   grid="0 (0.1) 1"/>
    <axis p1="r3" scale="6.9" label="theta" grid="0 (0.1) 1"/>
  </spacegrid>
      )XML",
      R"XML(
  <spacegrid coord="spherical">
    <origin p1="ion2"/>
    <axis p1="r1" scale="6.9" label="x"     grid="0 1"/>
    <axis p1="r2" scale="6.9" label="phi"   grid="0 1"/>
    <axis p1="r3" scale="6.9" label="theta" grid="0 1"/>
  </spacegrid>
      )XML",
      R"XML(
  <spacegrid coord="spherical">
    <origin p1="ion2"/>
    <axis p1="r2" scale="6.9" label="phi"   grid="0 1"/>
    <axis p1="r3" scale="6.9" label="theta" grid="0 1"/>
  </spacegrid>
      )XML",
      R"XML(
  <spacegrid coord="spherical">
    <origin p1="ion2"/>
    <axis p1="r6" scale="6.9" label="r"     grid="0 (0.1) 1"/>
    <axis p1="r2" scale="6.9" label="phi"   grid="1 -1"/>
    <axis p1="r3" scale="6.9" label="theta" grid="0 1"/>
  </spacegrid>
      )XML"};

  enum bad
  {
    LABELCART = 0,
    COORD,
    LABELSPHR,
    MISSINGAXIS,
    BADGRID
  };
};
} // namespace testing
} // namespace qmcplusplus

#endif
