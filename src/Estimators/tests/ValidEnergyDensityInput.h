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

#ifndef QMCPLUSPLUS_VALID_ENERGYDENSITY_INPUT_H
#define QMCPLUSPLUS_VALID_ENERGYDENSITY_INPUT_H

#include <array>
#include <string_view>

namespace qmcplusplus
{
namespace testing
{

class EnergyDensityInputs
{
public:
  enum class valid
  {
    CELL = 0,
    ION
  };

  static std::string_view getXml(valid val) { return xml[static_cast<std::size_t>(val)]; }
  auto begin() { return xml.begin(); }
  auto end() {return xml.end(); }

private:
  static constexpr std::array<std::string_view, 2> xml{
      R"XML(
<estimator type="EnergyDensity" name="EDcell" dynamic="e" static="ion">
   <spacegrid coord="cartesian">
     <origin p1="zero"/>
     <axis p1="a1" scale=".5" label="x" grid="-1 (.05) 1"/>
     <axis p1="a2" scale=".5" label="y" grid="-1 (.1) 1"/>
     <axis p1="a3" scale=".5" label="z" grid="-1 (.1) 1"/>
   </spacegrid>
</estimator>
      )XML",
      R"XML(
<estimator type="EnergyDensity" name="EDatom" dynamic="e" static="ion">
  <reference_points coord="cartesian">
    r1 1 0 0
    r2 0 1 0
    r3 0 0 1
  </reference_points>
  <spacegrid coord="spherical">
    <origin p1="ion1"/>
    <axis p1="r1" scale="6.9" label="r"     grid="0 1"/>
    <axis p1="r2" scale="6.9" label="phi"   grid="0 1"/>
    <axis p1="r3" scale="6.9" label="theta" grid="0 1"/>
  </spacegrid>
  <spacegrid coord="spherical">
    <origin p1="ion2"/>
    <axis p1="r1" scale="6.9" label="r"     grid="0 1"/>
    <axis p1="r2" scale="6.9" label="phi"   grid="0 1"/>
    <axis p1="r3" scale="6.9" label="theta" grid="0 1"/>
  </spacegrid>
</estimator>
)XML"};
};

class InvalidEnergyDensityInput
{
public:
  enum class invalid
  {
    BADNAME = 0,
  };

  static std::string_view getXml(invalid val) { return xml[static_cast<std::size_t>(val)]; }
  auto begin() { return xml.begin(); }
  auto end() {return xml.end(); }

private:
  static constexpr std::array<std::string_view, 1> xml{
      R"XML(
<estimator type="EnergyDensity" name="EDcell" dynamic="e" static="ion0">
</estimator>
      )XML"};
};

} // namespace testing
} // namespace qmcplusplus

#endif
