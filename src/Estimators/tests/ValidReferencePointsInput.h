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

#ifndef QMCPLUSPLUS_VALID_REFERENCEPOINTS_INPUT_H
#define QMCPLUSPLUS_VALID_REFERENCEPOINTS_INPUT_H

#include <array>
#include <string_view>

namespace qmcplusplus
{
namespace testing
{

struct ValidReferencePointsInputs
{
  static constexpr std::array<std::string_view, 2> xml{
      R"XML(
  <reference_points coord="cell">
    r1 1 0 0
    r2 0 1 0
    r3 0 0 1
  </reference_points>
)XML",
      R"XML(
  <reference_points coord="cartesian">
    r1 1 0 0
    r2 0 1 0
    r3 0 0 1
  </reference_points>
)XML"};

  enum valid
  {
    CELL = 0,
    CARTESIAN
  };
};

struct InvalidReferencePointsInputs
{
  static constexpr std::array<std::string_view, 2> xml{
      R"XML(
  <reference_points coord="cell">
    r1 1 0 0
    r2 0 1
    r3 0 0 1
  </reference_points>
)XML",
      R"XML(
  <reference_points coord="cartesian">
    r1 1 0 0
    r2 0 1 ab
  </reference_points>
)XML"};
  enum invalid
  {
    PARAM = 0,
    NUM
  };
};

} // namespace testing
} // namespace qmcplusplus

#endif
