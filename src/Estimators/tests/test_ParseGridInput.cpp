//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Peter  W. Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter  W. Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

/** \file
 *  This unit test was added long after the initial implementation of this grid spec
 *  the "correct" values are produced from the ported code as the original code is
 *  design does not easy allow a unit test.  It may not cover edge case behavior of the
 *  original code in QMCHamiltonian/SpaceGrid.cpp
 */
#include "catch.hpp"

#include "ParseGridInput.hpp"
#include "Message/UniformCommunicateError.h"
#include <StlPrettyPrint.hpp>
#include <NativeInitializerPrint.hpp>

namespace qmcplusplus
{

template<typename T>
std::ostream& operator<<(std::ostream& out, const AxisGrid<T>& rhs)
{
  out << "{";
  out << NativePrint(rhs.ndom_int) << ", //ndom_int\n";
  out << NativePrint(rhs.ndu_int) << ", //ndu_int\n";
  out << NativePrint(rhs.du_int) << ", //du_int\n";
  out << rhs.umin << ", //umin\n";
  out << rhs.umax << ", //umax\n";
  out << rhs.odu << ", //odu\n";
  out << NativePrint(rhs.gmap) << ", //gmap\n";
  out << NativePrint(rhs.ndu_per_interval) << ", //ndu_per_interval\n";
  out << rhs.dimensions << "} //dimensions\n";
  return out;
}

TEMPLATE_TEST_CASE("ParseGridInput::Good", "[estimators]", float, double)
{
  using Real = TestType;
  using namespace std::string_literals;

  std::istringstream grid_input_1("0 (0.1) 1"s);
  AxisGrid<Real> axis_grid_1 = parseGridInput<Real>(grid_input_1);
  AxisGrid<Real> expected_1{{
                                10,
                            }, //ndom_int
                            {
                                1,
                            }, //ndu_int
                            {
                                0.1,
                            },  //du_int
                            0,  //umin
                            1,  //umax
                            10, //odu
                            {
                                0,
                                1,
                                2,
                                3,
                                4,
                                5,
                                6,
                                7,
                                8,
                                9,
                            }, //gmap
                            {
                                1,
                                1,
                                1,
                                1,
                                1,
                                1,
                                1,
                                1,
                                1,
                                1,
                            },   //ndu_per_interval
                            10}; //dimensions

  CHECK(axis_grid_1 == expected_1);

  std::istringstream grid_input_2("0.1 (10) 0.2 (20) 0.4 (10) 0.8"s);
  AxisGrid<Real> axis_grid_2 = parseGridInput<Real>(grid_input_2);
  AxisGrid<Real> expected_2{{
                                10,
                                20,
                                10,
                            }, //ndom_int
                            {
                                1,
                                1,
                                4,
                            }, //ndu_int
                            {
                                0.01,
                                0.01,
                                0.04,
                            },   //du_int
                            0.1, //umin
                            0.8, //umax
                            100, //odu
                            {
                                0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17,
                                18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 30, 30, 30, 31, 31,
                                31, 31, 32, 32, 32, 32, 33, 33, 33, 33, 34, 34, 34, 34, 35, 35, 35, 35,
                                36, 36, 36, 36, 37, 37, 37, 37, 38, 38, 38, 38, 39, 39, 39, 39,
                            }, //gmap
                            {
                                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                            },   //ndu_per_interval
                            40}; //dimensions

  CHECK(axis_grid_2 == expected_2);
}

TEMPLATE_TEST_CASE("ParseGridInput::Bad", "[estimators]", float, double)
{
  using Real = TestType;
  using namespace std::string_literals;

  std::istringstream grid_input_1("1.1 (0.1) 1.5"s);
  CHECK_THROWS_AS(parseGridInput<Real>(grid_input_1), UniformCommunicateError);

  std::istringstream grid_input_2("0.8 (0.1) 0.5"s);
  CHECK_THROWS_AS(parseGridInput<Real>(grid_input_2), UniformCommunicateError);

  std::istringstream grid_input_3("0.8 (0.1) 1.5"s);
  CHECK_THROWS_AS(parseGridInput<Real>(grid_input_3), UniformCommunicateError);

  std::istringstream grid_input_4("0.1 (0.3) 0.2"s);
  CHECK_THROWS_AS(parseGridInput<Real>(grid_input_4), UniformCommunicateError);

  std::istringstream grid_input_5("0.1 (10) 0.2 (20) 0.35 (10) 0.73"s);
  CHECK_THROWS_AS(parseGridInput<Real>(grid_input_5), UniformCommunicateError);
}

TEMPLATE_TEST_CASE("ParseGridInput_constructors", "[estimators]", float, double)
{
  using Real = TestType;
  AxisGrid<Real> start_grid{{
                                10,
                            }, //ndom_int
                            {
                                1,
                            }, //ndu_int
                            {
                                0.1,
                            },  //du_int
                            0,  //umin
                            1,  //umax
                            10, //odu
                            {
                                0,
                                1,
                                2,
                                3,
                                4,
                                5,
                                6,
                                7,
                                8,
                                9,
                            }, //gmap
                            {
                                1,
                                1,
                                1,
                                1,
                                1,
                                1,
                                1,
                                1,
                                1,
                                1,
                            },   //ndu_per_interval
                            10}; //dimensions
  AxisGrid<Real> copied_grid(start_grid);
  AxisGrid<Real> assigned_grid(start_grid);

  CHECK(start_grid == copied_grid);
  CHECK(start_grid == assigned_grid);
}

} // namespace qmcplusplus
