//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include <iostream>
#include <string>
#include <array>
#include <unordered_map>
#include <strstream>
#include "NativeInitializerPrint.hpp"

/** \file
 *  Tests NativeInitializerPrint more strictly than necessary, what it really needs to produce
 *  is a text representation that can be pasted into test programs and produce the same
 *  Native representation as the object that was written out. That would require a second pass
 *  of compilation as part of the test program which I'd be impressed if someone could
 *  figure out how to integrate with our unite testing.
 */
namespace qmcplusplus
{

TEST_CASE("NativePrint::array", "[utilities][for_testing]")
{
  std::ostringstream oss;
  std::array<double, 4> test_array{1.1, 2.2, 3.3, 4.4};
  oss << NativePrint(test_array);
  CHECK(std::string("{ 1.1, 2.2, 3.3, 4.4,  }") == oss.str());
  std::array<double, 4> test_array2{1.11111111111111, 2.2222222222222, 3.333333333333333, 4.44444444444444};
  std::ostringstream oss2;
  oss2 << NativePrint(test_array2);
  CHECK(std::string("{ 1.111111111, 2.222222222, 3.333333333, 4.444444444,  }") == oss2.str());
}

TEST_CASE("NativePrint::unordered_map<std::string, std::vector<Vector<T>>>", "[utilities][for_testing]")
{
  std::ostringstream oss;
  std::unordered_map<std::string, std::vector<Vector<double>>> test_map{{"alpha", {{0.0, 1.1}, {2.2}}},
                                                                        {"beta", {{3.3}, {4.4}}},
                                                                        {"lambda", {{6.5, 3.6}, {3.2, 3.3}}}};
  oss << NativePrint(test_map);
  CHECK(std::string("{{{\"lambda\"}, {{6.5, 3.6, },{3.2, 3.3, }, }},\n{{\"beta\"}, {{3.3, },{4.4, }, "
                    "}},\n{{\"alpha\"}, {{0, 1.1, },{2.2, }, }},\n};") == oss.str());

  std::ostringstream oss2;
  std::unordered_map<std::string, std::vector<Vector<int>>> test_map2{{"alpha", {{0, 1}, {2}}},
                                                                      {"beta", {{3}, {4}}},
                                                                      {"lambda", {{6, 3}, {3, 3}}}};
  oss2 << NativePrint(test_map2);
  CHECK(std::string("{{{\"lambda\"}, {{6, 3, },{3, 3, }, }},\n{{\"beta\"}, {{3, },{4, }, "
                    "}},\n{{\"alpha\"}, {{0, 1, },{2, }, }},\n};") == oss2.str());
}
} // namespace qmcplusplus
