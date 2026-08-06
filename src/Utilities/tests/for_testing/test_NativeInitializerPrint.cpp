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
#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <iostream>
#include <string>
#include <array>
#include <unordered_map>
#include <sstream>
#include <vector>
#include "NativeInitializerPrint.hpp"

namespace qmcplusplus
{
namespace
{
/** Helper function to process the NativePrint test output
 *
 * NativePrint wraps newline-delimited map entries in an outer pair of braces. Sort 
 * the entry lines to avoid relying on the unspecified iteration order of std::unordered_map.
 */
std::vector<std::string> getSortedMapEntries(const std::string& output)
{
  REQUIRE(output.size() >= 4); // Protect following checks from oob access
  REQUIRE(output.front() == '{');
  REQUIRE(output.compare(output.size() - 3, 3, "\n};") == 0);

  std::istringstream entry_stream(output.substr(1, output.size() - 4));
  std::vector<std::string> entries;
  for (std::string entry; std::getline(entry_stream, entry);)
    entries.push_back(entry);

  std::sort(entries.begin(), entries.end());
  return entries;
}
} // namespace

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
  std::unordered_map<std::string, std::vector<Vector<double>>> test_map_with_doubles{{"alpha", {{0.0, 1.1}, {2.2}}},
                                                                                     {"beta", {{3.3}, {4.4}}},
                                                                                     {"lambda",
                                                                                      {{6.5, 3.6}, {3.2, 3.3}}}};
  oss << NativePrint(test_map_with_doubles);
  const std::vector<std::string> expected_double_entries{{R"({{"alpha"}, {{0, 1.1, },{2.2, }, }},)"},
                                                         {R"({{"beta"}, {{3.3, },{4.4, }, }},)"},
                                                         {R"({{"lambda"}, {{6.5, 3.6, },{3.2, 3.3, }, }},)"}};
  CHECK(expected_double_entries == getSortedMapEntries(oss.str()));

  std::ostringstream oss2;
  std::unordered_map<std::string, std::vector<Vector<int>>> test_map_with_ints{{"alpha", {{0, 1}, {2}}},
                                                                               {"beta", {{3}, {4}}},
                                                                               {"lambda", {{6, 3}, {3, 3}}}};
  oss2 << NativePrint(test_map_with_ints);
  const std::vector<std::string> expected_int_entries{{R"({{"alpha"}, {{0, 1, },{2, }, }},)"},
                                                      {R"({{"beta"}, {{3, },{4, }, }},)"},
                                                      {R"({{"lambda"}, {{6, 3, },{3, 3, }, }},)"}};
  CHECK(expected_int_entries == getSortedMapEntries(oss2.str()));
}
} // namespace qmcplusplus
