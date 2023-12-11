//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "ModernStringUtils.hpp"
#include <algorithm>
#include <cctype>

namespace qmcplusplus
{
using std::size_t;
using std::string_view;


std::string lowerCase(const std::string_view s)
{
  std::string lower_str{s};
  std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return lower_str;
}

namespace modernstrutil
{
std::vector<std::string_view> split(const string_view s, const string_view delimiters)
{
  std::vector<std::string_view> tokens;
  size_t right = 0;
  size_t left  = 0;
  while (true)
  {
    left = s.find_first_not_of(delimiters, right);
    if (left == s.npos)
      break;
    else
      right = s.find_first_of(delimiters, left);
    if (right == s.npos)
      right = s.size();
    size_t count = right - left;
    tokens.push_back(s.substr(left, count));
  }
  return tokens;
}

std::string_view strip(const string_view s)
{
  std::string_view delimiters = " \t\n\0";
  size_t left                 = s.find_first_not_of(delimiters, 0);
  if (left != s.npos) {
    size_t right = s.find_last_not_of(delimiters, s.npos);
    return s.substr(left, right - left + 1);
  }
  else
    return std::string_view{};
}
} // namespace modernstrutil


} // namespace qmcplusplus
