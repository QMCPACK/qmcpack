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

#ifndef QMCPLUSPLUS_MODERNSTRINGUTILS_HPP
#define QMCPLUSPLUS_MODERNSTRINGUTILS_HPP

#include <charconv>
#include <stdexcept>
#include <string>
#include <string_view>
#include <sstream>
#include <type_traits>
#include <vector>

namespace qmcplusplus
{

/** @ingroup C++17 string utility functions
 *  @{
 */
/** take string_view (or something that can be implicitly converted to one) and return lcase string.
 *  Don't call on anything but ASCII strings.
 *
 *  According to en.cppreference.com std::tolower is only defined for usigned char and EOF. the std::tolower conversion
 *  is based on the _c_ locale which to makes it unclear on what might happen to char > 127. 
 *  *  For tags, and keywords where we define explicitly define they are ASCII encoded and this is fine. 
 *  *  For other XML derived text this should never be used since we should assume that to be UTF-8 encoded.
 */
std::string lowerCase(const std::string_view s);

/** prevent clash with string_utils.h */
namespace modernstrutil
{
/** return string_view tokens */
std::vector<std::string_view> split(const std::string_view s, const std::string_view delimiters);
/** remove white space from each beginning and end of string_view */
std::string_view strip(const std::string_view s);
} // namespace modernstrutil

/** alternate to string2real
 *  calls c++ string to real conversion based on T's precision template<typename T>
 */
template<typename T>
inline T string2Real(const std::string_view svalue)
{
  static_assert(std::is_floating_point_v<T>);
  T result;
// full support for floating point from char not present until stdlibc++ aligned with gcc11
// there is still not from_char for floats as of libc++ 15
#if _GLIBCXX_RELEASE > 10
  auto [prt, ec] = std::from_chars(svalue.data(), svalue.data() + svalue.size(), result, std::chars_format::general);
  if (ec != std::errc())
    throw std::runtime_error("Could not convert from string to real value");
#else
  // atof must be given a null terminated string, string_view is not guaranteed to have a null terminator.
  std::string str_value(svalue);
  result = static_cast<T>(atof(str_value.c_str()));
#endif
  return result;
}

/** alternate to string2real
 *  calls c++ string to real conversion based on T's precision template<typename T>
 */
template<typename T>
inline T string2Int(const std::string_view svalue)
{
  static_assert(std::is_integral_v<T>);
  // full support for floating point from char not present until stdlibc++ aligned with gcc11
  // there is still not from_char for floats as of libc++ 15
  T result{};
#if _GLIBCXX_RELEASE > 10
  auto [prt, ec] = std::from_chars(svalue.data(), svalue.data() + svalue.size(), result);
  // std::errc() represents success or a value of 0 with respect to errc's error enumeration.
  if (ec != std::errc())
  {
    if (ec == std::errc::result_out_of_range)
      throw std::range_error("Value out of range for integral type parameter of string2Int");
    else
    {
      std::ostringstream msg;
      msg << "Could not convert from string " << std::string(svalue) << " to int!";
      throw std::runtime_error(msg.str());
    }
  }
#else
  // atof must be given a null terminated string, string_view is not guaranteed to have a null terminator.
  std::string str_value(svalue);
  if constexpr (std::is_same_v<T, int>)
    result = atoi(str_value.c_str());
  else if constexpr (std::is_same_v<T, long>)
    result = atol(str_value.c_str());
  else if constexpr (std::is_same_v<T, long long>)
    result = atol(str_value.c_str());
  else
    throw std::runtime_error("unsupported type for string to integral type conversion with pre v.10 stdlibc++!");
#endif
  return result;
}
/** @} */

} // namespace qmcplusplus

#endif
