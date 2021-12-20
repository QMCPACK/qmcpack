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

#ifndef QMCPLUSPLUS_MODERNSTRINGUTILS_HPP
#define QMCPLUSPLUS_MODERNSTRINGUTILS_HPP

#include <string>
#include <string_view>
#include <cctype>

namespace qmcplusplus
{

/** @ingroup C++17 string utility functions
 *  @{
 */
/** take string_view (or something that can be implicitly converted to one) and return lcase string.
 *  Don't call on anything but ASCII strings.
 *
 *  This would eventually replace the potentially unsafe to UTF-8 strtings tolower call in OhmmsElementBase.h
 *  which is used pretty randomly throughout the code.  The global function in Ohmmselement always operates inplace
 *  on the string buffer which is undesirable and should be explicitly evident in the code.
 *  i.e.
 *  my_string = strToLower(my_string);
 *
 *  According to en.cppreference.com std::tolower is only defined for usigned char and EOF. the std::tolower conversion
 *  is based on the _c_ locale which to makes it unclear on what might happen to char > 127. 
 *  *  For tags, and keywords where we define explicitly define they are ASCII encoded and this is fine. 
 *  *  For other XML derived text this should never be used since we should assume that to be UTF-8 encoded.
 */
inline std::string strToLower(const std::string_view s)
{
  std::string lower_str{s};
  std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return lower_str;
}
/** @} */


} // namespace qmcplusplus

#endif
