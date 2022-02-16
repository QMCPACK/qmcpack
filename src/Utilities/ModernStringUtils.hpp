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
/** @} */


} // namespace qmcplusplus

#endif
