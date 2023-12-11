//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File create by: Peter Doak
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_STLPRETTYPRINT_HPP
#define QMCPLUSPLUS_STLPRETTYPRINT_HPP

#include <iostream>
#include <vector>

namespace qmcplusplus
{

/** collapsed vector printout
 * [3, 3, 2, 2] is printed as [3(x2), 2(x2)]
 */
template<typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& rhs)
{
  out << "[";
  auto cursor = rhs.begin();
  while (cursor != rhs.end())
  {
    // each iteration handles one unique value
    const T ref_value = *cursor;
    size_t count      = 1;
    while (++cursor != rhs.end() && *cursor == ref_value)
      count++;
    out << ref_value;
    // identical elements are collapsed
    if (count > 1)
      out << "(x" << count << ")";
    // if not the last element, add a separator
    if (cursor != rhs.end())
      out << ", ";
  }
  out << "]";
  return out;
}

} // namespace qmcplusplus
#endif
