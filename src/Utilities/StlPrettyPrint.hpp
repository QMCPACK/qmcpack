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

template<typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& rhs)
{
  out << "[";
  if (!rhs.empty())
  {
    auto last = rhs.end();
    last--;
    copy(rhs.begin(), last, std::ostream_iterator<T>(out, ", "));
    out << *last;
  }
  out << "]";
  return out;
}

}
#endif
