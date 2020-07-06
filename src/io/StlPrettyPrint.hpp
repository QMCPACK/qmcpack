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
  copy(rhs.begin(), rhs.end(), std::ostream_iterator<double>(out, ", "));
  out << "]";
  out <<'\n';
  return out;
}

}
#endif
