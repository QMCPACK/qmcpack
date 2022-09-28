//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_TESTLISTENERFUNCTION_H
#define QMCPLUSPLUS_TESTLISTENERFUNCTION_H

namespace qmcplusplus
{
namespace testing
{

template<typename T>
auto getParticularListener(Matrix<T>& local_pots)
{
  return [&local_pots](const int walker_index, const std::string& name, const Vector<T>& inputV) {
    std::copy_n(inputV.begin(), inputV.size(), local_pots[walker_index]);
  };
}

} // namespace testing
} // namespace qmcplusplus

#endif
