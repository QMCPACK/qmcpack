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

#ifndef QMCPLUSPLUS_SETUPDIRACDETERMINANTRESOURCES_HPP
#define QMCPLUSPLUS_SETUPDIRACDETERMINANTRESOURCES_HPP

#include "QMCWaveFunctions/Fermion/DiracDeterminantBatched.h"

namespace qmcplusplus
{
// checks if class has a guardMultiWalkerRes
template<class T, typename = decltype(void(std::declval<T>().guardMultiWalkerRes({})))>
std::true_type has_guard_aux(T);
std::false_type has_guard_aux(...);
template<class V>
struct has_guard : decltype(has_guard_aux(std::declval<V>()))
{};

namespace testing
{
struct SetupDiracDetResources
{
  template<class DiracDet,
           std::enable_if_t<
               std::is_same<DiracDet,
                            DiracDeterminantBatched<
                                MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>>::value,
               int> = 0>
  void operator()(DiracDet& ddet)
  {
    ddet.guardMultiWalkerRes();
  }

  template<class DiracDet,
           std::enable_if_t<
               !std::is_same<DiracDet,
                             DiracDeterminantBatched<
                                 MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>>::value,
               int> = 0>
  void operator()(DiracDet& ddet, double = 0)
  {}
};
} // namespace testing
} // namespace qmcplusplus
#endif
