//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_MULTIEINSPLINEOFFLOADMAPPER_HPP
#define QMCPLUSPLUS_MULTIEINSPLINEOFFLOADMAPPER_HPP

#include "MultiBsplineOffloadMapperBase.hpp"

namespace qmcplusplus
{
template<typename T>
class MultiBsplineOffloadMapper : public MultiBsplineOffloadMapperBase<T>
{
  using HostBspline = MultiBsplineBase<T>;
  using MultiBsplineOffloadMapperBase<T>::host_bsplines_;
  using MultiBsplineOffloadMapperBase<T>::block_coefs_;

  /// map host coefficients to devices
  void mapToDevice();

public:
  MultiBsplineOffloadMapper(const HostBspline& host_bsplines);

  ~MultiBsplineOffloadMapper();
};

extern template class MultiBsplineOffloadMapper<float>;
extern template class MultiBsplineOffloadMapper<double>;
} // namespace qmcplusplus

#endif
