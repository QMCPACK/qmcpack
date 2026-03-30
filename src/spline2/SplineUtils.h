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


#ifndef QMCPLUSPLUS_SPLINE_UTILS_H
#define QMCPLUSPLUS_SPLINE_UTILS_H

#include "hdf/hdf_archive.h"
#include "spline2/MultiBsplineBase.hpp"
#include "spline2/MultiBspline1D.hpp"
#include "Message/Communicate.h"

namespace qmcplusplus
{
template<typename ST>
class SplineUtils
{
  static constexpr size_t my_index = 0;

public:
  static bool read(MultiBsplineBase<ST>& spline, hdf_archive& h5f);
  static bool write(MultiBsplineBase<ST>& spline, hdf_archive& h5f);

  static bool read(MultiBspline1D<ST>& spline, hdf_archive& h5f);
  static bool write(MultiBspline1D<ST>& spline, hdf_archive& h5f);

  static void gatherv(MultiBsplineBase<ST>& spline, const std::vector<int>& offset, Communicate& comm);
  static void bcast(MultiBsplineBase<ST>& spline, Communicate& comm);

  static void gatherv(MultiBspline1D<ST>& spline, size_t stride, const std::vector<int>& offset, Communicate& comm);
  static void bcast(MultiBspline1D<ST>& spline, Communicate& comm);
};

extern template class SplineUtils<float>;
extern template class SplineUtils<double>;


} // namespace qmcplusplus
#endif
