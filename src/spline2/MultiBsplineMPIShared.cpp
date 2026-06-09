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


#include "MultiBsplineMPIShared.hpp"

namespace qmcplusplus
{
template<typename T>
MultiBsplineMPIShared<T>::~MultiBsplineMPIShared()
{
  MPI_Win_free(&win);

  for (auto spline_m : Base::spline_blocks)
    delete spline_m;
}

template class MultiBsplineMPIShared<float>;
template class MultiBsplineMPIShared<double>;
} // namespace qmcplusplus
