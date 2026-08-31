//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "MultiBspline.hpp"

namespace qmcplusplus
{
template<typename T>
MultiBspline<T>::~MultiBspline()
{
  for (auto spline_m : Base::spline_blocks)
  {
    coefs_allocator.deallocate(spline_m->coefs, spline_m->coefs_size);
    delete spline_m;
  }
}

template class MultiBspline<float>;
template class MultiBspline<double>;
} // namespace qmcplusplus
