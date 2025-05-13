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
typename MultiBsplineBase<T>::SplineType* MultiBspline<T>::createImpl(const Ugrid grid[3],
                                                                      const typename Base::BoundaryCondition bc[3],
                                                                      int num_splines)
{
  static_assert(std::is_same<T, typename Alloc::value_type>::value, "MultiBspline and Alloc data types must agree!");
  if (getAlignedSize<T, Alloc::alignment>(num_splines) != num_splines)
    throw std::runtime_error("When creating the data space of MultiBspline, num_splines must be padded!\n");
  return myAllocator.allocateMultiBspline(grid[0], grid[1], grid[2], bc[0], bc[1], bc[2], num_splines);
}

template<typename T>
MultiBspline<T>::MultiBspline() = default;

template<typename T>
MultiBspline<T>::~MultiBspline()
{
  if (Base::spline_m != nullptr)
    myAllocator.destroy(Base::spline_m);
}

template class MultiBspline<float>;
template class MultiBspline<double>;
} // namespace qmcplusplus
