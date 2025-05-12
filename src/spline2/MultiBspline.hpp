//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file MultiBspline.hpp
 *
 * define classes MultiBspline
 * The evaluation functions are defined in MultiBsplineEval.hpp
 */
#ifndef QMCPLUSPLUS_MULTIEINSPLINE_HPP
#define QMCPLUSPLUS_MULTIEINSPLINE_HPP

#include "MultiBsplineBase.hpp"
#include "spline2/BsplineAllocator.hpp"
#include "CPU/SIMD/aligned_allocator.hpp"

namespace qmcplusplus
{
/** container class to hold a 3D multi spline pointer and BsplineAllocator
 * @tparam T the precision of splines
 */
template<typename T>
class MultiBspline : public MultiBsplineBase<T>
{
private:
  using Base  = MultiBsplineBase<T>;
  using Alloc = aligned_allocator<T>;
  ///use allocator
  BsplineAllocator<T, Alloc> myAllocator;

  typename Base::SplineType* createImpl(const Ugrid grid[3],
                                        const typename Base::BoundaryCondition bc[3],
                                        int num_splines) override;

public:
  MultiBspline();
  ~MultiBspline() override;
};

extern template class MultiBspline<float>;
extern template class MultiBspline<double>;
} // namespace qmcplusplus

#endif
