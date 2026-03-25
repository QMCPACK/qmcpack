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
  Alloc coefs_allocator;

public:
  template<typename BCT>
  MultiBspline(const Ugrid grid[3], const BCT& bc, size_t num_splines) : Base({0, num_splines})
  {
    Base::spline_blocks.resize(1, nullptr);
    auto* spline_m         = new typename Base::SplineType;
    Base::spline_blocks[0] = spline_m;
    Base::setMetaData(*spline_m, grid[0], grid[1], grid[2], Base::createBoundaryCondition(bc).data(), num_splines,
                      getAlignedSize<T, Alloc::alignment>(num_splines));
    spline_m->coefs = coefs_allocator.allocate(spline_m->coefs_size);
  }

  ~MultiBspline() override;
};

extern template class MultiBspline<float>;
extern template class MultiBspline<double>;
} // namespace qmcplusplus

#endif
