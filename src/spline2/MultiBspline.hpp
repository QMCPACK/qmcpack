//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
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
#ifndef QMCPLUSPLUS_MULTIEINSPLINE_COMMON_HPP
#define QMCPLUSPLUS_MULTIEINSPLINE_COMMON_HPP
#include <iostream>
#include <cstdlib>
#include <type_traits>
#include "config.h"
#include "spline2/BsplineAllocator.hpp"

namespace qmcplusplus
{
/** container class to hold a 3D multi spline pointer and BsplineAllocator
 * @tparam T the precision of splines
 * @tparam ALLOC memory allocator
 *
 * This class contains a pointer to a C object, copy and assign of this class is forbidden.
 */
template<typename T, typename ALLOC = aligned_allocator<T>>
class MultiBspline
{
private:
  ///define the einsplie object type
  using SplineType = typename bspline_traits<T, 3>::SplineType;
  ///define the real type
  using real_type = typename bspline_traits<T, 3>::real_type;
  ///define the boundary condition type
  using BoundaryCondition = typename bspline_traits<T, 3>::BCType;
  ///actual einspline multi-bspline object
  SplineType* spline_m;
  ///use allocator
  BsplineAllocator<T, ALLOC> myAllocator;

public:
  MultiBspline() : spline_m(nullptr) {}
  MultiBspline(const MultiBspline& in)            = delete;
  MultiBspline& operator=(const MultiBspline& in) = delete;

  ~MultiBspline()
  {
    if (spline_m != nullptr)
      myAllocator.destroy(spline_m);
  }

  SplineType* getSplinePtr() { return spline_m; }

  /** create the einspline as used in the builder
   * @tparam BCT boundary type
   * @param bc num_splines number of splines
   *
   * num_splines must be padded to the aligned size. The caller must be aware of padding and pad all result arrays.
   */
  template<typename BCT>
  inline void create(const Ugrid grid[3], const BCT& bc, int num_splines)
  {
    static_assert(std::is_same<T, typename ALLOC::value_type>::value, "MultiBspline and ALLOC data types must agree!");
    if (getAlignedSize<T, ALLOC::alignment>(num_splines) != num_splines)
      throw std::runtime_error("When creating the data space of MultiBspline, num_splines must be padded!\n");
    if (spline_m == nullptr)
    {
      BoundaryCondition xBC, yBC, zBC;
      xBC.lCode = bc[0].lCode;
      yBC.lCode = bc[1].lCode;
      zBC.lCode = bc[2].lCode;
      xBC.rCode = bc[0].rCode;
      yBC.rCode = bc[1].rCode;
      zBC.rCode = bc[2].rCode;
      xBC.lVal  = static_cast<T>(bc[0].lVal);
      yBC.lVal  = static_cast<T>(bc[1].lVal);
      zBC.lVal  = static_cast<T>(bc[2].lVal);
      xBC.rVal  = static_cast<T>(bc[0].rVal);
      yBC.rVal  = static_cast<T>(bc[1].rVal);
      zBC.rVal  = static_cast<T>(bc[2].rVal);
      spline_m  = myAllocator.allocateMultiBspline(grid[0], grid[1], grid[2], xBC, yBC, zBC, num_splines);
    }
    else
      throw std::runtime_error("MultiBspline::spline_m cannot be created twice!\n");
  }

  void flush_zero() const
  {
    if (spline_m != nullptr)
      std::fill(spline_m->coefs, spline_m->coefs + spline_m->coefs_size, T(0));
  }

  int num_splines() const { return (spline_m == nullptr) ? 0 : spline_m->num_splines; }

  size_t sizeInByte() const { return (spline_m == nullptr) ? 0 : spline_m->coefs_size * sizeof(T); }

  template<typename Allocator = ALLOC, typename = qmcplusplus::IsDualSpace<Allocator>>
  void finalize()
  {
    myAllocator.finalize(spline_m);
  }
};


/** copy a single spline to multi spline
   * @tparam SSDT single spline coefficient data type
   * @tparam MSDT multi spline coefficient data type
   * @param single UBspline_3d_(d,s)
   * @param multi multi_UBspline_3d_(d,s)
   * @param int index of single in multi
   */
template<typename SSDT, typename MSDT>
void copy_spline(const typename bspline_traits<SSDT, 3>::SingleSplineType& single,
                 typename bspline_traits<MSDT, 3>::SplineType& multi,
                 int i)
{
  if (single.x_grid.num != multi.x_grid.num || single.y_grid.num != multi.y_grid.num ||
      single.z_grid.num != multi.z_grid.num)
    throw std::runtime_error("Cannot copy a single spline to MultiSpline with a different grid!\n");

  intptr_t x_stride_in  = single.x_stride;
  intptr_t y_stride_in  = single.y_stride;
  intptr_t x_stride_out = multi.x_stride;
  intptr_t y_stride_out = multi.y_stride;
  intptr_t z_stride_out = multi.z_stride;
  const intptr_t istart = static_cast<intptr_t>(i);
  const intptr_t n0 = multi.x_grid.num + 3, n1 = multi.y_grid.num + 3, n2 = multi.z_grid.num + 3;
  for (intptr_t ix = 0; ix < n0; ++ix)
    for (intptr_t iy = 0; iy < n1; ++iy)
    {
      auto* __restrict__ out      = multi.coefs + ix * x_stride_out + iy * y_stride_out + istart;
      const auto* __restrict__ in = single.coefs + ix * x_stride_in + iy * y_stride_in;
      for (intptr_t iz = 0; iz < n2; ++iz)
        out[iz * z_stride_out] = in[iz];
    }
}

} // namespace qmcplusplus

#endif
