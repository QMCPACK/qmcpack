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
template<typename T,
         typename COEFS_ALLOC         = aligned_allocator<T>,
         typename MULTI_SPLINE_ALLOC  = aligned_allocator<typename bspline_traits<T, 3>::SplineType>,
         typename SINGLE_SPLINE_ALLOC = aligned_allocator<typename bspline_traits<T, 3>::SingleSplineType>>
class MultiBspline
{
private:
  ///define the einsplie object type
  using SplineType = typename bspline_traits<T, 3>::SplineType;
  ///define the real type
  using real_type = typename bspline_traits<T, 3>::real_type;
  ///actual einspline multi-bspline object
  SplineType* spline_m;
  ///use allocator
  BsplineAllocator<T, COEFS_ALLOC, MULTI_SPLINE_ALLOC, SINGLE_SPLINE_ALLOC> myAllocator;

public:
  MultiBspline() : spline_m(nullptr) {}
  MultiBspline(const MultiBspline& in) = delete;
  MultiBspline& operator=(const MultiBspline& in) = delete;

  ~MultiBspline()
  {
    if (spline_m != nullptr)
      myAllocator.destroy(spline_m);
  }

  SplineType* getSplinePtr() { return spline_m; }

  /** create the einspline as used in the builder
   * @tparam GT grid type
   * @tparam BCT boundary type
   * @param bc num_splines number of splines
   *
   * num_splines must be padded to the aligned size. The caller must be aware of padding and pad all result arrays.
   */
  template<typename GT, typename BCT>
  void create(GT& grid, BCT& bc, int num_splines)
  {
    static_assert(std::is_same<T, typename COEFS_ALLOC::value_type>::value, "MultiBspline and ALLOC data types must agree!");
    if (getAlignedSize<T, COEFS_ALLOC::alignment>(num_splines) != num_splines)
      throw std::runtime_error("When creating the data space of MultiBspline, num_splines must be padded!\n");
    if (spline_m == nullptr)
    {
      typename bspline_traits<T, 3>::BCType xBC, yBC, zBC;
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

  /** copy a single spline to the big table
   * @tparam SingleSpline single spline type
   * @param aSpline UBspline_3d_(d,s)
   * @param int index of aSpline
   */
  template<typename SingleSpline>
  void copy_spline(SingleSpline* aSpline, int i)
  {
    if (spline_m == nullptr)
      throw std::runtime_error("The internal storage of MultiBspline must be created first!\n");
    if (aSpline->x_grid.num != spline_m->x_grid.num || aSpline->y_grid.num != spline_m->y_grid.num ||
        aSpline->z_grid.num != spline_m->z_grid.num)
      throw std::runtime_error("Cannot copy a single spline to MultiSpline with a different grid!\n");

    const int BaseOffset[3] = {0, 0, 0};
    const int BaseN[3]      = {spline_m->x_grid.num + 3, spline_m->y_grid.num + 3, spline_m->z_grid.num + 3};
    myAllocator.copy(aSpline, spline_m, i, BaseOffset, BaseN);
  }
};

} // namespace qmcplusplus

#endif
