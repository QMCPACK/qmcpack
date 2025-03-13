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
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file MultiBsplineBase.hpp
 *
 * define classes MultiBsplineBase
 * The evaluation functions are defined in MultiBsplineBaseEval.hpp
 */
#ifndef QMCPLUSPLUS_MULTIEINSPLINEBASE_HPP
#define QMCPLUSPLUS_MULTIEINSPLINEBASE_HPP

#include "spline2/bspline_traits.hpp"

namespace qmcplusplus
{
/** container class to hold a 3D multi spline pointer
 * @tparam T the precision of splines
 *
 * This class contains a pointer to a C object, copy and assign of this class is forbidden.
 */
template<typename T>
class MultiBsplineBase
{
protected:
  ///define the einsplie object type
  using SplineType = typename bspline_traits<T, 3>::SplineType;
  ///define the real type
  using real_type = typename bspline_traits<T, 3>::real_type;
  ///define the boundary condition type
  using BoundaryCondition = typename bspline_traits<T, 3>::BCType;
  ///actual einspline multi-bspline object
  SplineType* spline_m;

  virtual SplineType* createImpl(const Ugrid grid[3], const BoundaryCondition bc[3], int num_splines) = 0;

public:
  MultiBsplineBase() : spline_m(nullptr) {}
  MultiBsplineBase(const MultiBsplineBase& in)            = delete;
  MultiBsplineBase& operator=(const MultiBsplineBase& in) = delete;

  virtual ~MultiBsplineBase() = default;

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
    if (spline_m == nullptr)
    {
      BoundaryCondition xyzBC[3];
      xyzBC[0].lCode = bc[0].lCode;
      xyzBC[1].lCode = bc[1].lCode;
      xyzBC[2].lCode = bc[2].lCode;
      xyzBC[0].rCode = bc[0].rCode;
      xyzBC[1].rCode = bc[1].rCode;
      xyzBC[2].rCode = bc[2].rCode;
      xyzBC[0].lVal  = static_cast<T>(bc[0].lVal);
      xyzBC[1].lVal  = static_cast<T>(bc[1].lVal);
      xyzBC[2].lVal  = static_cast<T>(bc[2].lVal);
      xyzBC[0].rVal  = static_cast<T>(bc[0].rVal);
      xyzBC[1].rVal  = static_cast<T>(bc[1].rVal);
      xyzBC[2].rVal  = static_cast<T>(bc[2].rVal);
      spline_m       = createImpl(grid, xyzBC, num_splines);
    }
    else
      throw std::runtime_error("MultiBsplineBase::spline_m cannot be created twice!\n");
  }

  void flush_zero() const
  {
    if (spline_m != nullptr)
      std::fill(spline_m->coefs, spline_m->coefs + spline_m->coefs_size, T(0));
  }

  int num_splines() const { return (spline_m == nullptr) ? 0 : spline_m->num_splines; }

  size_t sizeInByte() const { return (spline_m == nullptr) ? 0 : spline_m->coefs_size * sizeof(T); }

  virtual void finalize(){};
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
  static_assert(std::is_floating_point<SSDT>::value, "SSDT must be a float point type");
  static_assert(std::is_floating_point<MSDT>::value, "MSDT must be a float point type");

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
