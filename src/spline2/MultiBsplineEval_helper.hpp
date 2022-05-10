////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file MultiBsplineEval_helper.hpp
 *
 * header file for helper functions
 */
#ifndef SPLINE2_MULTIEINSPLINE_EVAL_HELPER_HPP
#define SPLINE2_MULTIEINSPLINE_EVAL_HELPER_HPP

#include <cmath>
#include <algorithm>
#include "bspline_traits.hpp"
#include "spline2/MultiBsplineData.hpp"

namespace spline2
{

/** break x into the integer part and residual part and apply bounds
 * @param x input coordinate
 * @param dx fractional part
 * @param ind integer part
 * @param nmax upper bound of the integer part
 *
 * x in the range of [0, nmax+1) will be split correctly.
 * x < 0, ind = 0, dx = 0
 * x >= nmax+1, ind = nmax, dx = 1 - epsilon
 *
 * Attention: nmax is not the number grid points but the maximum allowed grid index
 * For example, ng is the number of grid point.
 * the actual grid points indices are 0, 1, ..., ng - 1.
 * In a periodic/anti periodic spline, set nmax = ng - 1
 * In a natural boundary spline, set nmax = ng - 2
 * because the end point should be excluded and the last grid point has an index ng - 2.
 */
template<typename T, typename TRESIDUAL>
inline void getSplineBound(T x, TRESIDUAL& dx, int& ind, int nmax)
{
  // lower bound
  if (x < 0)
  {
    ind = 0;
    dx  = T(0);
  }
  else
  {
#if defined(__INTEL_LLVM_COMPILER) || defined(__INTEL_CLANG_COMPILER)
    T ipart = std::floor(x);
    dx = x - ipart;
#else
    T ipart;
    dx  = std::modf(x, &ipart);
#endif
    ind = static_cast<int>(ipart);
    // upper bound
    if (ind > nmax)
    {
      ind = nmax;
      dx  = T(1) - std::numeric_limits<T>::epsilon();
    }
  }
}

/** define computeLocationAndFractional: common to any implementation
 * compute the location of the spline grid point and residual coordinates
 * also it precomputes auxiliary array a, b and c
 */
template<typename T>
inline void computeLocationAndFractional(const typename qmcplusplus::bspline_traits<T, 3>::SplineType* restrict spline_m,
                            T x, T y, T z,
                            int& ix, int& iy, int& iz,
                            T a[4], T b[4], T c[4])
{
  x -= spline_m->x_grid.start;
  y -= spline_m->y_grid.start;
  z -= spline_m->z_grid.start;

  T tx, ty, tz;

  getSplineBound(x * spline_m->x_grid.delta_inv, tx, ix, spline_m->x_grid.num - 1);
  getSplineBound(y * spline_m->y_grid.delta_inv, ty, iy, spline_m->y_grid.num - 1);
  getSplineBound(z * spline_m->z_grid.delta_inv, tz, iz, spline_m->z_grid.num - 1);

  MultiBsplineData<T>::compute_prefactors(a, tx);
  MultiBsplineData<T>::compute_prefactors(b, ty);
  MultiBsplineData<T>::compute_prefactors(c, tz);
}

/** define computeLocationAndFractional: common to any implementation
 * compute the location of the spline grid point and residual coordinates
 * also it precomputes auxiliary array (a,b,c) (da,db,dc) (d2a,d2b,d2c)
 */
template<typename T>
inline void computeLocationAndFractional(const typename qmcplusplus::bspline_traits<T, 3>::SplineType* restrict spline_m,
                            T x, T y, T z,
                            int& ix, int& iy, int& iz,
                            T a[4], T b[4], T c[4],
                            T da[4], T db[4], T dc[4],
                            T d2a[4], T d2b[4], T d2c[4])
{
  x -= spline_m->x_grid.start;
  y -= spline_m->y_grid.start;
  z -= spline_m->z_grid.start;

  T tx, ty, tz;

  getSplineBound(x * spline_m->x_grid.delta_inv, tx, ix, spline_m->x_grid.num - 1);
  getSplineBound(y * spline_m->y_grid.delta_inv, ty, iy, spline_m->y_grid.num - 1);
  getSplineBound(z * spline_m->z_grid.delta_inv, tz, iz, spline_m->z_grid.num - 1);

  MultiBsplineData<T>::compute_prefactors(a, da, d2a, tx);
  MultiBsplineData<T>::compute_prefactors(b, db, d2b, ty);
  MultiBsplineData<T>::compute_prefactors(c, dc, d2c, tz);
}

} // namespace spline2

#endif
