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
#include "Numerics/SplineBound.hpp"

namespace spline2
{

/** define computeLocationAndFractional: common to any implementation
 * compute the location of the spline grid point and residual coordinates
 * also it precomputes auxiliary array a, b and c
 */
template<typename T>
inline void computeLocationAndFractional(
    const typename qmcplusplus::bspline_traits<T, 3>::SplineType* restrict spline_m,
    T x,
    T y,
    T z,
    int& ix,
    int& iy,
    int& iz,
    T a[4],
    T b[4],
    T c[4])
{
  x -= spline_m->x_grid.start;
  y -= spline_m->y_grid.start;
  z -= spline_m->z_grid.start;

  T tx, ty, tz;

  qmcplusplus::getSplineBound(x * spline_m->x_grid.delta_inv, spline_m->x_grid.num - 1, ix, tx);
  qmcplusplus::getSplineBound(y * spline_m->y_grid.delta_inv, spline_m->y_grid.num - 1, iy, ty);
  qmcplusplus::getSplineBound(z * spline_m->z_grid.delta_inv, spline_m->z_grid.num - 1, iz, tz);

  MultiBsplineData<T>::compute_prefactors(a, tx);
  MultiBsplineData<T>::compute_prefactors(b, ty);
  MultiBsplineData<T>::compute_prefactors(c, tz);
}

/** define computeLocationAndFractional: common to any implementation
 * compute the location of the spline grid point and residual coordinates
 * also it precomputes auxiliary array (a,b,c) (da,db,dc) (d2a,d2b,d2c)
 */
template<typename T>
inline void computeLocationAndFractional(
    const typename qmcplusplus::bspline_traits<T, 3>::SplineType* restrict spline_m,
    T x,
    T y,
    T z,
    int& ix,
    int& iy,
    int& iz,
    T a[4],
    T b[4],
    T c[4],
    T da[4],
    T db[4],
    T dc[4],
    T d2a[4],
    T d2b[4],
    T d2c[4])
{
  x -= spline_m->x_grid.start;
  y -= spline_m->y_grid.start;
  z -= spline_m->z_grid.start;

  T tx, ty, tz;

  qmcplusplus::getSplineBound(x * spline_m->x_grid.delta_inv, spline_m->x_grid.num - 1, ix, tx);
  qmcplusplus::getSplineBound(y * spline_m->y_grid.delta_inv, spline_m->y_grid.num - 1, iy, ty);
  qmcplusplus::getSplineBound(z * spline_m->z_grid.delta_inv, spline_m->z_grid.num - 1, iz, tz);

  MultiBsplineData<T>::compute_prefactors(a, da, d2a, tx);
  MultiBsplineData<T>::compute_prefactors(b, db, d2b, ty);
  MultiBsplineData<T>::compute_prefactors(c, dc, d2c, tz);
}

} // namespace spline2

#endif
