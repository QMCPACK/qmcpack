//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef SPLINE2_OFFLOAD_MULTIEINSPLINE_VALUE_HPP
#define SPLINE2_OFFLOAD_MULTIEINSPLINE_VALUE_HPP

#include "OMPTarget/OMPstd.hpp"
#include "spline2/MultiBsplineData.hpp"
#include "spline2/MultiBsplineEval_helper.hpp"

namespace spline2offload
{
/** define evaluate: common to any implementation */
template<typename T>
inline void evaluate_v_impl(const typename qmcplusplus::bspline_traits<T, 3>::SplineType* restrict spline_m,
                            T x,
                            T y,
                            T z,
                            T* restrict vals,
                            int first,
                            int last)
{
  x -= spline_m->x_grid.start;
  y -= spline_m->y_grid.start;
  z -= spline_m->z_grid.start;
  T tx, ty, tz;
  int ix, iy, iz;
  spline2::getSplineBound(x * spline_m->x_grid.delta_inv, tx, ix, spline_m->x_grid.num - 1);
  spline2::getSplineBound(y * spline_m->y_grid.delta_inv, ty, iy, spline_m->y_grid.num - 1);
  spline2::getSplineBound(z * spline_m->z_grid.delta_inv, tz, iz, spline_m->z_grid.num - 1);
  T a[4], b[4], c[4];

  spline2::MultiBsplineData<T>::compute_prefactors(a, tx);
  spline2::MultiBsplineData<T>::compute_prefactors(b, ty);
  spline2::MultiBsplineData<T>::compute_prefactors(c, tz);

  const intptr_t xs = spline_m->x_stride;
  const intptr_t ys = spline_m->y_stride;
  const intptr_t zs = spline_m->z_stride;

  const int num_splines = last - first;
  OMPstd::fill_n(vals, num_splines, T());

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
    {
      const T pre00              = a[i] * b[j];
      const T* restrict coefs    = spline_m->coefs + ((ix + i) * xs + (iy + j) * ys + iz * zs) + first;
      const T* restrict coefszs  = coefs + zs;
      const T* restrict coefs2zs = coefs + 2 * zs;
      const T* restrict coefs3zs = coefs + 3 * zs;
#ifdef ENABLE_OFFLOAD
#pragma omp for
#else
#pragma omp simd aligned(coefs, coefszs, coefs2zs, coefs3zs, vals : QMC_SIMD_ALIGNMENT)
#endif
      for (int n = 0; n < num_splines; n++)
        vals[n] += pre00 * (c[0] * coefs[n] + c[1] * coefszs[n] + c[2] * coefs2zs[n] + c[3] * coefs3zs[n]);
    }
}

/** evaluate value, gradients and hessian
 * @param ix location in x
 * @param iy location in y
 * @param iz location in z
 * @param index spline index
 * @param a interpolation parameter in x
 * @param b interpolation parameter in y
 * @param c interpolation parameter in z
 * @param vals value output location
 */
template<typename T>
inline void evaluate_v_impl_v2(const typename qmcplusplus::bspline_traits<T, 3>::SplineType* restrict spline_m,
                               int ix,
                               int iy,
                               int iz,
                               const int index,
                               const T a[4],
                               const T b[4],
                               const T c[4],
                               T* restrict vals)
{
  const intptr_t xs = spline_m->x_stride;
  const intptr_t ys = spline_m->y_stride;
  const intptr_t zs = spline_m->z_stride;

  T val = T();
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
    {
      const T* restrict coefs = spline_m->coefs + ((ix + i) * xs + (iy + j) * ys + iz * zs);
      val += a[i] * b[j] *
          (c[0] * coefs[index] + c[1] * coefs[index + zs] + c[2] * coefs[index + zs * 2] +
           c[3] * coefs[index + zs * 3]);
    }
  *vals = val;
}

} // namespace spline2offload
#endif
