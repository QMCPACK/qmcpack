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
#ifndef SPLINE2_MULTIEINSPLINE_VALUE_STD3_HPP
#define SPLINE2_MULTIEINSPLINE_VALUE_STD3_HPP

namespace spline2
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
  int ix, iy, iz;
  T a[4], b[4], c[4];

  computeLocationAndFractional(spline_m, x, y, z, ix, iy, iz, a, b, c);

  const intptr_t xs = spline_m->x_stride;
  const intptr_t ys = spline_m->y_stride;
  const intptr_t zs = spline_m->z_stride;

  constexpr T zero(0);
  const int num_splines = last - first;
  std::fill(vals, vals + num_splines, zero);

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
    {
      const T pre00              = a[i] * b[j];
      const T* restrict coefs    = spline_m->coefs + ((ix + i) * xs + (iy + j) * ys + iz * zs) + first;
      const T* restrict coefszs  = coefs + zs;
      const T* restrict coefs2zs = coefs + 2 * zs;
      const T* restrict coefs3zs = coefs + 3 * zs;
#pragma omp simd aligned(coefs, coefszs, coefs2zs, coefs3zs, vals: QMC_SIMD_ALIGNMENT)
      for (int n = 0; n < num_splines; n++)
        vals[n] += pre00 * (c[0] * coefs[n] + c[1] * coefszs[n] + c[2] * coefs2zs[n] + c[3] * coefs3zs[n]);
    }
}

} // namespace spline2
#endif
