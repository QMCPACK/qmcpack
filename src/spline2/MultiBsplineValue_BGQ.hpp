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
#ifndef SPLINE2_MULTIEINSPLINE_VALUE_BGQ_HPP
#define SPLINE2_MULTIEINSPLINE_VALUE_BGQ_HPP

#if defined(__xlC__)
#include <builtins.h>
#endif

namespace spline2
{
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

  vector4double vec_c0 = vec_splats(c[0]);
  vector4double vec_c1 = vec_splats(c[1]);
  vector4double vec_c2 = vec_splats(c[2]);
  vector4double vec_c3 = vec_splats(c[3]);

  const intptr_t xs = spline_m->x_stride;
  const intptr_t ys = spline_m->y_stride;
  const intptr_t zs = spline_m->z_stride;

  constexpr T zero(0);
  const int num_splines = last - first;
  std::fill(vals, vals + num_splines, zero);

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
    {
      const T pre00           = a[i] * b[j];
      vector4double vec_pre00 = vec_splats(pre00);
      T* restrict coefs0      = spline_m->coefs + ((ix + i) * xs + (iy + j) * ys + iz * zs) + first;
      T* restrict coefs1      = coefs0 + zs;
      T* restrict coefs2      = coefs0 + 2 * zs;
      T* restrict coefs3      = coefs0 + 3 * zs;
      for (int n = 0, p = 0; n < num_splines; n += 4, p += 4 * sizeof(T))
      {
        vector4double vec_coef0, vec_coef1, vec_coef2, vec_coef3, vec_val;

        __dcbt(&coefs0[n + 8]);
        __dcbt(&coefs1[n + 8]);
        __dcbt(&coefs2[n + 8]);
        __dcbt(&coefs3[n + 8]);

        vec_coef0 = vec_ld(p, coefs0);
        vec_coef1 = vec_ld(p, coefs1);
        vec_coef2 = vec_ld(p, coefs2);
        vec_coef3 = vec_ld(p, coefs3);
        vec_val   = vec_ld(p, vals);

        vec_coef0 = vec_mul(vec_c0, vec_coef0);
        vec_coef0 = vec_madd(vec_c1, vec_coef1, vec_coef0);
        vec_coef0 = vec_madd(vec_c2, vec_coef2, vec_coef0);
        vec_coef0 = vec_madd(vec_c3, vec_coef3, vec_coef0);
        vec_val   = vec_madd(vec_pre00, vec_coef0, vec_val);
        vec_st(vec_val, p, vals);
      }
    }
}

} // namespace spline2
#endif
