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
/**@file MultiBsplineStd.hpp
 *
 * Literal port of einspline/multi_bspline_eval_d_std3.cpp by Ye at ANL
 * Modified to handle both float/double and blocked operations
 * - template parameter T for the precision
 * - MUB spline object created by einspline allocators.
 * Function signatures modified anticipating its use by a class that can perform data parallel execution
 * - evaluate(...., int first, int last)
 */
#ifndef SPLINE2_MULTIEINSPLINE_VGLH_STD3_HPP
#define SPLINE2_MULTIEINSPLINE_VGLH_STD3_HPP

namespace spline2
{
template<typename T>
inline void evaluate_vgl_impl(const typename qmcplusplus::bspline_traits<T, 3>::SplineType* restrict spline_m,
                              T x,
                              T y,
                              T z,
                              T* restrict vals,
                              T* restrict grads,
                              T* restrict lapl,
                              size_t out_offset,
                              int first,
                              int last)
{
  int ix, iy, iz;
  T a[4], b[4], c[4], da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];

  computeLocationAndFractional(spline_m, x, y, z, ix, iy, iz, a, b, c, da, db, dc, d2a, d2b, d2c);

  const intptr_t xs = spline_m->x_stride;
  const intptr_t ys = spline_m->y_stride;
  const intptr_t zs = spline_m->z_stride;

  const int num_splines = last - first;

  T* restrict gx = grads;
  T* restrict gy = grads + out_offset;
  T* restrict gz = grads + 2 * out_offset;
  T* restrict lx = lapl;
  T* restrict ly = lapl + out_offset;
  T* restrict lz = lapl + 2 * out_offset;

  std::fill(vals, vals + num_splines, T());
  std::fill(gx, gx + num_splines, T());
  std::fill(gy, gy + num_splines, T());
  std::fill(gz, gz + num_splines, T());
  std::fill(lx, lx + num_splines, T());
  std::fill(ly, ly + num_splines, T());
  std::fill(lz, lz + num_splines, T());

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
    {
      const T pre20 = d2a[i] * b[j];
      const T pre10 = da[i] * b[j];
      const T pre00 = a[i] * b[j];
      const T pre01 = a[i] * db[j];
      const T pre02 = a[i] * d2b[j];

      const T* restrict coefs    = spline_m->coefs + ((ix + i) * xs + (iy + j) * ys + iz * zs) + first;
      const T* restrict coefszs  = coefs + zs;
      const T* restrict coefs2zs = coefs + 2 * zs;
      const T* restrict coefs3zs = coefs + 3 * zs;

#pragma omp simd aligned(coefs, coefszs, coefs2zs, coefs3zs, gx, gy, gz, lx, ly, lz, vals: QMC_SIMD_ALIGNMENT)
      for (int n = 0; n < num_splines; n++)
      {
        const T coefsv    = coefs[n];
        const T coefsvzs  = coefszs[n];
        const T coefsv2zs = coefs2zs[n];
        const T coefsv3zs = coefs3zs[n];

        T sum0 = c[0] * coefsv + c[1] * coefsvzs + c[2] * coefsv2zs + c[3] * coefsv3zs;
        T sum1 = dc[0] * coefsv + dc[1] * coefsvzs + dc[2] * coefsv2zs + dc[3] * coefsv3zs;
        T sum2 = d2c[0] * coefsv + d2c[1] * coefsvzs + d2c[2] * coefsv2zs + d2c[3] * coefsv3zs;
        gx[n] += pre10 * sum0;
        gy[n] += pre01 * sum0;
        gz[n] += pre00 * sum1;
        lx[n] += pre20 * sum0;
        ly[n] += pre02 * sum0;
        lz[n] += pre00 * sum2;
        vals[n] += pre00 * sum0;
      }
    }

  const T dxInv = spline_m->x_grid.delta_inv;
  const T dyInv = spline_m->y_grid.delta_inv;
  const T dzInv = spline_m->z_grid.delta_inv;

  const T dxInv2 = dxInv * dxInv;
  const T dyInv2 = dyInv * dyInv;
  const T dzInv2 = dzInv * dzInv;

#pragma omp simd aligned(gx, gy, gz, lx: QMC_SIMD_ALIGNMENT)
  for (int n = 0; n < num_splines; n++)
  {
    gx[n] *= dxInv;
    gy[n] *= dyInv;
    gz[n] *= dzInv;
    lx[n] = lx[n] * dxInv2 + ly[n] * dyInv2 + lz[n] * dzInv2;
  }
}

template<typename T>
inline void evaluate_vgh_impl(const typename qmcplusplus::bspline_traits<T, 3>::SplineType* restrict spline_m,
                              T x,
                              T y,
                              T z,
                              T* restrict vals,
                              T* restrict grads,
                              T* restrict hess,
                              size_t out_offset,
                              int first,
                              int last)
{
  int ix, iy, iz;
  T a[4], b[4], c[4], da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];

  computeLocationAndFractional(spline_m, x, y, z, ix, iy, iz, a, b, c, da, db, dc, d2a, d2b, d2c);

  const intptr_t xs = spline_m->x_stride;
  const intptr_t ys = spline_m->y_stride;
  const intptr_t zs = spline_m->z_stride;

  const int num_splines = last - first;

  T* restrict gx = grads;
  T* restrict gy = grads + out_offset;
  T* restrict gz = grads + 2 * out_offset;

  T* restrict hxx = hess;
  T* restrict hxy = hess + out_offset;
  T* restrict hxz = hess + 2 * out_offset;
  T* restrict hyy = hess + 3 * out_offset;
  T* restrict hyz = hess + 4 * out_offset;
  T* restrict hzz = hess + 5 * out_offset;

  std::fill(vals, vals + num_splines, T());
  std::fill(gx, gx + num_splines, T());
  std::fill(gy, gy + num_splines, T());
  std::fill(gz, gz + num_splines, T());
  std::fill(hxx, hxx + num_splines, T());
  std::fill(hxy, hxy + num_splines, T());
  std::fill(hxz, hxz + num_splines, T());
  std::fill(hyy, hyy + num_splines, T());
  std::fill(hyz, hyz + num_splines, T());
  std::fill(hzz, hzz + num_splines, T());

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
    {
      const T* restrict coefs    = spline_m->coefs + ((ix + i) * xs + (iy + j) * ys + iz * zs) + first;
      const T* restrict coefszs  = coefs + zs;
      const T* restrict coefs2zs = coefs + 2 * zs;
      const T* restrict coefs3zs = coefs + 3 * zs;

      const T pre20 = d2a[i] * b[j];
      const T pre10 = da[i] * b[j];
      const T pre00 = a[i] * b[j];
      const T pre11 = da[i] * db[j];
      const T pre01 = a[i] * db[j];
      const T pre02 = a[i] * d2b[j];

#pragma omp simd aligned(coefs, coefszs, coefs2zs, coefs3zs, gx, gy, gz, hxx, hxy, hxz, hyy, hyz, hzz, vals: QMC_SIMD_ALIGNMENT)
      for (int n = 0; n < num_splines; n++)
      {
        T coefsv    = coefs[n];
        T coefsvzs  = coefszs[n];
        T coefsv2zs = coefs2zs[n];
        T coefsv3zs = coefs3zs[n];

        T sum0 = c[0] * coefsv + c[1] * coefsvzs + c[2] * coefsv2zs + c[3] * coefsv3zs;
        T sum1 = dc[0] * coefsv + dc[1] * coefsvzs + dc[2] * coefsv2zs + dc[3] * coefsv3zs;
        T sum2 = d2c[0] * coefsv + d2c[1] * coefsvzs + d2c[2] * coefsv2zs + d2c[3] * coefsv3zs;

        hxx[n] += pre20 * sum0;
        hxy[n] += pre11 * sum0;
        hxz[n] += pre10 * sum1;
        hyy[n] += pre02 * sum0;
        hyz[n] += pre01 * sum1;
        hzz[n] += pre00 * sum2;
        gx[n] += pre10 * sum0;
        gy[n] += pre01 * sum0;
        gz[n] += pre00 * sum1;
        vals[n] += pre00 * sum0;
      }
    }

  const T dxInv = spline_m->x_grid.delta_inv;
  const T dyInv = spline_m->y_grid.delta_inv;
  const T dzInv = spline_m->z_grid.delta_inv;
  const T dxx   = dxInv * dxInv;
  const T dyy   = dyInv * dyInv;
  const T dzz   = dzInv * dzInv;
  const T dxy   = dxInv * dyInv;
  const T dxz   = dxInv * dzInv;
  const T dyz   = dyInv * dzInv;

#pragma omp simd aligned(gx, gy, gz, hxx, hxy, hxz, hyy, hyz, hzz: QMC_SIMD_ALIGNMENT)
  for (int n = 0; n < num_splines; n++)
  {
    gx[n] *= dxInv;
    gy[n] *= dyInv;
    gz[n] *= dzInv;
    hxx[n] *= dxx;
    hyy[n] *= dyy;
    hzz[n] *= dzz;
    hxy[n] *= dxy;
    hxz[n] *= dxz;
    hyz[n] *= dyz;
  }
}

} // namespace spline2
#endif
