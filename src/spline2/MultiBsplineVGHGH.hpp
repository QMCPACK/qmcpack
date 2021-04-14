//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
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
#ifndef SPLINE2_MULTIEINSPLINE_VGHGH_HPP
#define SPLINE2_MULTIEINSPLINE_VGHGH_HPP

namespace spline2
{
template<typename T>
inline void evaluate_vghgh_impl(const typename qmcplusplus::bspline_traits<T, 3>::SplineType* restrict spline_m,
                                T x,
                                T y,
                                T z,
                                T* restrict vals,
                                T* restrict grads,
                                T* restrict hess,
                                T* restrict ghess,
                                size_t out_offset,
                                int first,
                                int last)
{
  int ix, iy, iz;
  T tx, ty, tz;
  T a[4], b[4], c[4];
  T da[4], db[4], dc[4];
  T d2a[4], d2b[4], d2c[4];
  T d3a[4], d3b[4], d3c[4];

  x -= spline_m->x_grid.start;
  y -= spline_m->y_grid.start;
  z -= spline_m->z_grid.start;
  getSplineBound(x * spline_m->x_grid.delta_inv, tx, ix, spline_m->x_grid.num - 1);
  getSplineBound(y * spline_m->y_grid.delta_inv, ty, iy, spline_m->y_grid.num - 1);
  getSplineBound(z * spline_m->z_grid.delta_inv, tz, iz, spline_m->z_grid.num - 1);

  MultiBsplineData<T>::compute_prefactors(a, da, d2a, d3a, tx);
  MultiBsplineData<T>::compute_prefactors(b, db, d2b, d3b, ty);
  MultiBsplineData<T>::compute_prefactors(c, dc, d2c, d3c, tz);

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

  T* restrict gh_xxx = ghess;
  T* restrict gh_xxy = ghess + out_offset;
  T* restrict gh_xxz = ghess + 2 * out_offset;
  T* restrict gh_xyy = ghess + 3 * out_offset;
  T* restrict gh_xyz = ghess + 4 * out_offset;
  T* restrict gh_xzz = ghess + 5 * out_offset;
  T* restrict gh_yyy = ghess + 6 * out_offset;
  T* restrict gh_yyz = ghess + 7 * out_offset;
  T* restrict gh_yzz = ghess + 8 * out_offset;
  T* restrict gh_zzz = ghess + 9 * out_offset;

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

  std::fill(gh_xxx, gh_xxx + num_splines, T());
  std::fill(gh_xxy, gh_xxy + num_splines, T());
  std::fill(gh_xxz, gh_xxz + num_splines, T());
  std::fill(gh_xyy, gh_xyy + num_splines, T());
  std::fill(gh_xyz, gh_xyz + num_splines, T());
  std::fill(gh_xzz, gh_xzz + num_splines, T());
  std::fill(gh_yyy, gh_yyy + num_splines, T());
  std::fill(gh_yyz, gh_yyz + num_splines, T());
  std::fill(gh_yzz, gh_yzz + num_splines, T());
  std::fill(gh_zzz, gh_zzz + num_splines, T());

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

      const T pre30 = d3a[i] * b[j];
      const T pre21 = d2a[i] * db[j];
      const T pre12 = da[i] * d2b[j];
      const T pre03 = a[i] * d3b[j];


#pragma omp simd aligned(coefs, coefszs, coefs2zs, coefs3zs, gx, gy, gz, hxx, hxy, hxz, hyy, hyz, hzz, gh_xxx, gh_xxy, \
                         gh_xxz, gh_xyy, gh_xyz, gh_xzz, gh_yyy, gh_yyz, gh_yzz, gh_zzz, vals: QMC_SIMD_ALIGNMENT)
      for (int n = 0; n < num_splines; n++)
      {
        T coefsv    = coefs[n];
        T coefsvzs  = coefszs[n];
        T coefsv2zs = coefs2zs[n];
        T coefsv3zs = coefs3zs[n];

        T sum0 = c[0] * coefsv + c[1] * coefsvzs + c[2] * coefsv2zs + c[3] * coefsv3zs;
        T sum1 = dc[0] * coefsv + dc[1] * coefsvzs + dc[2] * coefsv2zs + dc[3] * coefsv3zs;
        T sum2 = d2c[0] * coefsv + d2c[1] * coefsvzs + d2c[2] * coefsv2zs + d2c[3] * coefsv3zs;
        T sum3 = d3c[0] * coefsv + d3c[1] * coefsvzs + d3c[2] * coefsv2zs + d3c[3] * coefsv3zs;

        gh_xxx[n] += pre30 * sum0;
        gh_xxy[n] += pre21 * sum0;
        gh_xxz[n] += pre20 * sum1;
        gh_xyy[n] += pre12 * sum0;
        gh_xyz[n] += pre11 * sum1;
        gh_xzz[n] += pre10 * sum2;
        gh_yyy[n] += pre03 * sum0;
        gh_yyz[n] += pre02 * sum1;
        gh_yzz[n] += pre01 * sum2;
        gh_zzz[n] += pre00 * sum3;

        hxx[n]  += pre20 * sum0;
        hxy[n]  += pre11 * sum0;
        hxz[n]  += pre10 * sum1;
        hyy[n]  += pre02 * sum0;
        hyz[n]  += pre01 * sum1;
        hzz[n]  += pre00 * sum2;
        gx[n]   += pre10 * sum0;
        gy[n]   += pre01 * sum0;
        gz[n]   += pre00 * sum1;
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

  const T dxxx = dxInv * dxInv * dxInv;
  const T dxxy = dxInv * dxInv * dyInv;
  const T dxxz = dxInv * dxInv * dzInv;
  const T dxyy = dxInv * dyInv * dyInv;
  const T dxyz = dxInv * dyInv * dzInv;
  const T dxzz = dxInv * dzInv * dzInv;
  const T dyyy = dyInv * dyInv * dyInv;
  const T dyyz = dyInv * dyInv * dzInv;
  const T dyzz = dyInv * dzInv * dzInv;
  const T dzzz = dzInv * dzInv * dzInv;

#pragma omp simd aligned(gx, gy, gz, hxx, hxy, hxz, hyy, hyz, hzz, gh_xxx, gh_xxy, gh_xxz, gh_xyy, gh_xyz, gh_xzz, \
                         gh_yyy, gh_yyz, gh_yzz, gh_zzz: QMC_SIMD_ALIGNMENT)
  for (int n = 0; n < num_splines; n++)
  {
    gx[n]  *= dxInv;
    gy[n]  *= dyInv;
    gz[n]  *= dzInv;
    hxx[n] *= dxx;
    hyy[n] *= dyy;
    hzz[n] *= dzz;
    hxy[n] *= dxy;
    hxz[n] *= dxz;
    hyz[n] *= dyz;

    gh_xxx[n] *= dxxx;
    gh_xxy[n] *= dxxy;
    gh_xxz[n] *= dxxz;
    gh_xyy[n] *= dxyy;
    gh_xyz[n] *= dxyz;
    gh_xzz[n] *= dxzz;
    gh_yyy[n] *= dyyy;
    gh_yyz[n] *= dyyz;
    gh_yzz[n] *= dyzz;
    gh_zzz[n] *= dzzz;
  }
}

} // namespace spline2
#endif
