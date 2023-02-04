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
#ifndef SPLINE2_OFFLOAD_MULTIEINSPLINE_VGLH_HPP
#define SPLINE2_OFFLOAD_MULTIEINSPLINE_VGLH_HPP

#include "OMPTarget/OMPstd.hpp"
#include "spline2/MultiBsplineData.hpp"
#include "spline2/MultiBsplineEval_helper.hpp"

namespace spline2offload
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
  x -= spline_m->x_grid.start;
  y -= spline_m->y_grid.start;
  z -= spline_m->z_grid.start;
  T tx, ty, tz;
  int ix, iy, iz;
  spline2::getSplineBound(x * spline_m->x_grid.delta_inv, tx, ix, spline_m->x_grid.num - 1);
  spline2::getSplineBound(y * spline_m->y_grid.delta_inv, ty, iy, spline_m->y_grid.num - 1);
  spline2::getSplineBound(z * spline_m->z_grid.delta_inv, tz, iz, spline_m->z_grid.num - 1);

  T a[4], b[4], c[4], da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];

  spline2::MultiBsplineData<T>::compute_prefactors(a, da, d2a, tx);
  spline2::MultiBsplineData<T>::compute_prefactors(b, db, d2b, ty);
  spline2::MultiBsplineData<T>::compute_prefactors(c, dc, d2c, tz);

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

  OMPstd::fill_n(vals, num_splines, T());
  OMPstd::fill_n(gx, num_splines, T());
  OMPstd::fill_n(gy, num_splines, T());
  OMPstd::fill_n(gz, num_splines, T());
  OMPstd::fill_n(lx, num_splines, T());
  OMPstd::fill_n(ly, num_splines, T());
  OMPstd::fill_n(lz, num_splines, T());

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
    {
      const T pre20 = d2a[i] * b[j];
      const T pre10 = da[i] * b[j];
      const T pre00 = a[i] * b[j];
      const T pre11 = da[i] * db[j];
      const T pre01 = a[i] * db[j];
      const T pre02 = a[i] * d2b[j];

      const T* restrict coefs    = spline_m->coefs + ((ix + i) * xs + (iy + j) * ys + iz * zs) + first;
      const T* restrict coefszs  = coefs + zs;
      const T* restrict coefs2zs = coefs + 2 * zs;
      const T* restrict coefs3zs = coefs + 3 * zs;

#ifdef ENABLE_OFFLOAD
#pragma omp for
#else
#pragma omp simd aligned(coefs, coefszs, coefs2zs, coefs3zs, gx, gy, gz, lx, ly, lz, vals : QMC_SIMD_ALIGNMENT)
#endif
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

#ifdef ENABLE_OFFLOAD
#pragma omp for
#else
#pragma omp simd aligned(gx, gy, gz, lx : QMC_SIMD_ALIGNMENT)
#endif
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
  T tx, ty, tz;
  T a[4], b[4], c[4], da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];

  x -= spline_m->x_grid.start;
  y -= spline_m->y_grid.start;
  z -= spline_m->z_grid.start;
  spline2::getSplineBound(x * spline_m->x_grid.delta_inv, tx, ix, spline_m->x_grid.num - 1);
  spline2::getSplineBound(y * spline_m->y_grid.delta_inv, ty, iy, spline_m->y_grid.num - 1);
  spline2::getSplineBound(z * spline_m->z_grid.delta_inv, tz, iz, spline_m->z_grid.num - 1);

  spline2::MultiBsplineData<T>::compute_prefactors(a, da, d2a, tx);
  spline2::MultiBsplineData<T>::compute_prefactors(b, db, d2b, ty);
  spline2::MultiBsplineData<T>::compute_prefactors(c, dc, d2c, tz);

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

  OMPstd::fill_n(vals, num_splines, T());
  OMPstd::fill_n(gx, num_splines, T());
  OMPstd::fill_n(gy, num_splines, T());
  OMPstd::fill_n(gz, num_splines, T());
  OMPstd::fill_n(hxx, num_splines, T());
  OMPstd::fill_n(hxy, num_splines, T());
  OMPstd::fill_n(hxz, num_splines, T());
  OMPstd::fill_n(hyy, num_splines, T());
  OMPstd::fill_n(hyz, num_splines, T());
  OMPstd::fill_n(hzz, num_splines, T());

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

#ifdef ENABLE_OFFLOAD
#pragma omp for
#else
#pragma omp simd aligned(coefs, coefszs, coefs2zs, coefs3zs, vals, gx, gy, gz, hxx, hyy, hzz, hxy, hxz, hyz \
                         : QMC_SIMD_ALIGNMENT)
#endif
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

#ifdef ENABLE_OFFLOAD
#pragma omp for
#else
#pragma omp simd aligned(gx, gy, gz, hxx, hyy, hzz, hxy, hxz, hyz : QMC_SIMD_ALIGNMENT)
#endif
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

/** evaluate value, gradients and hessian
 * @param ix location in x
 * @param iy location in y
 * @param iz location in z
 * @param index spline index
 * @param a interpolation parameter in x
 * @param b interpolation parameter in y
 * @param c interpolation parameter in z
 * @param da interpolation parameter, first derivative in x
 * @param db interpolation parameter, first derivative in y
 * @param dc interpolation parameter, first derivative in z
 * @param d2a interpolation parameter, second derivative in x
 * @param d2b interpolation parameter, second derivative in y
 * @param d2c interpolation parameter, second derivative in z
 * @param val_grads_hess Structure-of-Array output in strides of size out_offset
 * @param out_offset The stride of val_grads_hess
 */
template<typename T>
inline void evaluate_vgh_impl_v2(const typename qmcplusplus::bspline_traits<T, 3>::SplineType* restrict spline_m,
                                 int ix,
                                 int iy,
                                 int iz,
                                 const int index,
                                 const T a[4],
                                 const T b[4],
                                 const T c[4],
                                 const T da[4],
                                 const T db[4],
                                 const T dc[4],
                                 const T d2a[4],
                                 const T d2b[4],
                                 const T d2c[4],
                                 T* restrict val_grads_hess,
                                 const size_t out_offset)
{
  const intptr_t xs = spline_m->x_stride;
  const intptr_t ys = spline_m->y_stride;
  const intptr_t zs = spline_m->z_stride;

  T val = T();
  T gx  = T();
  T gy  = T();
  T gz  = T();
  T hxx = T();
  T hxy = T();
  T hxz = T();
  T hyy = T();
  T hyz = T();
  T hzz = T();

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
    {
      const T* restrict coefs = spline_m->coefs + ((ix + i) * xs + (iy + j) * ys + iz * zs);
      const T coefsv          = coefs[index];
      const T coefsvzs        = coefs[index + zs];
      const T coefsv2zs       = coefs[index + 2 * zs];
      const T coefsv3zs       = coefs[index + 3 * zs];


      const T pre20 = d2a[i] * b[j];
      const T pre10 = da[i] * b[j];
      const T pre00 = a[i] * b[j];
      const T pre11 = da[i] * db[j];
      const T pre01 = a[i] * db[j];
      const T pre02 = a[i] * d2b[j];

      T sum0 = c[0] * coefsv + c[1] * coefsvzs + c[2] * coefsv2zs + c[3] * coefsv3zs;
      T sum1 = dc[0] * coefsv + dc[1] * coefsvzs + dc[2] * coefsv2zs + dc[3] * coefsv3zs;
      T sum2 = d2c[0] * coefsv + d2c[1] * coefsvzs + d2c[2] * coefsv2zs + d2c[3] * coefsv3zs;

      hxx += pre20 * sum0;
      hxy += pre11 * sum0;
      hxz += pre10 * sum1;
      hyy += pre02 * sum0;
      hyz += pre01 * sum1;
      hzz += pre00 * sum2;
      gx += pre10 * sum0;
      gy += pre01 * sum0;
      gz += pre00 * sum1;
      val += pre00 * sum0;
    }

  const T dxInv = spline_m->x_grid.delta_inv;
  const T dyInv = spline_m->y_grid.delta_inv;
  const T dzInv = spline_m->z_grid.delta_inv;
  // put data back to the result vector
  val_grads_hess[0]              = val;
  val_grads_hess[out_offset]     = gx * dxInv;
  val_grads_hess[out_offset * 2] = gy * dyInv;
  val_grads_hess[out_offset * 3] = gz * dzInv;
  val_grads_hess[out_offset * 4] = hxx * dxInv * dxInv;
  val_grads_hess[out_offset * 5] = hxy * dxInv * dyInv;
  val_grads_hess[out_offset * 6] = hxz * dxInv * dzInv;
  val_grads_hess[out_offset * 7] = hyy * dyInv * dyInv;
  val_grads_hess[out_offset * 8] = hyz * dyInv * dzInv;
  val_grads_hess[out_offset * 9] = hzz * dzInv * dzInv;
}

} // namespace spline2offload
#endif
