//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file MultiBsplineBase.hpp
 *
 * The baseline implementation with grad and hessian in AoS structure.
 * LHS has 3-stride data access patterns
 */
#ifndef QMCPLUSPLUS_MULTIEINSPLINE_STD3_ENGINE_HPP
#define QMCPLUSPLUS_MULTIEINSPLINE_STD3_ENGINE_HPP

namespace qmcplusplus
{

  template<typename T>
    inline void 
    MultiBspline<T>::evaluate_vgl_impl(T x, T y, T z, 
        T* restrict vals, T* restrict grads, T* restrict lapl, int first, int last,size_t out_offset) const
    {
      x -= spline_m->x_grid.start;
      y -= spline_m->y_grid.start;
      z -= spline_m->z_grid.start;
      T tx,ty,tz;
      int ix,iy,iz;
      SplineBound<T>::get(x*spline_m->x_grid.delta_inv,tx,ix,spline_m->x_grid.num-1);
      SplineBound<T>::get(y*spline_m->y_grid.delta_inv,ty,iy,spline_m->y_grid.num-1);
      SplineBound<T>::get(z*spline_m->z_grid.delta_inv,tz,iz,spline_m->z_grid.num-1);

      T a[4],b[4],c[4],da[4],db[4],dc[4],d2a[4],d2b[4],d2c[4];

      MultiBsplineData<T>::compute_prefactors(a, da, d2a, tx);
      MultiBsplineData<T>::compute_prefactors(b, db, d2b, ty);
      MultiBsplineData<T>::compute_prefactors(c, dc, d2c, tz);

      const intptr_t xs = spline_m->x_stride;
      const intptr_t ys = spline_m->y_stride;
      const intptr_t zs = spline_m->z_stride;

      out_offset=(out_offset)?out_offset:spline_m->num_splines;
      const int num_splines=last-first;

      std::fill(vals,vals+num_splines,T());
      std::fill(grads,grads+3*num_splines,T());

      T lap3[num_splines*3];
      std::fill(lap3,lap3+3*num_splines,T());

      for (int i=0; i<4; i++)
        for (int j=0; j<4; j++){

          const T* restrict coefs = spline_m->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs)+first;

          const T pre20 = d2a[i]*  b[j];
          const T pre10 =  da[i]*  b[j];
          const T pre00 =   a[i]*  b[j];
          const T pre11 =  da[i]* db[j];
          const T pre01 =   a[i]* db[j];
          const T pre02 =   a[i]*d2b[j];
#if defined(ENABLE_EINSPLINE_SIMD)
#pragma omp simd
#endif
          for (int n=0; n<num_splines; n++) {
            T sum0 =   c[0] * coefs[n] +   c[1] * coefs[n+zs] +   c[2] * coefs[n+2*zs] +   c[3] * coefs[n+3*zs];
            T sum1 =  dc[0] * coefs[n] +  dc[1] * coefs[n+zs] +  dc[2] * coefs[n+2*zs] +  dc[3] * coefs[n+3*zs];
            T sum2 = d2c[0] * coefs[n] + d2c[1] * coefs[n+zs] + d2c[2] * coefs[n+2*zs] + d2c[3] * coefs[n+3*zs];
            grads[3*n+0] += pre10 * sum0;
            grads[3*n+1] += pre01 * sum0;
            grads[3*n+2] += pre00 * sum1;
            lap3[3*n+0] += pre20 * sum0;
            lap3[3*n+1] += pre02 * sum0;
            lap3[3*n+2] += pre00 * sum2;
            vals[n] += pre00 * sum0;
          }
        }

      const T dxInv = spline_m->x_grid.delta_inv;
      const T dyInv = spline_m->y_grid.delta_inv;
      const T dzInv = spline_m->z_grid.delta_inv;

      const T dxInv2 = dxInv*dxInv;
      const T dyInv2 = dyInv*dyInv;
      const T dzInv2 = dzInv*dzInv;

#if defined(ENABLE_EINSPLINE_SIMD)
#pragma omp simd
#endif
      for (int n=0; n<num_splines; n++) 
      {
        grads[3*n+0] *= dxInv;
        grads[3*n+1] *= dyInv;
        grads[3*n+2] *= dzInv;
        lapl[n] = lap3[3*n+0] + lap3[3*n+1] + lap3[3*n+2];
      }
    }

  template<typename T>
    inline void 
    MultiBspline<T>::evaluate_vgh_impl(T x, T y, T z, 
        T* restrict vals, T* restrict grads, T* restrict hess, int first, int last, size_t out_offset) const
    {
      x -= spline_m->x_grid.start;
      y -= spline_m->y_grid.start;
      z -= spline_m->z_grid.start;
      T tx,ty,tz;
      int ix,iy,iz;
      SplineBound<T>::get(x*spline_m->x_grid.delta_inv,tx,ix,spline_m->x_grid.num-1);
      SplineBound<T>::get(y*spline_m->y_grid.delta_inv,ty,iy,spline_m->y_grid.num-1);
      SplineBound<T>::get(z*spline_m->z_grid.delta_inv,tz,iz,spline_m->z_grid.num-1);
      T a[4],b[4],c[4],da[4],db[4],dc[4],d2a[4],d2b[4],d2c[4];

      MultiBsplineData<T>::compute_prefactors(a, da, d2a, tx);
      MultiBsplineData<T>::compute_prefactors(b, db, d2b, ty);
      MultiBsplineData<T>::compute_prefactors(c, dc, d2c, tz);

      const intptr_t xs = spline_m->x_stride;
      const intptr_t ys = spline_m->y_stride;
      const intptr_t zs = spline_m->z_stride;

      out_offset=(out_offset)?out_offset:spline_m->num_splines;
      const int num_splines=last-first;

      std::fill(vals,vals+num_splines,T());
      std::fill(grads,grads+3*num_splines,T());
      std::fill(hess,hess+6*num_splines,T());

      for (int i=0; i<4; i++)
        for (int j=0; j<4; j++){
          const T* restrict coefs = spline_m->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs)+first;
          const T pre20 = d2a[i]*  b[j];
          const T pre10 =  da[i]*  b[j];
          const T pre00 =   a[i]*  b[j];
          const T pre11 =  da[i]* db[j];
          const T pre01 =   a[i]* db[j];
          const T pre02 =   a[i]*d2b[j];

#if defined(ENABLE_EINSPLINE_SIMD)
#pragma omp simd
#endif
          for (int n=0; n<num_splines; n++) {
            T sum0 =   c[0] * coefs[n] +   c[1] * coefs[n+zs] +   c[2] * coefs[n+2*zs] +   c[3] * coefs[n+3*zs];
            T sum1 =  dc[0] * coefs[n] +  dc[1] * coefs[n+zs] +  dc[2] * coefs[n+2*zs] +  dc[3] * coefs[n+3*zs];
            T sum2 = d2c[0] * coefs[n] + d2c[1] * coefs[n+zs] + d2c[2] * coefs[n+2*zs] + d2c[3] * coefs[n+3*zs];
            hess[6*n+0] += pre20 * sum0;
            hess[6*n+1] += pre11 * sum0;
            hess[6*n+2] += pre10 * sum1;
            hess[6*n+3] += pre02 * sum0;
            hess[6*n+4] += pre01 * sum1;
            hess[6*n+5] += pre00 * sum2;
            grads[3*n+0] += pre10 * sum0;
            grads[3*n+1] += pre01 * sum0;
            grads[3*n+0] += pre00 * sum1;
            vals[n]+= pre00 * sum0;
          }
        }

      const T dxInv = spline_m->x_grid.delta_inv;
      const T dyInv = spline_m->y_grid.delta_inv;
      const T dzInv = spline_m->z_grid.delta_inv;
      const T dxx=dxInv*dxInv;
      const T dyy=dyInv*dyInv;
      const T dzz=dzInv*dzInv;
      const T dxy=dxInv*dyInv;
      const T dxz=dxInv*dzInv;
      const T dyz=dyInv*dzInv;


#if defined(ENABLE_EINSPLINE_SIMD)
#pragma omp simd
#endif
      for (int n=0; n<num_splines; n++)
      {
        grads[3*n+0]*=dxInv; 
        grads[3*n+1]*=dyInv; 
        grads[3*n+0]*=dzInv; 
        hess[6*n+0]*=dxx;
        hess[6*n+1]*=dxy;
        hess[6*n+2]*=dxz;
        hess[6*n+3]*=dyy;
        hess[6*n+4]*=dyz;
        hess[6*n+5]*=dzz;
      }

    }
}/** qmcplusplus namespace */
#endif

