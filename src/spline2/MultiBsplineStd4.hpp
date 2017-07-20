//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                     Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file MultiBsplineStd4.hpp
 *
 * Literal port of einspline/multi_bspline_eval_d_std3.cpp by Ye at ANL
 * Modified to handle both float/double and blocked operations
 * - template parameter T for the precision
 * - MUB spline object created by einspline allocators.
 * Function signatures modified anticipating its use by a class that can perform data parallel execution
 * - evaluate(...., int first, int last)
 */
#ifndef QMCPLUSPLUS_MULTIEINSPLINE_STD4_ENGINE_HPP
#define QMCPLUSPLUS_MULTIEINSPLINE_STD4_ENGINE_HPP

namespace qmcplusplus
{

  template<typename T>
    inline void 
    MultiBspline<T>::evaluate_vgl_impl(T x, T y, T z, 
        T* restrict vals, T* restrict grads, T* restrict lapl, int first, int last, size_t out_offset) const
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

      compute_prefactors(a, da, d2a, tx);
      compute_prefactors(b, db, d2b, ty);
      compute_prefactors(c, dc, d2c, tz);

      const intptr_t xs = spline_m->x_stride;
      const intptr_t ys = spline_m->y_stride;
      const intptr_t zs = spline_m->z_stride;

      out_offset=(out_offset)?out_offset:spline_m->num_splines;
      const int num_splines=last-first;

      ASSUME_ALIGNED(vals);
      T* restrict gx=grads;              ASSUME_ALIGNED(gx);
      T* restrict gy=grads+  out_offset; ASSUME_ALIGNED(gy);
      T* restrict gz=grads+2*out_offset; ASSUME_ALIGNED(gz);
      ASSUME_ALIGNED(lapl);

      const T dxInv = spline_m->x_grid.delta_inv;
      const T dyInv = spline_m->y_grid.delta_inv;
      const T dzInv = spline_m->z_grid.delta_inv;
      const T dxx=dxInv*dxInv;
      const T dyy=dyInv*dyInv;
      const T dzz=dzInv*dzInv;

      #pragma omp simd
      for (int n=0; n<num_splines; n++) {
        T val_local, gx_local, gy_local, gz_local;
        T hxx_local, hyy_local, hzz_local;
        val_local = gx_local = gy_local = gz_local
                  = hxx_local = hyy_local = hzz_local = (T)0.0;
        
        for (int i=0; i<4; i++)
          for (int j=0; j<4; j++){
            const T* restrict coefs = spline_m->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs) + first; ASSUME_ALIGNED(coefs);
            const T* restrict coefszs  = coefs+zs;       ASSUME_ALIGNED(coefszs);
            const T* restrict coefs2zs = coefs+2*zs;     ASSUME_ALIGNED(coefs2zs);
            const T* restrict coefs3zs = coefs+3*zs;     ASSUME_ALIGNED(coefs3zs);

            const T pre20 = d2a[i]*  b[j];
            const T pre10 =  da[i]*  b[j];
            const T pre00 =   a[i]*  b[j];
            const T pre01 =   a[i]* db[j];
            const T pre02 =   a[i]*d2b[j];

            const T sum0 =   c[0] * coefs[n] +   c[1] * coefszs[n] +   c[2] * coefs2zs[n] +   c[3] * coefs3zs[n];
            const T sum1 =  dc[0] * coefs[n] +  dc[1] * coefszs[n] +  dc[2] * coefs2zs[n] +  dc[3] * coefs3zs[n];
            const T sum2 = d2c[0] * coefs[n] + d2c[1] * coefszs[n] + d2c[2] * coefs2zs[n] + d2c[3] * coefs3zs[n];

            hxx_local += pre20 * sum0;
            hyy_local += pre02 * sum0;
            hzz_local += pre00 * sum2;
            gx_local += pre10 * sum0;
            gy_local += pre01 * sum0;
            gz_local += pre00 * sum1;
            val_local += pre00 * sum0;
          }

        lapl[n] =  hxx_local*dxx + hyy_local*dyy + hzz_local*dzz;
        gx[n] = gx_local*dxInv;
        gy[n] = gy_local*dyInv;
        gz[n] = gz_local*dzInv;
        vals[n] = val_local;
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

      compute_prefactors(a, da, d2a, tx);
      compute_prefactors(b, db, d2b, ty);
      compute_prefactors(c, dc, d2c, tz);

      const intptr_t xs = spline_m->x_stride;
      const intptr_t ys = spline_m->y_stride;
      const intptr_t zs = spline_m->z_stride;

      out_offset=(out_offset)?out_offset:spline_m->num_splines;
      const int num_splines=last-first;

      ASSUME_ALIGNED(vals);
      T* restrict gx=grads             ; ASSUME_ALIGNED(gx);
      T* restrict gy=grads  +out_offset; ASSUME_ALIGNED(gy);
      T* restrict gz=grads+2*out_offset; ASSUME_ALIGNED(gz);

      T* restrict hxx=hess             ; ASSUME_ALIGNED(hxx);
      T* restrict hxy=hess+  out_offset; ASSUME_ALIGNED(hxy);
      T* restrict hxz=hess+2*out_offset; ASSUME_ALIGNED(hxz);
      T* restrict hyy=hess+3*out_offset; ASSUME_ALIGNED(hyy);
      T* restrict hyz=hess+4*out_offset; ASSUME_ALIGNED(hyz);
      T* restrict hzz=hess+5*out_offset; ASSUME_ALIGNED(hzz);

      const T dxInv = spline_m->x_grid.delta_inv;
      const T dyInv = spline_m->y_grid.delta_inv;
      const T dzInv = spline_m->z_grid.delta_inv;
      const T dxx=dxInv*dxInv;
      const T dyy=dyInv*dyInv;
      const T dzz=dzInv*dzInv;
      const T dxy=dxInv*dyInv;
      const T dxz=dxInv*dzInv;
      const T dyz=dyInv*dzInv;

      #pragma omp simd
      for (int n=0; n<num_splines; n++) {
        T val_local, gx_local, gy_local, gz_local;
        T hxx_local, hxy_local, hxz_local;
        T hyy_local, hyz_local, hzz_local;
        val_local = gx_local = gy_local = gz_local
                  = hxx_local = hxy_local = hxz_local
                  = hyy_local = hyz_local = hzz_local = (T)0.0;
        
        for (int i=0; i<4; i++)
          for (int j=0; j<4; j++){
            const T* restrict coefs = spline_m->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs) + first; ASSUME_ALIGNED(coefs);
            const T* restrict coefszs  = coefs+zs;       ASSUME_ALIGNED(coefszs);
            const T* restrict coefs2zs = coefs+2*zs;     ASSUME_ALIGNED(coefs2zs);
            const T* restrict coefs3zs = coefs+3*zs;     ASSUME_ALIGNED(coefs3zs);

            const T pre20 = d2a[i]*  b[j];
            const T pre10 =  da[i]*  b[j];
            const T pre00 =   a[i]*  b[j];
            const T pre11 =  da[i]* db[j];
            const T pre01 =   a[i]* db[j];
            const T pre02 =   a[i]*d2b[j];

            const T sum0 =   c[0] * coefs[n] +   c[1] * coefszs[n] +   c[2] * coefs2zs[n] +   c[3] * coefs3zs[n];
            const T sum1 =  dc[0] * coefs[n] +  dc[1] * coefszs[n] +  dc[2] * coefs2zs[n] +  dc[3] * coefs3zs[n];
            const T sum2 = d2c[0] * coefs[n] + d2c[1] * coefszs[n] + d2c[2] * coefs2zs[n] + d2c[3] * coefs3zs[n];

            hxx_local += pre20 * sum0;
            hxy_local += pre11 * sum0;
            hxz_local += pre10 * sum1;
            hyy_local += pre02 * sum0;
            hyz_local += pre01 * sum1;
            hzz_local += pre00 * sum2;
            gx_local += pre10 * sum0;
            gy_local += pre01 * sum0;
            gz_local += pre00 * sum1;
            val_local += pre00 * sum0;
          }

        hxx[n] = hxx_local*dxx;
        hxy[n] = hxy_local*dxy;
        hxz[n] = hxz_local*dxz;
        hyy[n] = hyy_local*dyy;
        hyz[n] = hyz_local*dyz;
        hzz[n] = hzz_local*dzz;
        gx[n] = gx_local*dxInv;
        gy[n] = gy_local*dyInv;
        gz[n] = gz_local*dzInv;
        vals[n] = val_local;
      }
    }
}/** qmcplusplus namespace */
#endif

