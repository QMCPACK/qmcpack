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
#ifndef QMCPLUSPLUS_MULTIEINSPLINE_VALUE_HPP
#define QMCPLUSPLUS_MULTIEINSPLINE_VALUE_HPP

#if defined(__xlC__) && defined(BGQPX)
#include <builtins.h>
#endif

namespace qmcplusplus
{

#ifndef BGQPX
  /** define evaluate: common to any implementation */
  template<typename T>
    inline void
    MultiBspline<T>::evaluate_v_impl(T x, T y, T z, T* restrict vals, int first, int last) const
    {
      x -= spline_m->x_grid.start;
      y -= spline_m->y_grid.start;
      z -= spline_m->z_grid.start;
      T tx,ty,tz;
      int ix,iy,iz;
      SplineBound<T>::get(x*spline_m->x_grid.delta_inv,tx,ix,spline_m->x_grid.num-1);
      SplineBound<T>::get(y*spline_m->y_grid.delta_inv,ty,iy,spline_m->y_grid.num-1);
      SplineBound<T>::get(z*spline_m->z_grid.delta_inv,tz,iz,spline_m->z_grid.num-1);
      T a[4], b[4], c[4];

      MultiBsplineData<T>::compute_prefactors(a, tx);
      MultiBsplineData<T>::compute_prefactors(b, ty);
      MultiBsplineData<T>::compute_prefactors(c, tz);

      const intptr_t xs = spline_m->x_stride;
      const intptr_t ys = spline_m->y_stride;
      const intptr_t zs = spline_m->z_stride;

      CONSTEXPR T zero(0);
      const int num_splines=last-first;
      std::fill(vals,vals+num_splines,zero);

      for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
        {
          const T pre00 =  a[i]*b[j];
          const T* restrict coefs = spline_m->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs) + first;
          const T* restrict coefszs  = coefs+zs;
          const T* restrict coefs2zs = coefs+2*zs;
          const T* restrict coefs3zs = coefs+3*zs;
          #pragma omp simd aligned(coefs,coefszs,coefs2zs,coefs3zs,vals)
          for(int n=0; n<num_splines; n++)
            vals[n] += pre00*(c[0]*coefs[n] + c[1]*coefszs[n] + c[2]*coefs2zs[n] + c[3]*coefs3zs[n]);
        }
    }
#else
// this is only experimental, not protected for general use.
  template<typename T>
    inline void
    MultiBspline<T>::evaluate_v_impl(T x, T y, T z, T* restrict vals, int first, int last) const
    {
      x -= spline_m->x_grid.start;
      y -= spline_m->y_grid.start;
      z -= spline_m->z_grid.start;
      T tx,ty,tz;
      int ix,iy,iz;
      SplineBound<T>::get(x*spline_m->x_grid.delta_inv,tx,ix,spline_m->x_grid.num-1);
      SplineBound<T>::get(y*spline_m->y_grid.delta_inv,ty,iy,spline_m->y_grid.num-1);
      SplineBound<T>::get(z*spline_m->z_grid.delta_inv,tz,iz,spline_m->z_grid.num-1);
      T a[4], b[4], c[4];

      MultiBsplineData<T>::compute_prefactors(a, tx);
      MultiBsplineData<T>::compute_prefactors(b, ty);
      MultiBsplineData<T>::compute_prefactors(c, tz);

      vector4double vec_c0 = vec_splats(c[0]);
      vector4double vec_c1 = vec_splats(c[1]);
      vector4double vec_c2 = vec_splats(c[2]);
      vector4double vec_c3 = vec_splats(c[3]);

      const intptr_t xs = spline_m->x_stride;
      const intptr_t ys = spline_m->y_stride;
      const intptr_t zs = spline_m->z_stride;

      CONSTEXPR T zero(0);
      const int num_splines=last-first;
      std::fill(vals,vals+num_splines,zero);

      for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
        {
          const T pre00 =  a[i]*b[j];
          vector4double vec_pre00 = vec_splats(pre00);
          T* restrict coefs0 = spline_m->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs) + first;
          T* restrict coefs1 = coefs0 +   zs;
          T* restrict coefs2 = coefs0 + 2*zs;
          T* restrict coefs3 = coefs0 + 3*zs;
          for(int n=0, p=0; n<num_splines; n+=4, p+=4*sizeof(T))
          {
            vector4double vec_coef0, vec_coef1, vec_coef2, vec_coef3, vec_val;

            __dcbt(&coefs0[n+8]);
            __dcbt(&coefs1[n+8]);
            __dcbt(&coefs2[n+8]);
            __dcbt(&coefs3[n+8]);

            vec_coef0 = vec_ld(p, coefs0);
            vec_coef1 = vec_ld(p, coefs1);
            vec_coef2 = vec_ld(p, coefs2);
            vec_coef3 = vec_ld(p, coefs3);
            vec_val = vec_ld(p, vals);

            vec_coef0 = vec_mul(vec_c0, vec_coef0);
            vec_coef0 = vec_madd(vec_c1, vec_coef1, vec_coef0);
            vec_coef0 = vec_madd(vec_c2, vec_coef2, vec_coef0);
            vec_coef0 = vec_madd(vec_c3, vec_coef3, vec_coef0);
            vec_val = vec_madd(vec_pre00, vec_coef0, vec_val);
            vec_st(vec_val, p, vals);
          }
        }
    }
#endif

}
#endif
