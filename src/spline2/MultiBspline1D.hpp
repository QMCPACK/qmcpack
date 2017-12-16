//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file MultiBspline1D.hpp
 *
 * Literal port of einspline/multi_bspline_eval_d_std3.cpp by Ye at ANL
 * Modified to handle both float/double operations
 * - template parameter T for the precision
 * - MUB spline object created by einspline allocators.
 * Function signatures modified anticipating its use by a class that can perform data parallel execution
 * - evaluate(...., int first, int last)
 */
#ifndef QMCPLUSPLUS_MULTIEINSPLINE_1D_ENGINE_HPP
#define QMCPLUSPLUS_MULTIEINSPLINE_1D_ENGINE_HPP

namespace qmcplusplus
{

  template<typename T>
    inline void 
    MultiBspline1D<T>::evaluate_v_impl(T x, T* restrict vals) const
    {
      x -= spline_m.x_grid.start;
      T tx; int ix;
      SplineBound<T>::get(x*spline_m.x_grid.delta_inv,tx,ix,spline_m.x_grid.num-2);

      T a[4];
      MultiBsplineData<T>::compute_prefactors(a, tx);

      const intptr_t xs = spline_m.x_stride;
      const T dxInv = spline_m.x_grid.delta_inv;

      const T* restrict coefs0 = spline_m.coefs + ((ix  )*xs);
      const T* restrict coefs1 = spline_m.coefs + ((ix+1)*xs);
      const T* restrict coefs2 = spline_m.coefs + ((ix+2)*xs);
      const T* restrict coefs3 = spline_m.coefs + ((ix+3)*xs);

      #pragma omp simd aligned(vals,coefs0,coefs1,coefs2,coefs3)
      for (int n=0; n<spline_m.num_splines; n++)
        vals[n]  = a[0] * coefs0[n] + a[1] * coefs1[n] + a[2] * coefs2[n] + a[3] * coefs3[n];
    }

  template<typename T>
    inline void 
    MultiBspline1D<T>::evaluate_vgl_impl(T x, T* restrict vals, T* restrict grads, T* restrict lapl) const
    {
      x -= spline_m.x_grid.start;
      T tx; int ix;
      SplineBound<T>::get(x*spline_m.x_grid.delta_inv,tx,ix,spline_m.x_grid.num-2);

      T a[4], da[4], d2a[4];
      MultiBsplineData<T>::compute_prefactors(a, da, d2a, tx);

      const intptr_t xs = spline_m.x_stride;
      const T dxInv = spline_m.x_grid.delta_inv;

      const T* restrict coefs0 = spline_m.coefs + ((ix  )*xs);
      const T* restrict coefs1 = spline_m.coefs + ((ix+1)*xs);
      const T* restrict coefs2 = spline_m.coefs + ((ix+2)*xs);
      const T* restrict coefs3 = spline_m.coefs + ((ix+3)*xs);

      #pragma omp simd aligned(vals,grads,lapl,coefs0,coefs1,coefs2,coefs3)
      for (int n=0; n<spline_m.num_splines; n++)
      {
        const T coef_0=coefs0[n];
        const T coef_1=coefs1[n];
        const T coef_2=coefs2[n];
        const T coef_3=coefs3[n];
        vals[n]  =    a[0] * coef_0 +   a[1] * coef_1 +   a[2] * coef_2 +   a[3] * coef_3;
        grads[n] = ( da[0] * coef_0 +  da[1] * coef_1 +  da[2] * coef_2 +  da[3] * coef_3)*dxInv;
        lapl[n]  = (d2a[0] * coef_0 + d2a[1] * coef_1 + d2a[2] * coef_2 + d2a[3] * coef_3)*dxInv*dxInv;
      }
    }

}/** qmcplusplus namespace */
#endif

