//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file MultiBsplineStd2.hpp
 *
 * Literal port of einspline/multi_bspline_eval_d_std2.cpp by A. Benali at ANL
 * Modified to handle both float/double and blocked operations
 * - template parameter T for the precision
 * - MUB spline object created by einspline allocators.
 * Function signatures modified anticipating its use by a class that can perform data parallel execution
 * - evaluate(...., int first, int last)
 */
#ifndef QMCPLUSPLUS_MULTIEINSPLINE_STD2_ENGINE_HPP
#define QMCPLUSPLUS_MULTIEINSPLINE_STD2_ENGINE_HPP
//#include <spline/einspline_engine.hpp>
namespace qmcplusplus
{
#if 0
  template<typename T>
    template<typename MUB, typename PT>
    inline void 
    MultiBspline<T>::evaluate(const MUB* spline, PT x0, PT y0, PT z0, T* restrict vals)
    {
      T x = static_cast<T>(x0)-spline->x_grid.start;
      T y = static_cast<T>(y0)-spline->y_grid.start;
      T z = static_cast<T>(z0)-spline->z_grid.start;
      T tx,ty,tz;
      int ix,iy,iz;
      SplineBound<T>::get(x*spline->x_grid.delta_inv,tx,ix,spline->x_grid.num-1);
      SplineBound<T>::get(y*spline->y_grid.delta_inv,ty,iy,spline->y_grid.num-1);
      SplineBound<T>::get(z*spline->z_grid.delta_inv,tz,iz,spline->z_grid.num-1);
      T a[4], b[4], c[4];
      a[0] = ( ( A44[0]  * tx + A44[1] )  * tx + A44[2] )  * tx + A44[3]; 
      a[1] = ( ( A44[4]  * tx + A44[5] )  * tx + A44[6] )  * tx + A44[7]; 
      a[2] = ( ( A44[8]  * tx + A44[9] )  * tx + A44[10] ) * tx + A44[11]; 
      a[3] = ( ( A44[12] * tx + A44[13] ) * tx + A44[14] ) * tx + A44[15]; 

      b[0] = ( ( A44[0]  * ty + A44[1] )    * ty + A44[2] )  * ty + A44[3]; 
      b[1] = ( ( A44[4]  * ty + A44[5] )    * ty + A44[6] )  * ty + A44[7]; 
      b[2] = ( ( A44[8]  * ty + A44[9] )    * ty + A44[10] ) * ty + A44[11]; 
      b[3] = ( ( A44[12] * ty + A44[13] )   * ty + A44[14] ) * ty + A44[15]; 

      c[0] = ( ( A44[0]  * tz + A44[1] )  * tz + A44[2] )  * tz + A44[3]; 
      c[1] = ( ( A44[4]  * tz + A44[5] )  * tz + A44[6] )  * tz + A44[7]; 
      c[2] = ( ( A44[8]  * tz + A44[9] )  * tz + A44[10] ) * tz + A44[11]; 
      c[3] = ( ( A44[12] * tz + A44[13] ) * tz + A44[14] ) * tz + A44[15]; 
      intptr_t xs = spline->x_stride;
      intptr_t ys = spline->y_stride;
      intptr_t zs = spline->z_stride;

      T d[64];
      intptr_t mod_coefs[64];
      __assume_aligned(vals,64);

      {//precompute value & offsets
        int n=0;
        for ( int i=0; i<4; i++)
          for ( int j=0; j<4; j++) 
            for ( int k=0; k<4; k++) 
            {
              d[n] = a[i]*b[j]*c[k];
              mod_coefs[n] = ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
              n++;
            }
        std::fill(vals,vals+spline->num_splines,T());
      }

      //#pragma unroll(16)
      for(int i=0; i<64; ++i)
      {
        T s = d[i];
        const T* restrict p = spline->coefs+mod_coefs[i];
        __assume_aligned(p,64);
#pragma omp simd
        for(int n=0; n<spline->num_splines; ++n)
          vals[n]+=s*p[n];
      }
    }
#endif
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

      MultiBsplineData<T>::compute_prefactors(a, da, d2a, tx);
      MultiBsplineData<T>::compute_prefactors(b, db, d2b, ty);
      MultiBsplineData<T>::compute_prefactors(c, dc, d2c, tz);

      const intptr_t xs = spline_m->x_stride;
      const intptr_t ys = spline_m->y_stride;
      const intptr_t zs = spline_m->z_stride;
      const intptr_t offset=ix*xs+iy*ys+iz*zs+first;

      out_offset=(out_offset)?out_offset:spline_m->num_splines;
      const int num_splines=last-first;

      ASSUME_ALIGNED(vals);
      T* restrict gx=grads;              ASSUME_ALIGNED(gx);
      T* restrict gy=grads+  out_offset; ASSUME_ALIGNED(gy);
      T* restrict gz=grads+2*out_offset; ASSUME_ALIGNED(gz);
      T* restrict lx=lapl;               ASSUME_ALIGNED(lx);

      const T dxInv = spline_m->x_grid.delta_inv;
      const T dyInv = spline_m->y_grid.delta_inv;
      const T dzInv = spline_m->z_grid.delta_inv; 
      const T dxxInv = dxInv*dxInv;
      const T dxyInv = dxInv*dyInv;
      const T dxzInv = dxInv*dzInv;
      const T dyyInv = dyInv*dyInv;
      const T dyzInv = dyInv*dzInv;
      const T dzzInv = dzInv*dzInv;
      //const T* restrict coefs = spline->coefs + ix*xs + iy*ys + iz*zs;
//#pragma omp parallel for simd
#pragma omp simd
      for (int n=0; n<num_splines; n++) {
        const T* restrict coefs = spline_m->coefs + offset + n;
        T val = T(), grad0 = T(), grad1 = T(), grad2 = T();
        T lap0 = T(), lap1 = T(), lap2 = T(), lap3=T();

#pragma unroll(4)
        for (int i=0; i<4; i++) {
          const T pre0 =   a[i] *   b[0];
          const T pre1 =  da[i] *   b[0];
          const T pre2 = d2a[i] *   b[0];
          const T pre3 =   a[i] *  db[0];
          const T coef0 = coefs[i*xs];
          const T coef1 = coefs[i*xs + zs];
          const T coef2 = coefs[i*xs + 2*zs];
          const T coef3 = coefs[i*xs + 3*zs];
          const T sum0 = c[0] * coef0 + c[1] * coef1 + c[2] * coef2 + c[3] * coef3;
          const T sum1 = dc[0] * coef0 + dc[1] * coef1 + dc[2] * coef2 + dc[3] * coef3;

          const T pre01 =   a[i] *   b[1];
          const T pre11 =  da[i] *   b[1];
          const T pre21 = d2a[i] *   b[1];
          const T pre31 =   a[i] *  db[1];
          const T coef01 = coefs[i*xs + ys];
          const T coef11 = coefs[i*xs + ys + zs];
          const T coef21 = coefs[i*xs + ys + 2*zs];
          const T coef31 = coefs[i*xs + ys + 3*zs];
          const T sum01 = c[0] * coef01 + c[1] * coef11 + c[2] * coef21 + c[3] * coef31;
          const T sum11 = dc[0] * coef01 + dc[1] * coef11 + dc[2] * coef21 + dc[3] * coef31;

          const T pre02 =   a[i] *   b[2];
          const T pre12 =  da[i] *   b[2];
          const T pre22 = d2a[i] *   b[2];
          const T pre32 =   a[i] *  db[2];
          const T coef02 = coefs[i*xs + 2*ys];
          const T coef12 = coefs[i*xs + 2*ys + zs];
          const T coef22 = coefs[i*xs + 2*ys + 2*zs];
          const T coef32 = coefs[i*xs + 2*ys + 3*zs];
          const T sum02 = c[0] * coef02 + c[1] * coef12 + c[2] * coef22 + c[3] * coef32;
          const T sum12 = dc[0] * coef02 + dc[1] * coef12 + dc[2] * coef22 + dc[3] * coef32;

          const T pre03 =   a[i] *   b[3];
          const T pre13 =  da[i] *   b[3];
          const T pre23 = d2a[i] *   b[3];
          const T pre33 =   a[i] *  db[3];
          const T coef03 = coefs[i*xs + 3*ys];
          const T coef13 = coefs[i*xs + 3*ys + zs];
          const T coef23 = coefs[i*xs + 3*ys + 2*zs];
          const T coef33 = coefs[i*xs + 3*ys + 3*zs];
          const T sum03 = c[0] * coef03 + c[1] * coef13 + c[2] * coef23 + c[3] * coef33;
          const T sum13 = dc[0] * coef03 + dc[1] * coef13 + dc[2] * coef23 + dc[3] * coef33;

          val   +=   pre0 * sum0 + pre01 * sum01 + pre02 * sum02 + pre03 * sum03;
          lap2  += pre0 * (d2c[0] * coef0 +d2c[1] * coef1 + d2c[2] * coef2  + d2c[3] * coef3) + 
            pre01 * (d2c[0] * coef01 +d2c[1] * coef11 + d2c[2] * coef21  + d2c[3] * coef31) +
            pre02 * (d2c[0] * coef02 +d2c[1] * coef12 + d2c[2] * coef22  + d2c[3] * coef32) +
            pre03 * (d2c[0] * coef03 +d2c[1] * coef13 + d2c[2] * coef23  + d2c[3] * coef33);
          lap1 += (a[i] * d2b[0]) * sum0 + (a[i] * d2b[1]) * sum01 + (a[i] * d2b[2]) * sum02 + (a[i] * d2b[3]) * sum03;
          grad0 +=   pre1 * sum0 + pre11 * sum01 + pre12 * sum02 + pre13 * sum03;
          grad1 +=   pre3 * sum0 + pre33 * sum03 + pre31 * sum01 + pre32 * sum02;
          grad2 +=   pre0 * sum1 + pre01 * sum11 + pre02 * sum12 + pre03 * sum13;
          lap0  +=   pre2 * sum0 + pre21 * sum01 + pre22 * sum02 + pre23 * sum03;
        }

        vals[n] = val;
        gx[n] = grad0*dxInv;
        gy[n] = grad1*dyInv;
        gz[n] = grad2*dzInv;
        lx[n] = (lap0*dxxInv+lap1*dyyInv+lap2*dzzInv);
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
      const intptr_t offset=ix*xs+iy*ys+iz*zs+first;

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
      const T dxxInv = dxInv*dxInv;
      const T dxyInv = dxInv*dyInv;
      const T dxzInv = dxInv*dzInv;
      const T dyyInv = dyInv*dyInv;
      const T dyzInv = dyInv*dzInv;
      const T dzzInv = dzInv*dzInv;

      //#pragma omp parallel for simd
#pragma omp simd
      for (int n=0; n<num_splines; n++) {
        const T* restrict coefs = spline_m->coefs + offset+n;
        T val =T() , grad0 = T(), grad1 = T(), grad2 = T();
        T hess0 = T(), hess1 = T(), hess2 = T(), hess4 = T(), hess5 = T(), hess8 = T();
#pragma unroll(4)
        for (int i=0; i<4; i++) {

          const T pre0 =   a[i] *   b[0];
          const T pre1 =  da[i] *   b[0];
          const T pre2 = d2a[i] *   b[0];
          const T pre3 =   a[i] *  db[0];
          const T coef0 = coefs[i*xs];
          const T coef1 = coefs[i*xs + zs];
          const T coef2 = coefs[i*xs + 2*zs];
          const T coef3 = coefs[i*xs + 3*zs];
          const T sum0 = c[0] * coef0 + c[1] * coef1 + c[2] * coef2 + c[3] * coef3;
          const T sum1 = dc[0] * coef0 + dc[1] * coef1 + dc[2] * coef2 + dc[3] * coef3;

          const T pre01 =   a[i] *   b[1];
          const T pre11 =  da[i] *   b[1];
          const T pre21 = d2a[i] *   b[1];
          const T pre31 =   a[i] *  db[1];
          const T coef01 = coefs[i*xs + ys];
          const T coef11 = coefs[i*xs + ys + zs];
          const T coef21 = coefs[i*xs + ys + 2*zs];
          const T coef31 = coefs[i*xs + ys + 3*zs];
          const T sum01 = c[0] * coef01 + c[1] * coef11 + c[2] * coef21 + c[3] * coef31;
          const T sum11 = dc[0] * coef01 + dc[1] * coef11 + dc[2] * coef21 + dc[3] * coef31;

          const T pre02 =   a[i] *   b[2];
          const T pre12 =  da[i] *   b[2];
          const T pre22 = d2a[i] *   b[2];
          const T pre32 =   a[i] *  db[2];
          const T coef02 = coefs[i*xs + 2*ys];
          const T coef12 = coefs[i*xs + 2*ys + zs];
          const T coef22 = coefs[i*xs + 2*ys + 2*zs];
          const T coef32 = coefs[i*xs + 2*ys + 3*zs];
          const T sum02 = c[0] * coef02 + c[1] * coef12 + c[2] * coef22 + c[3] * coef32;
          const T sum12 = dc[0] * coef02 + dc[1] * coef12 + dc[2] * coef22 + dc[3] * coef32;

          const T pre03 =   a[i] *   b[3];
          const T pre13 =  da[i] *   b[3];
          const T pre23 = d2a[i] *   b[3];
          const T pre33 =   a[i] *  db[3];
          const T coef03 = coefs[i*xs + 3*ys];
          const T coef13 = coefs[i*xs + 3*ys + zs];
          const T coef23 = coefs[i*xs + 3*ys + 2*zs];
          const T coef33 = coefs[i*xs + 3*ys + 3*zs];
          const T sum03 = c[0] * coef03 + c[1] * coef13 + c[2] * coef23 + c[3] * coef33;
          const T sum13 = dc[0] * coef03 + dc[1] * coef13 + dc[2] * coef23 + dc[3] * coef33;

          val   +=   pre0 * sum0 + pre01 * sum01 + pre02 * sum02 + pre03 * sum03;
          hess8 += pre0 * (d2c[0] * coef0 +d2c[1] * coef1 + d2c[2] * coef2  + d2c[3] * coef3) + 
            pre01 * (d2c[0] * coef01 +d2c[1] * coef11 + d2c[2] * coef21  + d2c[3] * coef31) +
            pre02 * (d2c[0] * coef02 +d2c[1] * coef12 + d2c[2] * coef22  + d2c[3] * coef32) +
            pre03 * (d2c[0] * coef03 +d2c[1] * coef13 + d2c[2] * coef23  + d2c[3] * coef33);
          hess1 += (da[i] *  db[0]) * sum0 + (da[i] *  db[1]) * sum01 + (da[i] *  db[2]) * sum02 + (da[i] *  db[3]) * sum03;
          hess4 += (a[i] * d2b[0]) * sum0 + (a[i] * d2b[1]) * sum01 + (a[i] * d2b[2]) * sum02 + (a[i] * d2b[3]) * sum03;
          grad0 +=   pre1 * sum0 + pre11 * sum01 + pre12 * sum02 + pre13 * sum03;
          grad1 +=   pre3 * sum0 + pre33 * sum03 + pre31 * sum01 + pre32 * sum02;
          grad2 +=   pre0 * sum1 + pre01 * sum11 + pre02 * sum12 + pre03 * sum13;
          hess2 +=   pre1 * sum1 + pre11 * sum11 + pre12 * sum12 + pre13 * sum13;
          hess5 +=   pre3 * sum1 + pre31 * sum11 + pre32 * sum12 + pre33 * sum13;
          hess0 +=   pre2 * sum0 + pre21 * sum01 + pre22 * sum02 + pre23 * sum03;

        }

        vals[n] = val;
        gx[n] = grad0*dxInv;
        gy[n] = grad1*dyInv;
        gz[n] = grad2*dzInv;

        hxx[n]=hess0*dxxInv;
        hxy[n]=hess1*dxyInv;
        hxz[n]=hess2*dxzInv;
        hyy[n]=hess4*dxxInv;
        hyz[n]=hess5*dxyInv;
        hzz[n]=hess8*dxzInv;
      }
    }
}/** qmcplusplus namespace */
#endif

