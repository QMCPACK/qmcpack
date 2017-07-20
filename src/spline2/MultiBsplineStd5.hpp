//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file MultiBsplineStd5.hpp
 *
 * Experimental blocking version based on MultiBsplineStd.hpp
 */
#ifndef QMCPLUSPLUS_MULTIEINSPLINE_STD5_ENGINE_HPP
#define QMCPLUSPLUS_MULTIEINSPLINE_STD5_ENGINE_HPP

namespace qmcplusplus
{

#ifndef BGQPX

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
      T* restrict lx=lapl;               ASSUME_ALIGNED(lx);
      T* restrict ly=lapl+  out_offset;  ASSUME_ALIGNED(ly);
      T* restrict lz=lapl+2*out_offset;  ASSUME_ALIGNED(lz);

      std::fill(vals,vals+num_splines,T());
      std::fill(gx,gx+num_splines,T());
      std::fill(gy,gy+num_splines,T());
      std::fill(gz,gz+num_splines,T());
      std::fill(lx,lx+num_splines,T());
      std::fill(ly,ly+num_splines,T());
      std::fill(lz,lz+num_splines,T());

      for (int i=0; i<4; i++)
        for (int j=0; j<4; j++){

          const T pre20 = d2a[i]*  b[j];
          const T pre10 =  da[i]*  b[j];
          const T pre00 =   a[i]*  b[j];
          const T pre11 =  da[i]* db[j];
          const T pre01 =   a[i]* db[j];
          const T pre02 =   a[i]*d2b[j];

          const T* restrict coefs = spline_m->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs) + first; ASSUME_ALIGNED(coefs);
          const T* restrict coefszs  = coefs+zs;       ASSUME_ALIGNED(coefszs);
          const T* restrict coefs2zs = coefs+2*zs;     ASSUME_ALIGNED(coefs2zs);
          const T* restrict coefs3zs = coefs+3*zs;     ASSUME_ALIGNED(coefs3zs);

          #pragma noprefetch
          #pragma omp simd
          for (int n=0; n<num_splines; n++) {
            const T coefsv = coefs[n];
            const T coefsvzs = coefszs[n];
            const T coefsv2zs = coefs2zs[n];
            const T coefsv3zs = coefs3zs[n];

            T sum0 =   c[0] * coefsv +   c[1] * coefsvzs +   c[2] * coefsv2zs +   c[3] * coefsv3zs;
            T sum1 =  dc[0] * coefsv +  dc[1] * coefsvzs +  dc[2] * coefsv2zs +  dc[3] * coefsv3zs;
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

      const T dxInv2 = dxInv*dxInv;
      const T dyInv2 = dyInv*dyInv;
      const T dzInv2 = dzInv*dzInv;

      #pragma omp simd
      for (int n=0; n<num_splines; n++) 
      {
        gx[n] *= dxInv;
        gy[n] *= dyInv;
        gz[n] *= dzInv;
        lx[n] = lx[n]*dxInv2+ly[n]*dyInv2+lz[n]*dzInv2;
      }
    }

  template<typename T>
    inline void 
    MultiBspline<T>::evaluate_vgh_impl(T x, T y, T z, 
        T* restrict vals, T* restrict grads, T* restrict hess, int first, int last, size_t out_offset) const
    {

      int ix,iy,iz;
      T tx,ty,tz;
      T a[4],b[4],c[4],da[4],db[4],dc[4],d2a[4],d2b[4],d2c[4];

      x -= spline_m->x_grid.start;
      y -= spline_m->y_grid.start;
      z -= spline_m->z_grid.start;
      SplineBound<T>::get(x*spline_m->x_grid.delta_inv,tx,ix,spline_m->x_grid.num-1);
      SplineBound<T>::get(y*spline_m->y_grid.delta_inv,ty,iy,spline_m->y_grid.num-1);
      SplineBound<T>::get(z*spline_m->z_grid.delta_inv,tz,iz,spline_m->z_grid.num-1);

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

      std::fill(vals,vals+num_splines,T());
      std::fill(gx,gx+num_splines,T());
      std::fill(gy,gy+num_splines,T());
      std::fill(gz,gz+num_splines,T());
      std::fill(hxx,hxx+num_splines,T());
      std::fill(hxy,hxy+num_splines,T());
      std::fill(hxz,hxz+num_splines,T());
      std::fill(hyy,hyy+num_splines,T());
      std::fill(hyz,hyz+num_splines,T());
      std::fill(hzz,hzz+num_splines,T());

      const int iBlockSize = 64; 
      const int nBlocks = ((num_splines+iBlockSize-1)/iBlockSize); 
      int iLastBlockSize = (num_splines - nBlocks*iBlockSize );
      if ( iLastBlockSize == 0 ) {
        iLastBlockSize = iBlockSize; 
      }

      for ( int p = 0; p < nBlocks; p++ ) {

        const int iStartSpline = p*iBlockSize;
        const int iCurrBlockSize = (p==(nBlocks-1)) ? iLastBlockSize : iBlockSize;

        for (int i=0; i<4; i++) {
          for (int j=0; j<4; j++){
            const T* restrict coefs = spline_m->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs) + first; ASSUME_ALIGNED(coefs);
            const T* restrict coefszs  = coefs+zs + iStartSpline;       ASSUME_ALIGNED(coefszs);
            const T* restrict coefs2zs = coefs+2*zs + iStartSpline;     ASSUME_ALIGNED(coefs2zs);
            const T* restrict coefs3zs = coefs+3*zs + iStartSpline;     ASSUME_ALIGNED(coefs3zs);

            const T pre20 = d2a[i]*  b[j];
            const T pre10 =  da[i]*  b[j];
            const T pre00 =   a[i]*  b[j];
            const T pre11 =  da[i]* db[j];
            const T pre01 =   a[i]* db[j];
            const T pre02 =   a[i]*d2b[j];

#if defined(__AVX512F__) && defined(QMC_PREFETCH)
            const int iSplitPoint = std::min( iCurrBlockSize, (iBlockSize - 32));
#else
            const int iSplitPoint = iCurrBlockSize;
#endif
            //lamda function
            auto lamdaF = [&]( const int iFirst, const int iLast ) { 

              #pragma omp simd
              for (int n=iFirst; n<iLast; n++) {
            
                const int n1 = n+iStartSpline;
                T coefsv = coefs[n];
                T coefsvzs = coefszs[n];
                T coefsv2zs = coefs2zs[n];
                T coefsv3zs = coefs3zs[n];

                T sum0 =   c[0] * coefsv +   c[1] * coefsvzs +   c[2] * coefsv2zs +   c[3] * coefsv3zs;
                T sum1 =  dc[0] * coefsv +  dc[1] * coefsvzs +  dc[2] * coefsv2zs +  dc[3] * coefsv3zs;
                T sum2 = d2c[0] * coefsv + d2c[1] * coefsvzs + d2c[2] * coefsv2zs + d2c[3] * coefsv3zs;

                hxx[n1] += pre20 * sum0;
                hxy[n1] += pre11 * sum0;
                hxz[n1] += pre10 * sum1;
                hyy[n1] += pre02 * sum0;
                hyz[n1] += pre01 * sum1;
                hzz[n1] += pre00 * sum2;
                gx[n1] += pre10 * sum0;
                gy[n1] += pre01 * sum0;
                gz[n1] += pre00 * sum1;
                vals[n1]+= pre00 * sum0;
             }
           };

           lamdaF ( 0, iSplitPoint); 
#if defined(__AVX512F__) && defined(QMC_PREFETCH)
           {
             int pfi = (j==3) ? i+1:i;
             int pfj = (j+1)%4;
             T* restrict coefs = spline_m->coefs + ((ix+pfi)*xs + (iy+pfj)*ys + iz*zs) + first; ASSUME_ALIGNED(coefs);
             T* restrict coefszs  = coefs+zs + iStartSpline;       ASSUME_ALIGNED(coefszs);
             T* restrict coefs2zs = coefs+2*zs + iStartSpline;     ASSUME_ALIGNED(coefs2zs);
             T* restrict coefs3zs = coefs+3*zs + iStartSpline;     ASSUME_ALIGNED(coefs3zs);

             for ( int dist = 0; dist < 64; dist +=16 ) {
               _mm_prefetch((char const*)(coefs+dist),_MM_HINT_T1);
               _mm_prefetch((char const*)(coefszs+dist),_MM_HINT_T1);
               _mm_prefetch((char const*)(coefs2zs+dist),_MM_HINT_T1);
               _mm_prefetch((char const*)(coefs3zs+dist),_MM_HINT_T1);
             }
           }

           lamdaF ( iSplitPoint, iCurrBlockSize ); 
#endif
          }
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

      for ( int p = 0; p < nBlocks; p++ ) {
        const int iStartSpline = p*iBlockSize; 
        const int iCurrBlockSize = (p==(nBlocks-1)) ? iLastBlockSize : iBlockSize;
        #pragma omp simd
        for (int n=0; n<iCurrBlockSize; n++)
        {
          const int n1 = n+iStartSpline;  
          gx[n1]*=dxInv; 
          gy[n1]*=dyInv; 
          gz[n1]*=dzInv; 
          hxx[n1]*=dxx;
          hyy[n1]*=dyy;
          hzz[n1]*=dzz;
          hxy[n1]*=dxy;
          hxz[n1]*=dxz;
          hyz[n1]*=dyz;
        }
      }  
    }
#else

// this is only experimental, not protected for general use.
#ifdef __xlC__
#include <builtins.h>
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

      compute_prefactors(a, da, d2a, tx);
      compute_prefactors(b, db, d2b, ty);
      compute_prefactors(c, dc, d2c, tz);

      vector4double vec_c0 = vec_splats(c[0]);
      vector4double vec_c1 = vec_splats(c[1]);
      vector4double vec_c2 = vec_splats(c[2]);
      vector4double vec_c3 = vec_splats(c[3]);
      vector4double vec_dc0 = vec_splats(dc[0]);
      vector4double vec_dc1 = vec_splats(dc[1]);
      vector4double vec_dc2 = vec_splats(dc[2]);
      vector4double vec_dc3 = vec_splats(dc[3]);
      vector4double vec_d2c0 = vec_splats(d2c[0]);
      vector4double vec_d2c1 = vec_splats(d2c[1]);
      vector4double vec_d2c2 = vec_splats(d2c[2]);
      vector4double vec_d2c3 = vec_splats(d2c[3]);

      const intptr_t xs = spline_m->x_stride;
      const intptr_t ys = spline_m->y_stride;
      const intptr_t zs = spline_m->z_stride;

      out_offset=(out_offset)?out_offset:spline_m->num_splines;
      const int num_splines=last-first;

      ASSUME_ALIGNED(vals);
      T* restrict gx=grads;              ASSUME_ALIGNED(gx);
      T* restrict gy=grads+  out_offset; ASSUME_ALIGNED(gy);
      T* restrict gz=grads+2*out_offset; ASSUME_ALIGNED(gz);
      T* restrict lx=lapl;               ASSUME_ALIGNED(lx);
      T* restrict ly=lapl+  out_offset;  ASSUME_ALIGNED(ly);
      T* restrict lz=lapl+2*out_offset;  ASSUME_ALIGNED(lz);

      std::fill(vals,vals+num_splines,T());
      std::fill(gx,gx+num_splines,T());
      std::fill(gy,gy+num_splines,T());
      std::fill(gz,gz+num_splines,T());
      std::fill(lx,lx+num_splines,T());
      std::fill(ly,ly+num_splines,T());
      std::fill(lz,lz+num_splines,T());

      int n = 0;

      for (int i=0; i<4; i++)
        for (int j=0; j<4; j++){
          T* restrict coefs0 = spline_m->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs) + first;
          T* restrict coefs1 = coefs0 + zs;
          T* restrict coefs2 = coefs0 + 2*zs;
          T* restrict coefs3 = coefs0 + 3*zs;

          const T pre20 = d2a[i]*  b[j];
          const T pre10 =  da[i]*  b[j];
          const T pre00 =   a[i]*  b[j];
          const T pre01 =   a[i]* db[j];
          const T pre02 =   a[i]*d2b[j];

          vector4double vec_pre00 = vec_splats(pre00);
          vector4double vec_pre01 = vec_splats(pre01);
          vector4double vec_pre02 = vec_splats(pre02);
          vector4double vec_pre10 = vec_splats(pre10);
          vector4double vec_pre20 = vec_splats(pre20);

          n = 0;
          int val_p = 0;
          for (; n<num_splines; n+=4, val_p+=4*sizeof(T)) {

            __dcbt(&coefs0[n+8]);
            __dcbt(&coefs1[n+8]);
            __dcbt(&coefs2[n+8]);
            __dcbt(&coefs3[n+8]);
            __dcbt(&gx    [n+8]);
            __dcbt(&gy    [n+8]);
            __dcbt(&gz    [n+8]);
            //__dcbt(&lx    [n+8]);
            //__dcbt(&ly    [n+8]);
            //__dcbt(&lz    [n+8]);
            __dcbt(&vals  [n+8]);

            vector4double coef0 = vec_ld(0, &coefs0[n]);
            vector4double coef1 = vec_ld(0, &coefs1[n]);
            vector4double coef2 = vec_ld(0, &coefs2[n]);
            vector4double coef3 = vec_ld(0, &coefs3[n]);

            vector4double sum0, sum1, sum2;
            sum0 = vec_mul (vec_c0, coef0);
            sum0 = vec_madd(vec_c1, coef1, sum0);
            sum0 = vec_madd(vec_c2, coef2, sum0);
            sum0 = vec_madd(vec_c3, coef3, sum0);
            sum1 = vec_mul (vec_dc0, coef0);
            sum1 = vec_madd(vec_dc1, coef1, sum1);
            sum1 = vec_madd(vec_dc2, coef2, sum1);
            sum1 = vec_madd(vec_dc3, coef3, sum1);
            sum2 = vec_mul (vec_d2c0, coef0);
            sum2 = vec_madd(vec_d2c1, coef1, sum2);
            sum2 = vec_madd(vec_d2c2, coef2, sum2);
            sum2 = vec_madd(vec_d2c3, coef3, sum2);

            vector4double temp_vec;

            temp_vec = vec_ld(val_p, gx);
            temp_vec = vec_madd(vec_pre10, sum0, temp_vec);
            vec_st(temp_vec, val_p, gx);
            temp_vec = vec_ld(val_p, gy);
            temp_vec = vec_madd(vec_pre01, sum0, temp_vec);
            vec_st(temp_vec, val_p, gy);
            temp_vec = vec_ld(val_p, gz);
            temp_vec = vec_madd(vec_pre00, sum1, temp_vec);
            vec_st(temp_vec, val_p, gz);

            temp_vec = vec_ld(val_p, lx);
            temp_vec = vec_madd(vec_pre20, sum0, temp_vec);
            vec_st(temp_vec, val_p, lx);
            temp_vec = vec_ld(val_p, ly);
            temp_vec = vec_madd(vec_pre02, sum0, temp_vec);
            vec_st(temp_vec, val_p, ly);
            temp_vec = vec_ld(val_p, lz);
            temp_vec = vec_madd(vec_pre00, sum2, temp_vec);
            vec_st(temp_vec, val_p, lz);

            temp_vec = vec_ld(val_p, vals);
            temp_vec = vec_madd(vec_pre00, sum0, temp_vec);
            vec_st(temp_vec, val_p, vals);

          }
        }

      const T dxInv = spline_m->x_grid.delta_inv;
      const T dyInv = spline_m->y_grid.delta_inv;
      const T dzInv = spline_m->z_grid.delta_inv;

      const T dxInv2 = dxInv*dxInv;
      const T dyInv2 = dyInv*dyInv;
      const T dzInv2 = dzInv*dzInv;

      for (int n=0; n<num_splines; n++) 
      {
        gx[n] *= dxInv;
        gy[n] *= dyInv;
        gz[n] *= dzInv;
        lx[n] = lx[n]*dxInv2+ly[n]*dyInv2+lz[n]*dzInv2;
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

      vector4double vec_c0 = vec_splats(c[0]);
      vector4double vec_c1 = vec_splats(c[1]);
      vector4double vec_c2 = vec_splats(c[2]);
      vector4double vec_c3 = vec_splats(c[3]);
      vector4double vec_dc0 = vec_splats(dc[0]);
      vector4double vec_dc1 = vec_splats(dc[1]);
      vector4double vec_dc2 = vec_splats(dc[2]);
      vector4double vec_dc3 = vec_splats(dc[3]);
      vector4double vec_d2c0 = vec_splats(d2c[0]);
      vector4double vec_d2c1 = vec_splats(d2c[1]);
      vector4double vec_d2c2 = vec_splats(d2c[2]);
      vector4double vec_d2c3 = vec_splats(d2c[3]);

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

      std::fill(vals,vals+num_splines,T());
      std::fill(gx,gx+num_splines,T());
      std::fill(gy,gy+num_splines,T());
      std::fill(gz,gz+num_splines,T());
      std::fill(hxx,hxx+num_splines,T());
      std::fill(hxy,hxy+num_splines,T());
      std::fill(hxz,hxz+num_splines,T());
      std::fill(hyy,hyy+num_splines,T());
      std::fill(hyz,hyz+num_splines,T());
      std::fill(hzz,hzz+num_splines,T());

      int n = 0;

      for (int i=0; i<4; i++)
        for (int j=0; j<4; j++){
          T* restrict coefs0 = spline_m->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs) + first;
          T* restrict coefs1 = coefs0 + zs;
          T* restrict coefs2 = coefs0 + 2*zs;
          T* restrict coefs3 = coefs0 + 3*zs;

          const T pre20 = d2a[i]*  b[j];
          const T pre10 =  da[i]*  b[j];
          const T pre00 =   a[i]*  b[j];
          const T pre11 =  da[i]* db[j];
          const T pre01 =   a[i]* db[j];
          const T pre02 =   a[i]*d2b[j];

          vector4double vec_pre00 = vec_splats(pre00);
          vector4double vec_pre01 = vec_splats(pre01);
          vector4double vec_pre02 = vec_splats(pre02);
          vector4double vec_pre11 = vec_splats(pre11);
          vector4double vec_pre10 = vec_splats(pre10);
          vector4double vec_pre20 = vec_splats(pre20);

          n = 0;
          int val_p = 0;
          for (; n<num_splines; n+=4, val_p+=4*sizeof(T)) {

            __dcbt(&coefs0[n+8]);
            __dcbt(&coefs1[n+8]);
            __dcbt(&coefs2[n+8]);
            __dcbt(&coefs3[n+8]);
            __dcbt(&gx    [n+8]);
            __dcbt(&gy    [n+8]);
            __dcbt(&gz    [n+8]);
            //__dcbt(&hxx   [n+8]);
            //__dcbt(&hxy   [n+8]);
            //__dcbt(&hxz   [n+8]);
            //__dcbt(&hyy   [n+8]);
            //__dcbt(&hyz   [n+8]);
            //__dcbt(&hzz   [n+8]);
            __dcbt(&vals  [n+8]);

            vector4double coef0 = vec_ld(0, &coefs0[n]);
            vector4double coef1 = vec_ld(0, &coefs1[n]);
            vector4double coef2 = vec_ld(0, &coefs2[n]);
            vector4double coef3 = vec_ld(0, &coefs3[n]);

            vector4double sum0, sum1, sum2;
            sum0 = vec_mul (vec_c0, coef0);
            sum0 = vec_madd(vec_c1, coef1, sum0);
            sum0 = vec_madd(vec_c2, coef2, sum0);
            sum0 = vec_madd(vec_c3, coef3, sum0);
            sum1 = vec_mul (vec_dc0, coef0);
            sum1 = vec_madd(vec_dc1, coef1, sum1);
            sum1 = vec_madd(vec_dc2, coef2, sum1);
            sum1 = vec_madd(vec_dc3, coef3, sum1);
            sum2 = vec_mul (vec_d2c0, coef0);
            sum2 = vec_madd(vec_d2c1, coef1, sum2);
            sum2 = vec_madd(vec_d2c2, coef2, sum2);
            sum2 = vec_madd(vec_d2c3, coef3, sum2);

            vector4double temp_vec;

            temp_vec = vec_ld(val_p, hxx);
            temp_vec = vec_madd(vec_pre20, sum0, temp_vec);
            vec_st(temp_vec, val_p, hxx);
            temp_vec = vec_ld(val_p, hxy);
            temp_vec = vec_madd(vec_pre11, sum0, temp_vec);
            vec_st(temp_vec, val_p, hxy);
            temp_vec = vec_ld(val_p, hxz);
            temp_vec = vec_madd(vec_pre10, sum1, temp_vec);
            vec_st(temp_vec, val_p, hxz);
            temp_vec = vec_ld(val_p, hyy);
            temp_vec = vec_madd(vec_pre02, sum0, temp_vec);
            vec_st(temp_vec, val_p, hyy);
            temp_vec = vec_ld(val_p, hyz);
            temp_vec = vec_madd(vec_pre01, sum1, temp_vec);
            vec_st(temp_vec, val_p, hyz);
            temp_vec = vec_ld(val_p, hzz);
            temp_vec = vec_madd(vec_pre00, sum2, temp_vec);
            vec_st(temp_vec, val_p, hzz);

            temp_vec = vec_ld(val_p, gx);
            temp_vec = vec_madd(vec_pre10, sum0, temp_vec);
            vec_st(temp_vec, val_p, gx);
            temp_vec = vec_ld(val_p, gy);
            temp_vec = vec_madd(vec_pre01, sum0, temp_vec);
            vec_st(temp_vec, val_p, gy);
            temp_vec = vec_ld(val_p, gz);
            temp_vec = vec_madd(vec_pre00, sum1, temp_vec);
            vec_st(temp_vec, val_p, gz);

            temp_vec = vec_ld(val_p, vals);
            temp_vec = vec_madd(vec_pre00, sum0, temp_vec);
            vec_st(temp_vec, val_p, vals);

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

      for (int n=0; n<num_splines; n++)
      {
        gx[n]*=dxInv; 
        gy[n]*=dyInv; 
        gz[n]*=dzInv; 
        hxx[n]*=dxx;
        hyy[n]*=dyy;
        hzz[n]*=dzz;
        hxy[n]*=dxy;
        hxz[n]*=dxz;
        hyz[n]*=dyz;
      }
    }
#endif
}/** qmcplusplus namespace */
#endif

