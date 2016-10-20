//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp. 
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_MULTIEINSPLINE_VALUE_HPP
#define QMCPLUSPLUS_MULTIEINSPLINE_VALUE_HPP
namespace qmcplusplus
{
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

      const intptr_t xs = spline_m->x_stride;
      const intptr_t ys = spline_m->y_stride;
      const intptr_t zs = spline_m->z_stride;

      const T zero(0);
      const int num_splines=last-first;
      ASSUME_ALIGNED(vals);
      std::fill(vals,vals+num_splines,zero);

      constexpr int iSIMD=QMC_CLINE/sizeof(T);
      const int nBlocks = num_splines/iSIMD;
      //unnecessary
      //assert ( nBlocks*iSIMD == num_splines);

      for (int i=0; i<4; i++){
        for (int j=0; j<4; j++){

          const T pre00 =  a[i]*b[j];
          
          const T* restrict coefs = spline_m->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs)+first; ASSUME_ALIGNED(coefs);
          const T* restrict coefszs  = coefs+zs;       ASSUME_ALIGNED(coefszs);
          const T* restrict coefs2zs = coefs+2*zs;     ASSUME_ALIGNED(coefs2zs);
          const T* restrict coefs3zs = coefs+3*zs;     ASSUME_ALIGNED(coefs3zs);

          for(int n=0; n<nBlocks; n++) {
            int nnBlocks = n*iSIMD;
#if defined(__AVX512F__) && defined(QMC_PREFETCH)
            {

                int pfi = (j==3) ? i+1 : i;
                int pfj = ((j+1)%4);

                const T* restrict coefs = spline_m->coefs + ((ix+pfi)*xs + (iy+pfj)*ys + iz*zs)+first; ASSUME_ALIGNED(coefs);
                const T* restrict coefszs  = coefs+zs;       ASSUME_ALIGNED(coefszs);
                const T* restrict coefs2zs = coefs+2*zs;     ASSUME_ALIGNED(coefs2zs);
                const T* restrict coefs3zs = coefs+3*zs;     ASSUME_ALIGNED(coefs3zs);

                _mm_prefetch((char const*)(coefs+nnBlocks),_MM_HINT_T1);
                _mm_prefetch((char const*)(coefszs+nnBlocks),_MM_HINT_T1);
                _mm_prefetch((char const*)(coefs2zs+nnBlocks),_MM_HINT_T1);
                _mm_prefetch((char const*)(coefs3zs+nnBlocks),_MM_HINT_T1);
            }
#endif

#pragma omp simd
            for ( int m = 0; m < iSIMD; m++ ) {
              int mIdx = nnBlocks + m;
              vals[mIdx] += pre00*(c[0]*coefs[mIdx] + c[1]*coefszs[mIdx] + c[2]*coefs2zs[mIdx] + c[3]*coefs3zs[mIdx]);
            }

          }


        } // j loop
      } // i loop 
      
    }
}
#endif
