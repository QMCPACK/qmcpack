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


/**@file dirac_computeGL.h
 * Utility functions to compute gradient and lpalacian using VectorSoaContainer
 */
#ifndef QMCPLUSPLUS_DET_AUX_H
#define QMCPLUSPLUS_DET_AUX_H

//use blas::gemv to compute GL(D,N)*psiV(N) 
#define QMC_USE_GEMV_FOR_GL 1

namespace qmcplusplus
{

#if QMC_USE_GEMV_FOR_GL
  template<typename T, typename T2>
  inline void computeGL(T* row, VectorSoaContainer<T,4>& gl_v, TinyVector<T2,3>& grad, T2& lap)
  {
    CONSTEXPR T czero(0);
    CONSTEXPR T cone(1);
    int four=4;
    int na=gl_v.size();
    int lda=gl_v.capacity();
    T y[]={czero,czero,czero,czero};
    BLAS::gemv('T',na,four,cone,gl_v.data(),lda,row,1,czero,y,1);
    grad[0]=y[0]; grad[1]=y[1]; grad[2]=y[2]; lap=y[3];
  }

  template<typename T>
  inline TinyVector<T,3> computeG(T* row, VectorSoaContainer<T,4>& gl_v)
  {
    constexpr T czero(0);
    constexpr T cone(1);
    int three=3;
    int na=gl_v.size();
    int lda=gl_v.capacity();
    T y[]={czero,czero,czero,czero};
    BLAS::gemv('T',na,three,cone,gl_v.data(),lda,row,1,czero,y,1);
    return TinyVector<T,3>(y[0],y[1],y[2]);
  }

  template<typename T>
  inline void computeGL(T* row, VectorSoaContainer<T,5>& gl_v, TinyVector<T,3>& grad, T& lap)
  {
    CONSTEXPR T czero(0);
    CONSTEXPR T cone(1);
    int four=4;
    int na=gl_v.size();
    int lda=gl_v.capacity();
    T y[]={czero,czero,czero,czero};
    BLAS::gemv('T',na,four,cone,gl_v.data(1),lda,row,1,czero,y,1);
    grad[0]=y[0]; grad[1]=y[1]; grad[2]=y[2]; lap=y[3];
  }

  template<typename T>
  inline TinyVector<T,3> computeG(T* row, VectorSoaContainer<T,5>& gl_v)
  {
    constexpr T czero(0);
    constexpr T cone(1);
    int three=3;
    int na=gl_v.size();
    int lda=gl_v.capacity();
    T y[]={czero,czero,czero,czero};
    BLAS::gemv('T',na,three,cone,gl_v.data(1),lda,row,1,czero,y,1);
    return TinyVector<T,3>(y[0],y[1],y[2]);
  }

#else

#if QMC_COMPLEX
#error "Cannot do complex yet with compute GL. Use GEMV\n"
#else
  ///real version using simd: only for testing
  template<typename T, typename T2>
  inline void computeGL(T* row, VectorSoaContainer<T,4>& gl_v, TinyVector<T2,3>& grad, T2& lap)
  {
    constexpr T czero(0);
    constexpr T cone(1);
    const T* restrict gx_p=gl_v.data(0);
    const T* restrict gy_p=gl_v.data(1);
    const T* restrict gz_p=gl_v.data(2);
    const T* restrict l_p=gl_v.data(3);
    T gx=czero, gy=czero, gz=czero,l=czero;
    const int n=gl_v.size();
#pragma omp simd reduction(+:gx,gy,gz,l)
    for(size_t i=0; i<n; ++i)
    {
      gx +=row[i]*gx_p[i];
      gy +=row[i]*gy_p[i];
      gz +=row[i]*gz_p[i];
      l+=row[i]*l_p[i];
    }
    grad[0]=gx;
    grad[1]=gy;
    grad[2]=gz;
    lap=l;
  }

  template<typename T>
  inline TinyVector<T,3> computeG(T* row, VectorSoaContainer<T,4>& gl_v)
  {
    constexpr T czero(0);
    constexpr T cone(1);
    const T* restrict gx_p=gl_v.data(0);
    const T* restrict gy_p=gl_v.data(1);
    const T* restrict gz_p=gl_v.data(2);
    T gx=czero, gy=czero, gz=czero;
    const int n=gl_v.size();
#pragma omp simd reduction(+:gx,gy,gz)
    for(size_t i=0; i<n; ++i)
    {
      gx+=row[i]*gx_p[i];
      gy+=row[i]*gy_p[i];
      gz+=row[i]*gz_p[i];
    }
    return TinyVector<T,3>(gx,gy,gz);
  }

  ///real version using simd: only for testing
  template<typename T>
  inline void computeGL(T* row, VectorSoaContainer<T,5>& gl_v, TinyVector<T,3>& grad, T& lap)
  {
    constexpr T czero(0);
    constexpr T cone(1);
    const T* restrict gx_p=gl_v.data(1);
    const T* restrict gy_p=gl_v.data(2);
    const T* restrict gz_p=gl_v.data(3);
    const T* restrict l_p=gl_v.data(4);
    lap=czero;
    T gx=czero, gy=czero, gz=czero,l=czero;
    const int n=gl_v.size();
#pragma omp simd reduction(+:gx,gy,gz,l)
    for(size_t i=0; i<n; ++i)
    {
      gx +=row[i]*gx_p[i];
      gy +=row[i]*gy_p[i];
      gz +=row[i]*gz_p[i];
      l+=row[i]*l_p[i];
    }
    grad[0]=gx;
    grad[1]=gy;
    grad[2]=gz;
    lap=l;
  }

  template<typename T>
  inline TinyVector<T,3> computeG(T* row, VectorSoaContainer<T,5>& gl_v)
  {
    constexpr T czero(0);
    constexpr T cone(1);
    const T* restrict gx_p=gl_v.data(1);
    const T* restrict gy_p=gl_v.data(2);
    const T* restrict gz_p=gl_v.data(3);
    T gx=czero, gy=czero, gz=czero;
    const int n=gl_v.size();
#pragma omp simd reduction(+:gx,gy,gz)
    for(size_t i=0; i<n; ++i)
    {
      gx+=row[i]*gx_p[i];
      gy+=row[i]*gy_p[i];
      gz+=row[i]*gz_p[i];
    }
    return TinyVector<T,3>(gx,gy,gz);
  }
#endif //QMC_COMPLEX

#endif

}
#endif
