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
/** @file blas1.hpp
 */
#ifndef QMCPLUSPLUS_SIMD_BLAS1_HPP
#define QMCPLUSPLUS_SIMD_BLAS1_HPP

//#include <Numerics/OhmmsBlas.h>
//#include <OhmmsPETE/TinyVector.h>

namespace qmcplusplus {

  //going to replace simd
  namespace blas {

    /** copy function using memcpy 
     * @param source starting address of the source
     * @param target starting address of the target
     * @param n size of the data to copy
     */
    template<typename T1, typename T2>
      inline void copy(const T1* restrict source, T2* restrict target, size_t n)
      {
        std::copy(source,source+n,target);
      }

    /** copy complex to two real containers */
    template<typename T1, typename T2>
      inline void copy(const std::complex<T1>* restrict source, T2* restrict target_r, T2* restrict target_i, size_t n)
      {
        const T1* restrict iptr=reinterpret_cast<T1*>(source);
        ASSUMED_ALIGNED(source); ASSUMED_ALIGNED(iptr);
        for(int i=0,i2=0; i<n; ++i,i2+=2)
        {
          target_r[i] =static_cast<T2>(iptr[i2+0]);
          target_i[i] =static_cast<T2>(iptr[i2+1]);
        }
      }

    /** dot product
     */
    template<typename T, typename TSUM>
      inline TSUM dot(const T* restrict a, const T* restrict b, int n, TSUM res) 
      {
        ASSUMED_ALIGNED(a); ASSUMED_ALIGNED(b);
        #pragma omp simd reduction(+:res)
        for(int i=0; i<n; i++) res += a[i]*b[i];
        return res;
      }

    /** dot product of complex and real
     */
    template<typename T, typename TSUM>
      inline std::complex<TSUM> 
      dot(const std::complex<T>* restrict a, const T* restrict b, int n, const std::complex<TSUM>& res) 
      {
        TSUM res_r=res.real(),res_i=res.imag();
        const T* restrict iptr=reinterpret_cast<T*>(a); ASSUMED_ALIGNED(iptr); 
        ASSUMED_ALIGNED(b);
        for(int i=0,i2=0; i<n; i++,i2+=2) {res_r+=b[i]*iptr[i2]; res_i+=b[i]*iptr[i2+1];}
        return std::complex<T>(res_r,res_i);
      }

    template<typename T, typename TSUM>
      inline std::complex<TSUM> 
      dot(const T* restrict a, const std::complex<T>* restrict b, int n, const std::complex<TSUM>& res) 
      {
        return dot(b,a,n,res);
      }

    template<typename T>
      inline void axpy(T alpha, const T* restrict in, T* restrict out, int n)
      {
        #pragma omp simd aligned(in,out)
        for(int i=0; i<n; ++i)
          out[i]=alpha*in[i];
      }

    template<typename T>
      inline void scal(T alpha, T* restrict inout, int n)
      {
        #pragma omp simd aligned(inout)
        for(int i=0; i<n; ++i)
          inout[i]*=alpha;
      }

    template<typename T>
      inline void accumulate_phases(const int& n
          , const std::complex<T> * restrict x_, const std::complex<T> * restrict y_
          , T& rN, T& iN, T& riN)
      {
        const T* restrict x=reinterpret_cast<T*>(x_);
        const T* restrict y=reinterpret_cast<T*>(y_);
        for(int i=0,i2=0; i<n; ++i,i2+=2)
        {
          T tr=x[i2]*y[i2  ]-x[i2+1]*y[i2+1];
          T ti=x[i2]*y[i2+1]+x[i2+1]*y[i2  ];
          //T tr=x[i].real()*y[i].real()-x[i].imag()*y[i].imag();
          //T ti=x[i].real()*y[i].imag()+x[i].imag()*y[i].real();
          rN += tr*tr;
          iN += ti*ti;
          riN += tr*ti;
        }//
      }
  } //simd namepsace

}
#endif

