//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


/** @file inner_product.h
 *
 * Inline dot and gemv functions. 
 * On modern processors with good compilers, there is no need to use blas1, e.g., ddot, but 
 * strange things can happen on some machines.
 * SIMD specialization can be implemented here.
 */
#ifndef QMCPLUSPLUS_INNER_PRODUCT_HPP
#define QMCPLUSPLUS_INNER_PRODUCT_HPP

#include <Numerics/OhmmsBlas.h>
#include <OhmmsPETE/TinyVector.h>

namespace qmcplusplus {

  namespace simd {

    /** copy function using memcpy 
     * @param target starting address of the target
     * @param source starting address of the source
     * @param n size of the data to copy
     */
    template<typename T1, typename T2>
      inline void copy(T1* restrict target, const T2* restrict source, size_t n)
      {
        for(size_t i=0; i<n;  ++i) target[i]=static_cast<T1>(source[i]);
      }

    /** copy function using memcpy 
     * @param target starting address of the target
     * @param source starting address of the source
     * @param n size of the data to copy
     */
    template<typename T>
      inline void copy(T* restrict target, const T* restrict source, size_t n)
      {
        memcpy(target,source,sizeof(T)*n);
      }

    /** copy complex to two real containers */
    template<typename T1, typename T2>
      inline void copy(T1* restrict target_r, T1* restrict target_i, const std::complex<T2>* restrict source, size_t n)
      {
        for(int i=0; i<n; ++i)
        {
          *target_r++=static_cast<T1>(source[i].real());
          *target_i++=static_cast<T1>(source[i].imag());
        }
      }

    /** dot product
     * @param a starting address of an array of type T
     * @param b starting address of an array of type T
     * @param n size
     * @param res initial value, default=0.0
     * @return  \f$ res = \sum_i a[i] b[i]\f$
     *
     * same as inner_product(a,a+n,b,0.0)
     * The arguments of this inline function are different from BLAS::dot
     * This is more efficient than BLAS::dot due to the function overhead,
     * if a compiler knows how to inline.
     */
    template<typename T>
      inline T dot(const T* restrict a, const T* restrict b, int n, T res=T()) 
      {
        for(int i=0; i<n; i++) res += a[i]*b[i];
        return res;
      }

    /** dot product of mixed type
     * @param a starting address of complex
     * @param b starting address of real
     * @param n size
     * @param res initial value of the sum
     */
    template<typename T>
      inline std::complex<T> 
      dot(const std::complex<T>* restrict a, const T* restrict b, int n, const std::complex<T>& res=T()) 
      {
        T res_r=res.real(),res_i=res.imag();
        for(int i=0; i<n; i++) {res_r+=b[i]*a[i].real(); res_i+=b[i]*a[i].imag();}
        return std::complex<T>(res_r,res_i);
      }

    /** dot product of mixed type
     * @param a starting address of real
     * @param b starting address of complex
     * @param n size
     * @param res initial value of the sum
     */
    template<typename T>
      inline std::complex<T> 
      dot(const T* restrict a, const std::complex<T>* restrict b, int n,  const std::complex<T>& res=T())
      {
        T res_r=res.real(),res_i=res.imag();
        for(int i=0; i<n; i++) {res_r+=a[i]*b[i].real(); res_i+=a[i]*b[i].imag();}
        return std::complex<T>(res_r,res_i);
      }

    /** inline dot product 
     * @param a starting address of an array of type T
     * @param b starting address of an array of type TinyVector<T,D>
     * @param n size
     * @return \f$ {\bf v} = \sum_i a[i] {\bf b}[i]\f$
     */
    template<class T, unsigned D>
      inline Tensor<T,D>
      dot(const T* a, const Tensor<T,D>* b, int n)
      {
        Tensor<T,D> res;
//#if defined(USE_GEMV_FOR_G_DOT_V)
//        BLAS::gemv(n,3,b->data(),a,res.data());
//#else
        for(int i=0; i<n; i++) res += a[i]*b[i];
//#endif
        return res;
      }
      
     template<class T, unsigned D>
      inline TinyVector<T,D>
      dot(const T* a, const TinyVector<T,D>* b, int n)
      {
        TinyVector<T,D> res;
#if defined(USE_GEMV_FOR_G_DOT_V)
        BLAS::gemv(n,3,b->data(),a,res.data());
#else
        for(int i=0; i<n; i++) res += a[i]*b[i];
#endif
        return res;
      }

    /** x*y dot product of two vectors using the same arugment list for blas::dot
     * @param n size
     * @param x starting address of x 
     * @param incx stride of x
     * @param y starting address of y 
     * @param incx stride of y
     * @param return \f$\sum_i x[i+=incx]*y[i+=incy]\f$
     */
    template<typename T>
      inline 
      T dot(int n, const T* restrict x, int incx, const T* restrict y, int incy) 
      {
        const int xmax=incx*n;
        T res=0.0;
        for(int ic=0,jc=0; ic<xmax; ic+=incx,jc+=incy) res += x[ic]*y[jc];
        return res;
      }

    template<typename T>
      inline void accumulate_phases(const int& n
          , const std::complex<T> * restrict x, const std::complex<T> * restrict y
          , T& rN, T& iN, T& riN)
      {
        for(int i=0; i<n; ++i)
        {
          T tr=x[i].real()*y[i].real()-x[i].imag()*y[i].imag();
          T ti=x[i].real()*y[i].imag()+x[i].imag()*y[i].real();
          rN += tr*tr;
          iN += ti*ti;
          riN += tr*ti;
        }//
      }

    /** transpose of A(m,n) to B(n,m) 
     * @param A starting address, A(m,lda)
     * @param m number of A rows
     * @param lda stride of A's row
     * @param B starting address B(n,ldb)
     * @param n number of B rows
     * @param ldb stride of B's row
     *
     * Blas-like interface
     */
    template<typename T, typename TO>
      inline void transpose(const T* restrict A, size_t m, size_t lda,  TO* restrict B, size_t n, size_t ldb)
      {
        for(size_t i=0; i < n; ++i)
          for(size_t j=0; j < m; ++j)
            B[i*ldb+j]=A[j*lda+i];
      }

  } //simd namepsace
}
#endif
