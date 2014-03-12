//////////////////////////////////////////////////////////////////
// (c) Copyright 2011- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
      inline void copy(T1* restrict target_r, T1* restrict target_i, const complex<T2>* restrict source, size_t n)
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
      inline T dot(const T* restrict a, const T* restrict b, int n, T res=0.0) 
      {
#if defined(USE_BLAS_DOT)
        return BLAS::dot(n,a,1,b,1);
#else
        for(int i=0; i<n; i++) res += a[i]*b[i];
        return res;
#endif
      }

    /** dot product specialized for double */
    inline double dot(const double* restrict a, const double* restrict b, int n, double res=0.0) 
    {
#if defined(USE_BLAS_DOT)
      return ddot(n,a,1,b,1);
#else
      for(int i=0;i<n; ++i) res += a[i]*b[i];
      return res;
#endif
    }

    /** dot product of mixed type
     * @param a starting address of complex
     * @param b starting address of real
     * @param n size
     * @param res initial value of the sum
     */
    template<typename T>
      inline std::complex<T> 
      dot(const std::complex<T>* restrict a, const T* restrict b, int n, const std::complex<T>& res=0.0) 
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
      dot(const T* restrict a, const std::complex<T>* restrict b, int n,  const std::complex<T>& res=0.0)
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

    //template<typename T>
    //  inline void gemv(const Matrix<T>& a, const T* restrict v, T* restrict b)
    //  {
    //    MatrixOperators::product(a,v,b);
    //  }

    //template<typename T,unsigned D>
    //  inline void gemv(const Matrix<T>& a, const TinyVector<T,D>* restrict v, TinyVector<T,D>* restrict b)
    //  {
    //    MatrixOperators::product(a,v,b);
    //  }

    //template<typename T,unsigned D>
    //  inline void gemv(const Matrix<T>& a, const Tensor<T,D>* restrict v, Tensor<T,D>* restrict b)
    //  {
    //    MatrixOperators::product(a,v,b);
    //  }


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

    /** transpose 2D data
     * @param in starting address of n*m
     * @param out starting address of the transposed data
     * @param n rows of in
     * @param m cols of in
     *
     * in(n,m) -> out(m,n)
     */
    template<typename T>
      inline void transpose(const T* restrict in, T* restrict out, int n, int m)
      {
        for(int i=0; i<m; ++i)
          for(int j=0,jj=i; j<n; ++j,jj+=m) *out++ = in[jj];
      }

  } //simd namepsace

//    template<typename T>
//    inline T dot(const TinyVector<T,2>& a, const TinyVector<T,2>& b)
//    {
//      return a[0]*b[0]+a[1]*b[1];
//    }
//
//    template<typename T>
//    inline std::complex<T> dot(const TinyVector<complex<T>,2>& a, const TinyVector<complex<T>,2>& b) 
//    { 
//      return a[0]*b[0]+a[1]*b[1];
//    }
//
//    template<typename T>
//    inline T dot(const TinyVector<T,3>& a, const TinyVector<T,3>& b)
//    {
//      return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
//    }
//
//    template<typename T>
//    inline std::complex<T> dot(const TinyVector<complex<T>,3>& a, const TinyVector<complex<T>,3>& b) 
//    { 
//      return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; 
//    }
//
//    template<typename T>
//    inline T dot(const TinyVector<T,4>& a, const TinyVector<T,4>& b)
//    {
//      return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3]; 
//    }
//
//    template<typename T>
//    inline std::complex<T> dot(const TinyVector<complex<T>,4>& a, const TinyVector<complex<T>,4>& b) 
//    { 
//      return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3]; 
//    }
//
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5077 $   $Date: 2010-12-09 03:14:51 -0600 (Thu, 09 Dec 2010) $
 * $Id: inner_product.hpp 5077 2010-12-09 09:14:51Z jmcminis $ 
 ***************************************************************************/
