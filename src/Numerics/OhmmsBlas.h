//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_NUMERIC_BLAS_H
#define QMCPLUSPLUS_NUMERIC_BLAS_H

//generic header for blas routines
#include "Numerics/Blasf.h"

/** Interfaces to blas library
 *
 *   static data members to facilitate C++/Fortran blas interface
 *   static member functions to use blas functions
 *   - inline static void axpy
 *   - inline static double norm2
 *   - inline static float norm2
 *   - inline static void symv
 *   - inline static double dot
 *   - inline static float dot
 *
 *  Arguments (float/double/complex\<float\>/complex\<double\>) determine
 *  which BLAS routines are actually used.
 *  Note that symv can be call in many ways.
 */
struct BLAS
{

  static const int INCX = 1;
  static const int INCY = 1;
  static const char UPLO = 'L';
  static const char TRANS = 'T';
  static const char NOTRANS = 'N';
  static const float sone;
  static const float szero;
  static const double done;
  static const double dzero;
  static const complex<float> cone;
  static const complex<float> czero;
  static const complex<double> zone;
  static const complex<double> zzero;

  inline static
  void axpy(int n, double x, const double* a, double* b)
  {
    daxpy(n, x, a, INCX, b, INCY);
  }

  inline static
  void axpy(int n, double x, const double* a, int incx, double* b, int incy)
  {
    daxpy(n, x, a, incx, b, incy);
  }

  inline static
  void axpy(int n, const double* a, double* b)
  {
    daxpy(n, done, a, INCX, b, INCY);
  }

  inline static
  void axpy(int n, float x, const float* a, int incx, float* b, int incy)
  {
    saxpy(n, x, a, incx, b, incy);
  }

  inline static
  void axpy(int n, const complex<double> x, const complex<double>* a, int incx,
            complex<double>* b, int incy)
  {
    zaxpy(n, x, a, incx, b, incy);
  }

  inline static
  double norm2(int n, const double* a, int incx =1)
  {
    return dnrm2(n, a, incx);
  }

  inline static
  double norm2(int n, const complex<double>* a, int incx =1)
  {
    return dznrm2(n, a, incx);
  }

  inline static
  float norm2(int n, const float* a, int incx =1)
  {
    return snrm2(n, a, incx);
  }

  inline static
  void scal(int n, double alpha, double* x)
  {
    dscal(n,alpha,x,INCX);
  }

  //inline static
  //void gemv(char trans, int n, int m, const double* amat, const double* x, double* y) {
  //  dgemv(trans, n, m, done, amat, n, x, INCX, dzero, y, INCY);
  //}
  inline static
  void gemv(int n, int m, const double* restrict amat, const double* restrict x, double* restrict y)
  {
    dgemv(NOTRANS, m, n, done, amat, m, x, INCX, dzero, y, INCY);
  }

  inline static
  void gemv_trans(int n, int m, const double* restrict amat, const double* restrict x, double* restrict y)
  {
    dgemv(TRANS, m, n, done, amat, m, x, INCX, dzero, y, INCY);
  }

  inline static
  void gemv_trans(int n, int m, const std::complex<double>* restrict amat, const std::complex<double>* restrict x
                  , std::complex<double>* restrict y)
  {
    zgemv(TRANS, m, n, done, amat, m, x, INCX, dzero, y, INCY);
  }

  inline static
  void gemv(int n, int m,
            const std::complex<double>* restrict amat,
            const std::complex<double>* restrict x,
            std::complex<double>* restrict y)
  {
    zgemv(NOTRANS, m, n, zone, amat, m, x, INCX, zzero, y, INCY);
  }

  inline static
  void gemv(char trans_in, int n, int m
            ,double alpha, const double* restrict amat, int lda
            , const double* x, int incx, double beta
            , double* y, int incy)
  {
    dgemv(trans_in, n, m, alpha, amat, lda, x, incx, beta, y, incy);
  }

  inline static
  void gemv(char trans_in, int n, int m
            , const complex<double>& alpha, const complex<double>* restrict amat, int lda
            , const complex<double>* restrict x, int incx, const complex<double>& beta
            , complex<double>* y, int incy)
  {
    zgemv(trans_in, n, m, alpha, amat, lda, x, incx, beta, y, incy);
  }

  inline static
  void gemm (char Atrans, char Btrans, int M, int N, int K, double alpha,
             const double *A, int lda, const double* restrict B, int ldb,
             double beta, double* restrict C, int ldc)
  {
    dgemm (Atrans, Btrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  }

  inline static
  void gemm (char Atrans, char Btrans, int M, int N, int K, complex<double> alpha,
             const complex<double> *A, int lda, const complex<double>* restrict B, int ldb,
             complex<double> beta, complex<double>* restrict C, int ldc)
  {
    zgemm (Atrans, Btrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  }


//   inline static
//   void symv(char uplo, int n, const double alpha, double* a, int lda,
//             double* x, const int incx, const double beta, double* y,
//             const int incy) {
//     dsymv(&uplo,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy);
//   }

//   inline static
//   void symv(char uplo, int n, const complex<double> alpha,
//             complex<double>* a, int lda, complex<double>* x, const int incx,
//             const complex<double> beta, complex<double>* y, const int incy) {
//     zsymv(&uplo,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy);
//   }

//   inline static
//   void symv(const char uplo, int n, const float alpha, float* a, int lda,
//             float* x, const int incx, const float beta, float* y,
//             const int incy) {
//     ssymv(&uplo,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy);
//   }

//   inline static
//   void symv(const char uplo, int n, const complex<float> alpha,
//             complex<float>* a, int lda, complex<float>* x, const int incx,
//             const complex<float> beta, complex<float>* y, const int incy) {
//     csymv(&uplo,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy);
//   }

//   inline static
//   void symv(int n, double alpha, double* a, double* x, double* y) {
//     dsymv(&UPLO,&n,&alpha,a,&n,x,&INCX,&dzero,y,&INCY);
//   }

//   inline static
//   void symv(int n, const double* a, const double* x, double* y) {
//     dsymv(&UPLO,&n,&done,a,&n,x,&INCX,&dzero,y,&INCY);
//   }

//   inline static
//   void symv(int n, float alpha, float* a, float* x, float* y) {
//     ssymv(&UPLO,&n,&alpha,a,&n,x,&INCX,&szero,y,&INCY);
//   }

//   inline static
//   void symv(int n, float* a, float* x, float* y) {
//     ssymv(&UPLO,&n,&sone,a,&n,x,&INCX,&szero,y,&INCY);
//   }

//   inline static void
//   symv(int n, complex<double> alpha, complex<double>* a, complex<double>* x,
//        complex<double>* y) {
//     zsymv(&UPLO,&n,&alpha,a,&n,x,&INCX,&zzero,y,&INCY);
//   }

//   inline static
//   void symv(int n, complex<double>* a, complex<double>* x, complex<double>* y) {
//     zsymv(&UPLO,&n,&zone,a,&n,x,&INCX,&zzero,y,&INCY);
//   }

//   inline static void
//   symv(int n, complex<float> alpha, complex<float>* a, complex<float>* x,
//        complex<float>* y) {
//     csymv(&UPLO,&n,&alpha,a,&n,x,&INCX,&czero,y,&INCY);
//   }

//   inline static
//   void symv(int n, complex<float>* a, complex<float>* x, complex<float>* y) {
//     csymv(&UPLO,&n,&cone,a,&n,x,&INCX,&czero,y,&INCY);
//   }

  //template<typename T>
  inline static
  double dot(int n, const double* restrict a, const int incx, const double* restrict b, const int incy )
  {
    return ddot(n,a,incx,b,incy);
  }

  inline static
  complex<double> dot(int n, const complex<double>* restrict a, const int incx, const complex<double>* restrict b, const int incy )
  {
    return zdotu(n,a,incx,b,incy);
  }

  template<typename T>
  inline static
  T dot(int n, const T* restrict a, const T* restrict b)
  {
    T res=0.0;
    for(int i=0; i<n; ++i)
      res += a[i]*b[i];
    return res;
  }

  template<typename T>
  inline static
  complex<T> dot(int n, const complex<T>* restrict a, const T* restrict b)
  {
    complex<T> res=0.0;
    for(int i=0; i<n; ++i)
      res += a[i]*b[i];
    return res;
  }

  template<typename T>
  inline static
  complex<T> dot(int n, const complex<T>* restrict a, const complex<T>* restrict b)
  {
    complex<T> res=0.0;
    for(int i=0; i<n; ++i)
      res += a[i]*b[i];
    return res;
  }


  template<typename T>
  inline static
  complex<T> dot(int n, const T* restrict a, const complex<T>* restrict b)
  {
    complex<T> res=0.0;
    for(int i=0; i<n; ++i)
      res += a[i]*b[i];
    return res;
  }

  template<typename T>
  inline static
  void copy(int n, const T* restrict a, T* restrict b)
  {
    memcpy(b,a,sizeof(T)*n);
  }

  /** copy using memcpy(target,source,size)
   * @param target starting address of the targe
   * @param source starting address of the source
   * @param number of elements to copy
   */
  template<typename T>
  inline static
  void copy(T* restrict target, const T* restrict source, int n)
  {
    memcpy(target,source,sizeof(T)*n);
  }

  template<typename T>
  inline static
  void copy(int n, const complex<T>* restrict a, T* restrict b)
  {
    for(int i=0; i<n; ++i)
      b[i]=a[i].real();
  }

  template<typename T>
  inline static
  void copy(int n, const T* restrict a, complex<T>* restrict b)
  {
    for(int i=0; i<n; ++i)
      b[i]=a[i];
  }

  template<typename T>
  inline static
  void copy(int n, const T* restrict x, int incx, T* restrict y, int incy)
  {
    const int xmax=incx*n;
    for(int ic=0,jc=0; ic<xmax; ic+=incx,jc+=incy)
      y[jc]=x[ic];
  }

  /*
    inline static
    void copy(int n, double x, double* a) {
      dinit(n,x,a,INCX);
    }

    inline static
    void copy(int n, const complex<double>* restrict a, complex<double>* restrict b)
    {
      zcopy(n,a,INCX,b,INCY);
    }

    inline static
    void copy(int n, const complex<double>* restrict a, int ia, complex<double>* restrict b, int ib) {
      zcopy(n,a,ia,b,ib);
    }
  */

  inline static
  void ger(int m, int n, double alpha
           , const double* x, int incx
           , const double* y , int incy
           , double* a, int lda)
  {
    dger(&m,&n,&alpha,x,&incx,y,&incy,a,&lda);
  }

  inline static
  void ger(int m, int n, const complex<double>& alpha
           , const complex<double>* x
           , int incx, const complex<double>* y, int incy
           , complex<double>* a, int lda)
  {
    zgeru(&m,&n,&alpha,x,&incx,y,&incy,a,&lda);
  }

};
#endif // OHMMS_BLAS_H
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

