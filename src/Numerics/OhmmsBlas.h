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
#ifndef OHMMSNUMERIC_BLAS_H
#define OHMMSNUMERIC_BLAS_H

//generic header for blas routines
#include "Numerics/Blasf.h"

/** Interfaces to blas library
 *
 *   static data members to facilitate C++/Fortran blas interface
 *   static member functions to use blas functions
 *
 *     inline static void axpy
 *     inline static double norm2
 *     inline static float norm2
 *     inline static void symv
 *     inline static double dot
 *     inline static float dot
 *
 *  Arguments (float/double/complex<float>/complex<double>) determine
 *  which BLAS routines are actually used.
 *  Note that symv can be call in many ways.
 */
struct BLAS {

  static const int INCX = 1;
  static const int INCY = 1;
  static const char UPLO = 'L';
  static const float sone;
  static const float szero;
  static const double done;
  static const double dzero;
  static const complex<float> cone;
  static const complex<float> czero;
  static const complex<double> zone;
  static const complex<double> zzero;

  inline static 
  void axpy(int n, double x, const double* a, double* b){
    daxpy(n, x, a, INCX, b, INCY);
  }

  inline static 
  void axpy(int n, double x, const double* a, int incx, double* b, int incy){
    daxpy(n, x, a, incx, b, incy);
  }

  inline static 
  void axpy(int n, const double* a, double* b){
    daxpy(n, done, a, INCX, b, INCY);
  }

  inline static 
  void axpy(int n, float x, const float* a, int incx, float* b, int incy){
    saxpy(n, x, a, incx, b, incy);
  }

  inline static 
  void axpy(int n, const complex<double> x, const complex<double>* a, int incx, 
            complex<double>* b, int incy){
    zaxpy(n, x, a, incx, b, incy);
  }

  inline static 
  double norm2(int n, const double* a, int incx =1) {
    return dnrm2(n, a, incx);
  }

  inline static 
  double norm2(int n, const complex<double>* a, int incx =1) {
    return dznrm2(n, a, incx);
  }

  inline static 
  float norm2(int n, const float* a, int incx =1) {
    return snrm2(n, a, incx);
  }

  inline static 
  void scal(int n, double alpha, double* x) {
    dscal(n,alpha,x,INCX);
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

  inline static
  double dot(int n, const double* a, const double* b) {
    return ddot(n,a,INCX,b,INCY);
  }

  inline static
  float dot(int n, const float* a, const float* b) {
    return sdot(n,a,INCX,b,INCY);
  }

  inline static
  void copy(int n, const double* a, double* b) {
    dcopy(n,a,INCX,b,INCY);
  }

  inline static
  void copy(int n, const double* a, int ia, double* b, int ib) {
    dcopy(n,a,ia,b,ib);
  }

/*
  inline static
  void copy(int n, double x, double* a) {
    dinit(n,x,a,INCX);
  }
*/

  inline static
  void copy(int n, const complex<double>* a, complex<double>* b) {
    zcopy(n,a,INCX,b,INCY);
  }

  inline static
  void copy(int n, const complex<double>* a, int ia, complex<double>* b, int ib) {
    zcopy(n,a,ia,b,ib);
  }
};
#endif // OHMMS_BLAS_H
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

