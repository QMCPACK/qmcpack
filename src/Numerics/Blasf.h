//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
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
#ifndef OHMMS_BLAS_FUNCTIONDEFS_H
#define OHMMS_BLAS_FUNCTIONDEFS_H

#include <complex>
using namespace std;

#ifdef ADD_
#define daxpy daxpy_
#define saxpy saxpy_
#define zaxpy zaxpy_
#define dnrm2 dnrm2_
#define snrm2 snrm2_
#define dznrm2 dznrm2_
#define dsymv dsymv_
#define ssymv ssymv_
#define csymv csymv_
#define zsymv zsymv_
#define ddot  ddot_
#define sdot  sdot_
#define dscal  dscal_
#define dcopy dcopy_
#define zcopy zcopy_
#define dsyrk  dsyrk_
#define dsymm  dsymm_
#define dgemm  dgemm_
#define dgemv  dgemv_
#define dsyr2k dsyr2k_
#define dgetrf dgetrf_
#define dgetri dgetri_
#define zgetrf zgetrf_
#define zgetri zgetri_
#endif 

// declaring Fortran interfaces
extern "C" {

  double ddot(const int& n, const double *dx, const int& incx, const double *dy, const int &incy);

  float sdot(const int& n, const float *dx, const int& incx,  const float *dy, const int &incy);

  void daxpy(const int& n, const double& da,  
	     const double *dx, const int& incx, double *dy, const int& incy);

  void saxpy(const int& n, const float& da,  
	     const float *dx, const int& incx, float *dy, const int& incy);

  void zaxpy(const int& n, const complex<double>&  da,  const complex<double> *dx,
             const int& incx, complex<double> *dy, const int& incy);


  double dnrm2(const int& n, const double *dx, const int& incx);
  float  snrm2(const int& n, const float *dx, const int& incx);
  double dznrm2(const int& n, const  complex<double> *dx, const int& incx);


  double dscal(const int& n, const double&, double* x, const int&);

  void  dsymv(const char& uplo, const int& n, const double& alpha, 
              const double& a, const int& lda, const double* x, const int& incx, 
              const double& beta, double *y, const int& incy);

  void  ssymv(const char& uplo, const int& n, const float& alpha, 
              const float& a, const int& lda, 
	      const float* x, const int& incx, 
              const float& beta, float *y, const int& incy);

  void  zsymv(const char& uplo, const int& n, const complex<double>& alpha,
              complex<double>* a, const int& lda, complex<double>* x, 
              const int& incx, const complex<double>& beta, 
              complex<double> *y, const int& incy);

  void  csymv(const char& uplo, const int& n, const complex<float>& alpha, 
              complex<float>* a, const int& lda, complex<float>* x, 
              const int& incx, const complex<float>& beta, 
              complex<float> *y, const int& incy);

  void zcopy(const int& n, const complex<double>* x, 
	     const int& incx, complex<double>* y, const int& incy);
  void dcopy(const int& n, const double*, const int& , double *, const int&);


  void dsyr2k(const char&, const char&, const int&, const int&,
	      const double&, const double*, const int&,
	      const double*, const int&, 
	      const double&, double*, const int&);

  void dsymm(const char&, const char&, const int&, const int&, const double&, 
	     const double* A, const int& lda,
	     const double* B, const int& ldb,
             const double& beta, double* C, const int& ldc );

  void dgemm(const char&, const char&,
	     const int&, const int&, const int&,
	     const double&, const double*, const int&, const double*, const int&,
	     const double&, double*, const int&);

  void dgemv(const char& trans, const int& nr, const int& nc, 
	     const double& alpha, const double* amat, const int& lda, 
             const double* bv, const int& incx,
	     const double& beta, double* cv, const int& incy);

  void dsyrk(const char&, const char&, const int&, const int&,
             const double&, const double*, const int&,
             const double&, double*,const int&);

  void dgetrf(const int& n, const int& m, double* a, const int& n0, 
	      int* piv, int& st);

  void zgetrf(const int& n, const int& m, complex<double>* a, const int& n0, 
	      int* piv, int& st);

  void dgetri(const int& n, double* a, const int& n0, 
	      int* piv, double* work, const int&, int& st);


  void zgetri(const int& n, complex<double>* a, const int& n0, 
	      int* piv, complex<double>* work, const int&, int& st);
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
