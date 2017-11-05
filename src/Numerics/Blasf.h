//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef OHMMS_BLAS_FUNCTIONDEFS_H
#define OHMMS_BLAS_FUNCTIONDEFS_H

#include <complex>
#include <cstring>

#ifdef ADD_
#define caxpy caxpy_
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
#define dscal dscal_
#define sscal sscal_
#define dcopy dcopy_
#define zcopy zcopy_
#define dsyrk  dsyrk_
#define dsymm  dsymm_
#define dgemm  dgemm_
#define sgemm  sgemm_
#define zgemm  zgemm_
#define cgemm  cgemm_
#define dgemv  dgemv_
#define sgemv  sgemv_
#define zgemv  zgemv_
#define cgemv  cgemv_
#define dsyr2k dsyr2k_
#define dgetrf dgetrf_
#define dgetri dgetri_
#define sgetrf sgetrf_
#define sgetri sgetri_
#define zgetrf zgetrf_
#define zgetri zgetri_
#define cgetrf cgetrf_
#define cgetri cgetri_
#define dgesvd dgesvd_
#define sgesvd sgesvd_
#define dgeev dgeev_
#define sgeev sgeev_
#define dggev dggev_
#define sggev sggev_
#define dger dger_
#define sger sger_
#define zgeru zgeru_
#define cgeru cgeru_

#define dgeqrf dgeqrf_
#define dormqr dormqr_
#define dgghrd dgghrd_
#define dhgeqz dhgeqz_
#define dtgexc dtgexc_
#define dtgevc dtgevc_

#define dsyevr dsyevr_
#define zheevr zheevr_
#define zhegvx zhegvx_
#define zgeqrf zgeqrf_
#define zungqr zungqr_

#define cgeqrf cgeqrf_
#define sgeqrf sgeqrf_
#define dorgqr dorgqr_
#define sorgqr sorgqr_
#define cungqr cungqr_
#define zgelqf zgelqf_
#define dgelqf dgelqf_
#define cgelqf cgelqf_
#define sgelqf sgelqf_
#define dorglq dorglq_
#define sorglq sorglq_
#define zunglq zunglq_
#define cunglq cunglq_

#if defined(HAVE_MKL)
#define dzgemv  dzgemv_
#define scgemv  scgemv_
#define dzgemm  dzgemm_
#define scgemm  scgemm_
#endif

#endif

// declaring Fortran interfaces
extern "C" {

  void caxpy(const int& n, const std::complex<float>&  da,  const std::complex<float> *dx,
             const int& incx, std::complex<float> *dy, const int& incy);

  void daxpy(const int& n, const double& da,
             const double *dx, const int& incx, double *dy, const int& incy);

  void saxpy(const int& n, const float& da,
             const float *dx, const int& incx, float *dy, const int& incy);

  void zaxpy(const int& n, const std::complex<double>&  da,  const std::complex<double> *dx,
             const int& incx, std::complex<double> *dy, const int& incy);


  double dnrm2(const int& n, const double *dx, const int& incx);
  float  snrm2(const int& n, const float *dx, const int& incx);
  double dznrm2(const int& n, const  std::complex<double> *dx, const int& incx);


  double dscal(const int& n, const double&, double* x, const int&);
  double sscal(const int& n, const float&, float* x, const int&);

  void  dsymv(const char& uplo, const int& n, const double& alpha,
              const double& a, const int& lda, const double* x, const int& incx,
              const double& beta, double *y, const int& incy);

  void  ssymv(const char& uplo, const int& n, const float& alpha,
              const float& a, const int& lda,
              const float* x, const int& incx,
              const float& beta, float *y, const int& incy);

  void  zsymv(const char& uplo, const int& n, const std::complex<double>& alpha,
              std::complex<double>* a, const int& lda, std::complex<double>* x,
              const int& incx, const std::complex<double>& beta,
              std::complex<double> *y, const int& incy);

  void  csymv(const char& uplo, const int& n, const std::complex<float>& alpha,
              std::complex<float>* a, const int& lda, std::complex<float>* x,
              const int& incx, const std::complex<float>& beta,
              std::complex<float> *y, const int& incy);

  void zcopy(const int& n, const std::complex<double>* x,
             const int& incx, std::complex<double>* y, const int& incy);
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

  void sgemm(const char&, const char&,
             const int&, const int&, const int&,
             const float&, const float*, const int&, const float*, const int&,
             const float&, float*, const int&);

  void zgemm(const char&, const char&,
             const int&, const int&, const int&,
             const std::complex<double>&, const std::complex<double>*, const int&, const std::complex<double>*, const int&,
             const std::complex<double>&, std::complex<double>*, const int&);

  void cgemm(const char&, const char&,
             const int&, const int&, const int&,
             const std::complex<float>&, const std::complex<float>*, const int&, const std::complex<float>*, const int&,
             const std::complex<float>&, std::complex<float>*, const int&);

  void dgemv(const char& trans, const int& nr, const int& nc,
             const double& alpha, const double* amat, const int& lda,
             const double* bv, const int& incx,
             const double& beta, double* cv, const int& incy);

  void sgemv(const char& trans, const int& nr, const int& nc,
             const float& alpha, const float* amat, const int& lda,
             const float* bv, const int& incx,
             const float& beta, float* cv, const int& incy);

  void zgemv(const char& trans, const int& nr, const int& nc,
             const std::complex<double>& alpha, const std::complex<double>* amat, const int& lda,
             const std::complex<double>* bv, const int& incx,
             const std::complex<double>& beta, std::complex<double>* cv, const int& incy);

  void cgemv(const char& trans, const int& nr, const int& nc,
             const std::complex<float>& alpha, const std::complex<float>* amat, const int& lda,
             const std::complex<float>* bv, const int& incx,
             const std::complex<float>& beta, std::complex<float>* cv, const int& incy);

#if defined(HAVE_MKL)

  void dzgemm(const char&, const char&,
             const int&, const int&, const int&,
             const std::complex<double>&, const double*, const int&, const std::complex<double>*, const int&,
             const std::complex<double>&, std::complex<double>*, const int&);

  void scgemm(const char&, const char&,
             const int&, const int&, const int&,
             const std::complex<float>&, const float*, const int&, const std::complex<float>*, const int&,
             const std::complex<float>&, std::complex<float>*, const int&);

  void dzgemv(const char& trans, const int& nr, const int& nc,
             const std::complex<double>& alpha, const double* amat, const int& lda,
             const std::complex<double>* bv, const int& incx,
             const std::complex<double>& beta, std::complex<double>* cv, const int& incy);

  void scgemv(const char& trans, const int& nr, const int& nc,
             const std::complex<float>& alpha, const float* amat, const int& lda,
             const std::complex<float>* bv, const int& incx,
             const std::complex<float>& beta, std::complex<float>* cv, const int& incy);

#endif

  void dsyrk(const char&, const char&, const int&, const int&,
             const double&, const double*, const int&,
             const double&, double*,const int&);

  void dgetrf(const int& n, const int& m, double* a, const int& n0,
              int* piv, int& st);

  void sgetrf(const int& n, const int& m, float* a, const int& n0,
              int* piv, int& st);

  void zgetrf(const int& n, const int& m, std::complex<double>* a, const int& n0,
              int* piv, int& st);

  void cgetrf(const int& n, const int& m, std::complex<float>* a, const int& n0,
              int* piv, int& st);

  void dgetri(const int& n, double* a, const int& n0,
              int* piv, double* work, const int&, int& st);

  void sgetri(const int& n, float* a, const int& n0,
              int* piv, float* work, const int&, int& st);

  void zgetri(const int& n, std::complex<double>* a, const int& n0,
              int* piv, std::complex<double>* work, const int&, int& st);

  void cgetri(const int& n, std::complex<float>* a, const int& n0,
              int* piv, std::complex<float>* work, const int&, int& st);

  void dgesvd(char *JOBU, char* JOBVT, int *M, int *N,
              double *A, int *LDA, double *S, double *U,
              int *LDU, double *VT, int *LDVT, double *work,
              int *LWORK, int *INFO);

  void sgesvd(char *JOBU, char* JOBVT, int *M, int *N,
              float *A, int *LDA, float *S, float *U,
              int *LDU, float *VT, int *LDVT, float *work,
              int *LWORK, int *INFO);

  void dgeev(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, double *ALPHAR, double *ALPHAI,
             double *VL, int *LDVL, double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO );

  void sgeev(char *JOBVL, char *JOBVR, int *N, float *A, int *LDA, float *ALPHAR, float *ALPHAI,
             float *VL, int *LDVL, float *VR, int *LDVR, float *WORK, int *LWORK, int *INFO );

  void dggev(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, double *B, int *LDB,double *ALPHAR, double *ALPHAI,
             double *BETA, double *VL, int *LDVL, double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO );

  void sggev(char *JOBVL, char *JOBVR, int *N, float *A, int *LDA, float *B, int *LDB,float *ALPHAR, float *ALPHAI,
             float *BETA, float *VL, int *LDVL, float *VR, int *LDVR, float *WORK, int *LWORK, int *INFO );

  void dsyevr (char &JOBZ, char &RANGE, char &UPLO, int &N, double *A, int &LDA, double &VL, double &VU, int &IL, 
             int &IU, double &ABSTOL, int &M, double *W, double* Z, int &LDZ, int* ISUPPZ, double *WORK, 
             int &LWORK, int* IWORK, int &LIWORK, int &INFO);

  void zheevr (char &JOBZ, char &RANGE, char &UPLO, int &N, std::complex<double> *A, int &LDA, double &VL, double &VU, 
             int &IL, int &IU, double &ABSTOL, int &M, double *W, std::complex<double>* Z, int &LDZ, int* ISUPPZ, 
             std::complex<double> *WORK, int &LWORK, double* RWORK, int &LRWORK, int* IWORK, int &LIWORK, int &INFO);

  void zhegvx (int&, char &JOBZ, char &RANGE, char &UPLO, int &N, std::complex<double> *A, int &LDA, std::complex<double> *B, 
             int &LDB, double &VL, double &VU, int &IL, int &IU, double &ABSTOL, int &M, double *W, std::complex<double>* Z, 
             int &LDZ, std::complex<double> *WORK, int &LWORK, double* RWORK, int* IWORK, int* IFAIL, int &INFO);

  void zgeqrf( const int &M, const int &N, std::complex<double> *A, const int &LDA, std::complex<double> *TAU, std::complex<double> *WORK, const int &LWORK, int &INFO );

  void cgeqrf( const int &M, const int &N, std::complex<float> *A, const int &LDA, std::complex<float> *TAU, std::complex<float> *WORK, const int &LWORK, int &INFO );

  void dgeqrf( const int &M, const int &N, double *A, const int &LDA, double *TAU, double *WORK, const int &LWORK, int &INFO );

  void sgeqrf( const int &M, const int &N, float *A, const int &LDA, float *TAU, float *WORK, const int &LWORK, int &INFO );


  void zungqr( const int &M, const int &N, const int &K, std::complex<double> *A, const int &LDA, std::complex<double> *TAU, std::complex<double> *WORK, const int &LWORK, int &INFO );

  void cungqr( const int &M, const int &N, const int &K, std::complex<float> *A, const int &LDA, std::complex<float> *TAU, std::complex<float> *WORK, const int &LWORK, int &INFO );

  void dorgqr( const int &M, const int &N, const int &K, double *A, const int &LDA, double *TAU, double *WORK, const int &LWORK, int &INFO );

  void sorgqr( const int &M, const int &N, const int &K, float *A, const int &LDA, float *TAU, float *WORK, const int &LWORK, int &INFO );

  void zgelqf( const int &M, const int &N, std::complex<double> *A, const int &LDA, std::complex<double> *TAU, std::complex<double> *WORK, const int &LWORK, int &INFO );

  void cgelqf( const int &M, const int &N, std::complex<float> *A, const int &LDA, std::complex<float> *TAU, std::complex<float> *WORK, const int &LWORK, int &INFO );

  void dgelqf( const int &M, const int &N, double *A, const int &LDA, double *TAU, double *WORK, const int &LWORK, int &INFO );

  void sgelqf( const int &M, const int &N, float *A, const int &LDA, float *TAU, float *WORK, const int &LWORK, int &INFO );


  void zunglq( const int &M, const int &N, const int &K, std::complex<double> *A, const int &LDA, std::complex<double> *TAU, std::complex<double> *WORK, const int &LWORK, int &INFO );

  void cunglq( const int &M, const int &N, const int &K, std::complex<float> *A, const int &LDA, std::complex<float> *TAU, std::complex<float> *WORK, const int &LWORK, int &INFO );

  void dorglq( const int &M, const int &N, const int &K, double *A, const int &LDA, double *TAU, double *WORK, const int &LWORK, int &INFO );

  void sorglq( const int &M, const int &N, const int &K, float *A, const int &LDA, float *TAU, float *WORK, const int &LWORK, int &INFO );


  void dger(const int* m, const int* n, const double* alpha
            , const double* x, const int* incx, const double* y, const int* incy
            , double* a, const int* lda);

  void sger(const int* m, const int* n, const float* alpha
            , const float* x, const int* incx, const float* y, const int* incy
            , float* a, const int* lda);

  void zgeru(const int* m, const int* n, const std::complex<double>* alpha
             , const std::complex<double>* x, const int* incx, const std::complex<double>* y, const int* incy
             , std::complex<double>* a, const int* lda);

  void cgeru(const int* m, const int* n, const std::complex<float>* alpha
             , const std::complex<float>* x, const int* incx, const std::complex<float>* y, const int* incy
             , std::complex<float>* a, const int* lda);

  void dormqr( const char *SIDE, const char *TRANS, const int *M, const int *N, const int *K, const double *A, const int * LDA, const double *TAU, double *C, const int *LDC, double *WORK, int *LWORK, int *INFO );

  void dgghrd(const char *COMPQ, const char *COMPZ, const int *N, const int *ILO, const int *IHI, double *A, const int *LDA, double *B, const int *LDB, double *Q, const int *LDQ, double *Z, const int *LDZ, int *INFO );

  void dhgeqz(const char *JOB, const char *COMPQ, const char *COMPZ, const int *N, const int *ILO, const int *IHI, double *H, const int *LDH, double *T, const int *LDT, double *ALPHAR, double *ALPHAI, double *BETA, double *Q, const int *LDQ, double *Z, const int *LDZ, double *WORK, int *LWORK, int *INFO);

  void dtgexc(const bool *WANTQ, const bool *WANTZ, const int *N, double *A, const int *LDA, double *B, const int *LDB, double *Q, const int *LDQ, double *Z, const int *LDZ, int *IFST, int *ILST, double *WORK, int *LWORK, int *INFO);

  void dtgevc(const char *SIDE, const char *HOWMNY, const bool *SELECT, const int *N, double *S, const int *LDS, double *P, const int *LDP, double *VL, const int *LDVL, double *VR, const int *LDVR, const int *MM, int *M, double *WORK, int *INFO );
}
#endif

