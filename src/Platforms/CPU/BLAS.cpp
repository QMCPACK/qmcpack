//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "BLAS.hpp"
#include "Blasf.h"
#include <algorithm>
#include <cassert>

namespace BLAS
{
  void axpy(int n, double x, const double* a, double* b) { daxpy(n, x, a, INCX, b, INCY); }

  void axpy(int n, double x, const double* a, int incx, double* b, int incy) { daxpy(n, x, a, incx, b, incy); }

  void axpy(int n, const double* a, double* b) { daxpy(n, done, a, INCX, b, INCY); }

  void axpy(int n, float x, const float* a, int incx, float* b, int incy) { saxpy(n, x, a, incx, b, incy); }

  void axpy(int n, float x, const float* a, float* b) { saxpy(n, x, a, INCX, b, INCY); }

  void axpy(int n,
            const std::complex<float> x,
            const std::complex<float>* a,
            int incx,
            std::complex<float>* b,
            int incy)
  {
    caxpy(n, x, a, incx, b, incy);
  }

  void axpy(int n,
            const std::complex<double> x,
            const std::complex<double>* a,
            int incx,
            std::complex<double>* b,
            int incy)
  {
    zaxpy(n, x, a, incx, b, incy);
  }

  float norm2(int n, const float* a, int incx) { return snrm2(n, a, incx); }

  float norm2(int n, const std::complex<float>* a, int incx) { return scnrm2(n, a, incx); }

  double norm2(int n, const double* a, int incx) { return dnrm2(n, a, incx); }

  double norm2(int n, const std::complex<double>* a, int incx) { return dznrm2(n, a, incx); }

  void scal(int n, float alpha, float* x, int incx) { sscal(n, alpha, x, incx); }

  void scal(int n, std::complex<float> alpha, std::complex<float>* x, int incx) { cscal(n, alpha, x, incx); }

  void scal(int n, double alpha, double* x, int incx) { dscal(n, alpha, x, incx); }

  void scal(int n, std::complex<double> alpha, std::complex<double>* x, int incx) { zscal(n, alpha, x, incx); }

  void scal(int n, double alpha, std::complex<double>* x, int incx) { zdscal(n, alpha, x, incx); }

  void scal(int n, float alpha, std::complex<float>* x, int incx) { csscal(n, alpha, x, incx); }

  void gemv(int n, int m, const double* restrict amat, const double* restrict x, double* restrict y)
  {
    dgemv(NOTRANS, m, n, done, amat, m, x, INCX, dzero, y, INCY);
  }

  void gemv(int n, int m, const float* restrict amat, const float* restrict x, float* restrict y)
  {
    sgemv(NOTRANS, m, n, done, amat, m, x, INCX, dzero, y, INCY);
  }

  void gemv(int n,
            int m,
            const std::complex<double>* restrict amat,
            const std::complex<double>* restrict x,
            std::complex<double>* restrict y)
  {
    zgemv(NOTRANS, m, n, zone, amat, m, x, INCX, zzero, y, INCY);
  }

  void gemv(int n,
            int m,
            const std::complex<float>* restrict amat,
            const std::complex<float>* restrict x,
            std::complex<float>* restrict y)
  {
    cgemv(NOTRANS, m, n, cone, amat, m, x, INCX, czero, y, INCY);
  }

  void gemv_trans(int n, int m, const double* restrict amat, const double* restrict x, double* restrict y)
  {
    dgemv(TRANS, m, n, done, amat, m, x, INCX, dzero, y, INCY);
  }

  void gemv_trans(int n, int m, const float* restrict amat, const float* restrict x, float* restrict y)
  {
    sgemv(TRANS, m, n, done, amat, m, x, INCX, dzero, y, INCY);
  }

  void gemv_trans(int n,
                  int m,
                  const std::complex<double>* restrict amat,
                  const std::complex<double>* restrict x,
                  std::complex<double>* restrict y)
  {
    zgemv(TRANS, m, n, done, amat, m, x, INCX, dzero, y, INCY);
  }

  void gemv_trans(int n,
                  int m,
                  const std::complex<float>* restrict amat,
                  const std::complex<float>* restrict x,
                  std::complex<float>* restrict y)
  {
    cgemv(TRANS, m, n, done, amat, m, x, INCX, dzero, y, INCY);
  }

  void gemv(char trans_in,
            int n,
            int m,
            double alpha,
            const double* restrict amat,
            int lda,
            const double* x,
            int incx,
            double beta,
            double* y,
            int incy)
  {
    assert(lda >= std::max(1, n));
    dgemv(trans_in, n, m, alpha, amat, lda, x, incx, beta, y, incy);
  }

  void gemv(char trans_in,
            int n,
            int m,
            float alpha,
            const float* restrict amat,
            int lda,
            const float* x,
            int incx,
            float beta,
            float* y,
            int incy)
  {
    assert(lda >= std::max(1, n));
    sgemv(trans_in, n, m, alpha, amat, lda, x, incx, beta, y, incy);
  }

  void gemv(char trans_in,
            int n,
            int m,
            const std::complex<double>& alpha,
            const std::complex<double>* restrict amat,
            int lda,
            const std::complex<double>* restrict x,
            int incx,
            const std::complex<double>& beta,
            std::complex<double>* y,
            int incy)
  {
    assert(lda >= std::max(1, n));
    zgemv(trans_in, n, m, alpha, amat, lda, x, incx, beta, y, incy);
  }

  void gemv(char trans_in,
            int n,
            int m,
            const std::complex<float>& alpha,
            const std::complex<float>* restrict amat,
            int lda,
            const std::complex<float>* restrict x,
            int incx,
            const std::complex<float>& beta,
            std::complex<float>* y,
            int incy)
  {
    assert(lda >= std::max(1, n));
    cgemv(trans_in, n, m, alpha, amat, lda, x, incx, beta, y, incy);
  }

  void gemm(char Atrans,
            char Btrans,
            int M,
            int N,
            int K,
            double alpha,
            const double* A,
            int lda,
            const double* restrict B,
            int ldb,
            double beta,
            double* restrict C,
            int ldc)
  {
    assert(lda >= std::max(1, (Atrans == 'N' || Atrans == 'n') ? M : K));
    assert(ldb >= std::max(1, (Btrans == 'N' || Btrans == 'n') ? K : N));
    assert(ldc >= std::max(1, M));
    dgemm(Atrans, Btrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  }

  void gemm(char Atrans,
            char Btrans,
            int M,
            int N,
            int K,
            float alpha,
            const float* A,
            int lda,
            const float* restrict B,
            int ldb,
            float beta,
            float* restrict C,
            int ldc)
  {
    assert(lda >= std::max(1, (Atrans == 'N' || Atrans == 'n') ? M : K));
    assert(ldb >= std::max(1, (Btrans == 'N' || Btrans == 'n') ? K : N));
    assert(ldc >= std::max(1, M));
    sgemm(Atrans, Btrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  }

  void gemm(char Atrans,
            char Btrans,
            int M,
            int N,
            int K,
            std::complex<double> alpha,
            const std::complex<double>* A,
            int lda,
            const std::complex<double>* restrict B,
            int ldb,
            std::complex<double> beta,
            std::complex<double>* restrict C,
            int ldc)
  {
    assert(lda >= std::max(1, (Atrans == 'N' || Atrans == 'n') ? M : K));
    assert(ldb >= std::max(1, (Btrans == 'N' || Btrans == 'n') ? K : N));
    assert(ldc >= std::max(1, M));
    zgemm(Atrans, Btrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  }

  void gemm(char Atrans,
            char Btrans,
            int M,
            int N,
            int K,
            std::complex<float> alpha,
            const std::complex<float>* A,
            int lda,
            const std::complex<float>* restrict B,
            int ldb,
            std::complex<float> beta,
            std::complex<float>* restrict C,
            int ldc)
  {
    assert(lda >= std::max(1, (Atrans == 'N' || Atrans == 'n') ? M : K));
    assert(ldb >= std::max(1, (Btrans == 'N' || Btrans == 'n') ? K : N));
    assert(ldc >= std::max(1, M));
    cgemm(Atrans, Btrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  }

  void ger(int m, int n, double alpha, const double* x, int incx, const double* y, int incy, double* a, int lda)
  {
    assert(lda >= std::max(1, m));
    dger(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
  }

  void ger(int m, int n, float alpha, const float* x, int incx, const float* y, int incy, float* a, int lda)
  {
    assert(lda >= std::max(1, m));
    sger(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
  }

  void ger(int m,
           int n,
           const std::complex<double>& alpha,
           const std::complex<double>* x,
           int incx,
           const std::complex<double>* y,
           int incy,
           std::complex<double>* a,
           int lda)
  {
    assert(lda >= std::max(1, m));
    zgeru(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
  }

  void ger(int m,
           int n,
           const std::complex<float>& alpha,
           const std::complex<float>* x,
           int incx,
           const std::complex<float>* y,
           int incy,
           std::complex<float>* a,
           int lda)
  {
    assert(lda >= std::max(1, m));
    cgeru(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
  }
} // namespace BLAS

void LAPACK::heev(char& jobz,
                  char& uplo,
                  int& n,
                  std::complex<float>* a,
                  int& lda,
                  float* w,
                  std::complex<float>* work,
                  int& lwork,
                  float* rwork,
                  int& info)
{
  assert(lda >= std::max(1, n));
  cheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
}

void LAPACK::heev(char& jobz,
                  char& uplo,
                  int& n,
                  std::complex<double>* a,
                  int& lda,
                  double* w,
                  std::complex<double>* work,
                  int& lwork,
                  double* rwork,
                  int& info)
{
  assert(lda >= std::max(1, n));
  zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
}

void LAPACK::gesvd(const char& jobu,
                   const char& jobvt,
                   const int& m,
                   const int& n,
                   float* a,
                   const int& lda,
                   float* s,
                   float* u,
                   const int& ldu,
                   float* vt,
                   const int& ldvt,
                   float* work,
                   const int& lwork,
                   int& info)
{
  assert(lda >= std::max(1, m));
  sgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
}

void LAPACK::gesvd(const char& jobu,
                   const char& jobvt,
                   const int& m,
                   const int& n,
                   double* a,
                   const int& lda,
                   double* s,
                   double* u,
                   const int& ldu,
                   double* vt,
                   const int& ldvt,
                   double* work,
                   const int& lwork,
                   int& info)
{
  assert(lda >= std::max(1, m));
  dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
}

void LAPACK::gesvd(const char& jobu,
                   const char& jobvt,
                   const int& m,
                   const int& n,
                   std::complex<float>* a,
                   const int& lda,
                   float* s,
                   std::complex<float>* u,
                   const int& ldu,
                   std::complex<float>* vt,
                   const int& ldvt,
                   std::complex<float>* work,
                   const int& lwork,
                   float* rwork,
                   int& info)
{
  assert(lda >= std::max(1, m));
  cgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
}

void LAPACK::gesvd(const char& jobu,
                   const char& jobvt,
                   const int& m,
                   const int& n,
                   std::complex<double>* a,
                   const int& lda,
                   double* s,
                   std::complex<double>* u,
                   const int& ldu,
                   std::complex<double>* vt,
                   const int& ldvt,
                   std::complex<double>* work,
                   const int& lwork,
                   double* rwork,
                   int& info)
{
  assert(lda >= std::max(1, m));
  zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
}

void LAPACK::geev(char* jobvl,
                  char* jobvr,
                  int* n,
                  double* a,
                  int* lda,
                  double* alphar,
                  double* alphai,
                  double* vl,
                  int* ldvl,
                  double* vr,
                  int* ldvr,
                  double* work,
                  int* lwork,
                  int* info)
{
  assert(*lda >= std::max(1, *n));
  dgeev(jobvl, jobvr, n, a, lda, alphar, alphai, vl, ldvl, vr, ldvr, work, lwork, info);
}

void LAPACK::geev(char* jobvl,
                  char* jobvr,
                  int* n,
                  float* a,
                  int* lda,
                  float* alphar,
                  float* alphai,
                  float* vl,
                  int* ldvl,
                  float* vr,
                  int* ldvr,
                  float* work,
                  int* lwork,
                  int* info)
{
  assert(*lda >= std::max(1, *n));
  sgeev(jobvl, jobvr, n, a, lda, alphar, alphai, vl, ldvl, vr, ldvr, work, lwork, info);
}

void LAPACK::geev(char* jobvl,
                  char* jobvr,
                  int* n,
                  std::complex<double>* a,
                  int* lda,
                  std::complex<double>* alpha,
                  std::complex<double>* vl,
                  int* ldvl,
                  std::complex<double>* vr,
                  int* ldvr,
                  std::complex<double>* work,
                  int* lwork,
                  double* rwork,
                  int* info)
{
  assert(*lda >= std::max(1, *n));
  zgeev(jobvl, jobvr, n, a, lda, alpha, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}

void LAPACK::geev(char* jobvl,
                  char* jobvr,
                  int* n,
                  std::complex<float>* a,
                  int* lda,
                  std::complex<float>* alpha,
                  std::complex<float>* vl,
                  int* ldvl,
                  std::complex<float>* vr,
                  int* ldvr,
                  std::complex<float>* work,
                  int* lwork,
                  float* rwork,
                  int* info)
{
  assert(*lda >= std::max(1, *n));
  cgeev(jobvl, jobvr, n, a, lda, alpha, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}

void LAPACK::ggev(char* jobvl,
                  char* jobvr,
                  int* n,
                  double* a,
                  int* lda,
                  double* b,
                  int* ldb,
                  double* alphar,
                  double* alphai,
                  double* beta,
                  double* vl,
                  int* ldvl,
                  double* vr,
                  int* ldvr,
                  double* work,
                  int* lwork,
                  int* info)
{
  assert(*lda >= std::max(1, *n));
  assert(*ldb >= std::max(1, *n));
  dggev(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info);
}

void LAPACK::ggev(char* jobvl,
                  char* jobvr,
                  int* n,
                  float* a,
                  int* lda,
                  float* b,
                  int* ldb,
                  float* alphar,
                  float* alphai,
                  float* beta,
                  float* vl,
                  int* ldvl,
                  float* vr,
                  int* ldvr,
                  float* work,
                  int* lwork,
                  int* info)
{
  assert(*lda >= std::max(1, *n));
  assert(*ldb >= std::max(1, *n));
  sggev(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info);
}

void LAPACK::hevr(char& JOBZ,
                  char& RANGE,
                  char& UPLO,
                  int& N,
                  float* A,
                  int& LDA,
                  float& VL,
                  float& VU,
                  int& IL,
                  int& IU,
                  float& ABSTOL,
                  int& M,
                  float* W,
                  float* Z,
                  int& LDZ,
                  int* ISUPPZ,
                  float* WORK,
                  int& LWORK,
                  float* RWORK,
                  int& LRWORK,
                  int* IWORK,
                  int& LIWORK,
                  int& INFO)
{
  assert(LDA >= std::max(1, N));
  if (WORK)
    WORK[0] = 0;
  ssyevr(JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, RWORK, LRWORK, IWORK, LIWORK,
         INFO);
}

void LAPACK::hevr(char& JOBZ,
                  char& RANGE,
                  char& UPLO,
                  int& N,
                  double* A,
                  int& LDA,
                  double& VL,
                  double& VU,
                  int& IL,
                  int& IU,
                  double& ABSTOL,
                  int& M,
                  double* W,
                  double* Z,
                  int& LDZ,
                  int* ISUPPZ,
                  double* WORK,
                  int& LWORK,
                  double* RWORK,
                  int& LRWORK,
                  int* IWORK,
                  int& LIWORK,
                  int& INFO)
{
  assert(LDA >= std::max(1, N));
  if (WORK)
    WORK[0] = 0;
  dsyevr(JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, RWORK, LRWORK, IWORK, LIWORK,
         INFO);
}

void LAPACK::hevr(char& JOBZ,
                  char& RANGE,
                  char& UPLO,
                  int& N,
                  std::complex<float>* A,
                  int& LDA,
                  float& VL,
                  float& VU,
                  int& IL,
                  int& IU,
                  float& ABSTOL,
                  int& M,
                  float* W,
                  std::complex<float>* Z,
                  int& LDZ,
                  int* ISUPPZ,
                  std::complex<float>* WORK,
                  int& LWORK,
                  float* RWORK,
                  int& LRWORK,
                  int* IWORK,
                  int& LIWORK,
                  int& INFO)
{
  assert(LDA >= std::max(1, N));
  cheevr(JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, RWORK, LRWORK,
         IWORK, LIWORK, INFO);
}

void LAPACK::hevr(char& JOBZ,
                  char& RANGE,
                  char& UPLO,
                  int& N,
                  std::complex<double>* A,
                  int& LDA,
                  double& VL,
                  double& VU,
                  int& IL,
                  int& IU,
                  double& ABSTOL,
                  int& M,
                  double* W,
                  std::complex<double>* Z,
                  int& LDZ,
                  int* ISUPPZ,
                  std::complex<double>* WORK,
                  int& LWORK,
                  double* RWORK,
                  int& LRWORK,
                  int* IWORK,
                  int& LIWORK,
                  int& INFO)
{
  assert(LDA >= std::max(1, N));
  zheevr(JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, RWORK, LRWORK,
         IWORK, LIWORK, INFO);
}

void LAPACK::getrf(const int& n, const int& m, double* a, const int& n0, int* piv, int& st)
{
  dgetrf(n, m, a, n0, piv, st);
}

void LAPACK::getrf(const int& n, const int& m, float* a, const int& n0, int* piv, int& st)
{
  sgetrf(n, m, a, n0, piv, st);
}

void LAPACK::getrf(const int& n, const int& m, std::complex<double>* a, const int& n0, int* piv, int& st)
{
  zgetrf(n, m, a, n0, piv, st);
}

void LAPACK::getrf(const int& n, const int& m, std::complex<float>* a, const int& n0, int* piv, int& st)
{
  cgetrf(n, m, a, n0, piv, st);
}

void LAPACK::getri(int n,
                   float* restrict a,
                   int n0,
                   int const* restrict piv,
                   float* restrict work,
                   int const& n1,
                   int& status)
{
  sgetri(n, a, n0, piv, work, n1, status);
}

void LAPACK::getri(int n,
                   double* restrict a,
                   int n0,
                   int const* restrict piv,
                   double* restrict work,
                   int const& n1,
                   int& status)
{
  dgetri(n, a, n0, piv, work, n1, status);
}

void LAPACK::getri(int n,
                   std::complex<float>* restrict a,
                   int n0,
                   int const* restrict piv,
                   std::complex<float>* restrict work,
                   int const& n1,
                   int& status)
{
  cgetri(n, a, n0, piv, work, n1, status);
}

void LAPACK::getri(int n,
                   std::complex<double>* restrict a,
                   int n0,
                   int const* restrict piv,
                   std::complex<double>* restrict work,
                   int const& n1,
                   int& status)
{
  zgetri(n, a, n0, piv, work, n1, status);
}

void LAPACK::geqrf(int M,
                   int N,
                   std::complex<double>* A,
                   const int LDA,
                   std::complex<double>* TAU,
                   std::complex<double>* WORK,
                   int LWORK,
                   int& INFO)
{
  assert(LDA >= std::max(1, M));
  zgeqrf(M, N, A, LDA, TAU, WORK, LWORK, INFO);
}

void LAPACK::geqrf(int M, int N, double* A, const int LDA, double* TAU, double* WORK, int LWORK, int& INFO)
{
  assert(LDA >= std::max(1, M));
  dgeqrf(M, N, A, LDA, TAU, WORK, LWORK, INFO);
}

void LAPACK::geqrf(int M,
                   int N,
                   std::complex<float>* A,
                   const int LDA,
                   std::complex<float>* TAU,
                   std::complex<float>* WORK,
                   int LWORK,
                   int& INFO)
{
  assert(LDA >= std::max(1, M));
  cgeqrf(M, N, A, LDA, TAU, WORK, LWORK, INFO);
}

void LAPACK::geqrf(int M, int N, float* A, const int LDA, float* TAU, float* WORK, int LWORK, int& INFO)
{
  assert(LDA >= std::max(1, M));
  sgeqrf(M, N, A, LDA, TAU, WORK, LWORK, INFO);
}

void LAPACK::gelqf(int M,
                   int N,
                   std::complex<double>* A,
                   const int LDA,
                   std::complex<double>* TAU,
                   std::complex<double>* WORK,
                   int LWORK,
                   int& INFO)
{
  assert(LDA >= std::max(1, M));
  zgelqf(M, N, A, LDA, TAU, WORK, LWORK, INFO);
}

void LAPACK::gelqf(int M, int N, double* A, const int LDA, double* TAU, double* WORK, int LWORK, int& INFO)
{
  assert(LDA >= std::max(1, M));
  dgelqf(M, N, A, LDA, TAU, WORK, LWORK, INFO);
}

void LAPACK::gelqf(int M,
                   int N,
                   std::complex<float>* A,
                   const int LDA,
                   std::complex<float>* TAU,
                   std::complex<float>* WORK,
                   int LWORK,
                   int& INFO)
{
  assert(LDA >= std::max(1, M));
  cgelqf(M, N, A, LDA, TAU, WORK, LWORK, INFO);
}

void LAPACK::gelqf(int M, int N, float* A, const int LDA, float* TAU, float* WORK, int LWORK, int& INFO)
{
  assert(LDA >= std::max(1, M));
  sgelqf(M, N, A, LDA, TAU, WORK, LWORK, INFO);
}

void LAPACK::gqr(int M,
                 int N,
                 int K,
                 std::complex<double>* A,
                 const int LDA,
                 std::complex<double>* TAU,
                 std::complex<double>* WORK,
                 int LWORK,
                 int& INFO)
{
  assert(LDA >= std::max(1, M));
  zungqr(M, N, K, A, LDA, TAU, WORK, LWORK, INFO);
}

void LAPACK::gqr(int M, int N, int K, double* A, const int LDA, double* TAU, double* WORK, int LWORK, int& INFO)
{
  assert(LDA >= std::max(1, M));
  dorgqr(M, N, K, A, LDA, TAU, WORK, LWORK, INFO);
}

void LAPACK::gqr(int M,
                 int N,
                 int K,
                 std::complex<float>* A,
                 const int LDA,
                 std::complex<float>* TAU,
                 std::complex<float>* WORK,
                 int LWORK,
                 int& INFO)
{
  assert(LDA >= std::max(1, M));
  cungqr(M, N, K, A, LDA, TAU, WORK, LWORK, INFO);
}

void LAPACK::gqr(int M, int N, int K, float* A, const int LDA, float* TAU, float* WORK, int LWORK, int& INFO)
{
  assert(LDA >= std::max(1, M));
  sorgqr(M, N, K, A, LDA, TAU, WORK, LWORK, INFO);
}

void LAPACK::glq(int M,
                 int N,
                 int K,
                 std::complex<double>* A,
                 const int LDA,
                 std::complex<double>* TAU,
                 std::complex<double>* WORK,
                 int LWORK,
                 int& INFO)
{
  assert(LDA >= std::max(1, M));
  zunglq(M, N, K, A, LDA, TAU, WORK, LWORK, INFO);
}

void LAPACK::glq(int M, int N, int K, double* A, const int LDA, double* TAU, double* WORK, int LWORK, int& INFO)
{
  assert(LDA >= std::max(1, M));
  dorglq(M, N, K, A, LDA, TAU, WORK, LWORK, INFO);
}

void LAPACK::glq(int M,
                 int N,
                 int K,
                 std::complex<float>* A,
                 const int LDA,
                 std::complex<float>* TAU,
                 std::complex<float>* WORK,
                 int LWORK,
                 int& INFO)
{
  assert(LDA >= std::max(1, M));
  cunglq(M, N, K, A, LDA, TAU, WORK, LWORK, INFO);
}

void LAPACK::glq(int M, int N, int K, float* A, const int LDA, float* TAU, float* WORK, int const LWORK, int& INFO)
{
  assert(LDA >= std::max(1, M));
  sorglq(M, N, K, A, LDA, TAU, WORK, LWORK, INFO);
}

void LAPACK::potrf(const char& UPLO, const int& N, float* A, const int& LDA, int& INFO)
{
  assert(LDA >= std::max(1, N));
  spotrf(UPLO, N, A, LDA, INFO);
}

void LAPACK::potrf(const char& UPLO, const int& N, double* A, const int& LDA, int& INFO)
{
  assert(LDA >= std::max(1, N));
  dpotrf(UPLO, N, A, LDA, INFO);
}

void LAPACK::potrf(const char& UPLO, const int& N, std::complex<float>* A, const int& LDA, int& INFO)
{
  assert(LDA >= std::max(1, N));
  cpotrf(UPLO, N, A, LDA, INFO);
}

void LAPACK::potrf(const char& UPLO, const int& N, std::complex<double>* A, const int& LDA, int& INFO)
{
  assert(LDA >= std::max(1, N));
  zpotrf(UPLO, N, A, LDA, INFO);
}
