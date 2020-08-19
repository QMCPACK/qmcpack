////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Ken Esler, kpesler@gmail.com,
//    University of Illinois at Urbana-Champaign
// Miguel Morales, moralessilva2@llnl.gov,
//    Lawrence Livermore National Laboratory
// Jeongnim Kim, jeongnim.kim@gmail.com,
//    University of Illinois at Urbana-Champaign
// Jeremy McMinnis, jmcminis@gmail.com,
//    University of Illinois at Urbana-Champaign
// Mark A. Berrill, berrillma@ornl.gov,
//    Oak Ridge National Laboratory
// Alfredo A. Correa, correaa@llnl.gov
//    Lawrence Livermore National Laboratory
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
//
// File created by:
// Jeongnim Kim, jeongnim.kim@gmail.com,
//    University of Illinois at Urbana-Champaign
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_BLAS_CPU_H
#define AFQMC_BLAS_CPU_H

// generic header for blas routines
#include "AFQMC/Numerics/detail/CPU/Blasf.h"
#include "AFQMC/config.0.h"

#if defined(HAVE_MKL)
inline CBLAS_TRANSPOSE cblas_operation(char Op)
{
  if (Op == 'N')
    return CblasNoTrans;
  else if (Op == 'T')
    return CblasTrans;
  else if (Op == 'H' || Op == 'C')
    return CblasConjTrans;
  else
    throw std::runtime_error("unknown cblas_peration option");
}
#endif


namespace ma
{
/*
  static const int INCX     = 1;
  static const int INCY     = 1;
  static const char UPLO    = 'L';
  static const char TRANS   = 'T';
  static const char NOTRANS = 'N';
  static const float sone   = 1.0;
  static const float szero  = 0.0;
  static const double done  = 1.0;
  static const double dzero = 0.0;
  static const std::complex<float> cone = std::complex<float>(1.0,0.0);
  static const std::complex<float> czero = std::complex<float>(0.0,0.0);
  static const std::complex<double> zone = std::complex<double>(1.0,0.0);
  static const std::complex<double> zzero = std::complex<double>(0.0,0.0);
*/

inline static void axpy(int n, double x, const double* a, double* b) { daxpy(n, x, a, INCX, b, INCY); }

inline static void axpy(int n, double x, const double* a, int incx, double* b, int incy)
{
  daxpy(n, x, a, incx, b, incy);
}

inline static void axpy(int n, const double* a, double* b) { daxpy(n, done, a, INCX, b, INCY); }

inline static void axpy(int n, float x, const float* a, int incx, float* b, int incy) { saxpy(n, x, a, incx, b, incy); }

inline static void axpy(int n,
                        const std::complex<float> x,
                        const std::complex<float>* a,
                        int incx,
                        std::complex<float>* b,
                        int incy)
{
  caxpy(n, x, a, incx, b, incy);
}

inline static void axpy(int n,
                        const std::complex<double> x,
                        const std::complex<double>* a,
                        int incx,
                        std::complex<double>* b,
                        int incy)
{
  zaxpy(n, x, a, incx, b, incy);
}

inline static void axpy(int n, const float x, const float* a, int incx, double* b, int incy)
{
  for (int i = 0; i < n; ++i, a += incx, b += incy)
    (*b) += static_cast<double>(x * (*a));
}

inline static void axpy(int n,
                        const std::complex<float> x,
                        const std::complex<float>* a,
                        int incx,
                        std::complex<double>* b,
                        int incy)
{
  for (int i = 0; i < n; ++i, a += incx, b += incy)
    (*b) += static_cast<std::complex<double>>(x * (*a));
}

inline static double norm2(int n, const double* a, int incx = 1) { return dnrm2(n, a, incx); }

inline static double norm2(int n, const std::complex<double>* a, int incx = 1) { return dznrm2(n, a, incx); }

inline static float norm2(int n, const float* a, int incx = 1) { return snrm2(n, a, incx); }

inline static void scal(int n, float alpha, float* x, int incx = 1) { sscal(n, alpha, x, incx); }

inline static void scal(int n, std::complex<float> alpha, std::complex<float>* x, int incx = 1)
{
  cscal(n, alpha, x, incx);
}

inline static void scal(int n, double alpha, double* x, int incx = 1) { dscal(n, alpha, x, incx); }

inline static void scal(int n, std::complex<double> alpha, std::complex<double>* x, int incx = 1)
{
  zscal(n, alpha, x, incx);
}

inline static void scal(int n, double alpha, std::complex<double>* x, int incx = 1) { zdscal(n, alpha, x, incx); }

inline static void scal(int n, float alpha, std::complex<float>* x, int incx = 1) { csscal(n, alpha, x, incx); }

// inline static
// void gemv(char trans, int n, int m, const double* amat, const double* x,
// double* y) {
//  dgemv(trans, n, m, done, amat, n, x, INCX, dzero, y, INCY);
//}
inline static void gemv(int n, int m, const double* restrict amat, const double* restrict x, double* restrict y)
{
  dgemv(NOTRANS, m, n, done, amat, m, x, INCX, dzero, y, INCY);
}

inline static void gemv(int n, int m, const float* restrict amat, const float* restrict x, float* restrict y)
{
  sgemv(NOTRANS, m, n, done, amat, m, x, INCX, dzero, y, INCY);
}

inline static void gemv_trans(int n, int m, const double* restrict amat, const double* restrict x, double* restrict y)
{
  dgemv(TRANS, m, n, done, amat, m, x, INCX, dzero, y, INCY);
}

inline static void gemv_trans(int n, int m, const float* restrict amat, const float* restrict x, float* restrict y)
{
  sgemv(TRANS, m, n, done, amat, m, x, INCX, dzero, y, INCY);
}

inline static void gemv_trans(int n,
                              int m,
                              const std::complex<double>* restrict amat,
                              const std::complex<double>* restrict x,
                              std::complex<double>* restrict y)
{
  zgemv(TRANS, m, n, done, amat, m, x, INCX, dzero, y, INCY);
}

inline static void gemv_trans(int n,
                              int m,
                              const std::complex<float>* restrict amat,
                              const std::complex<float>* restrict x,
                              std::complex<float>* restrict y)
{
  cgemv(TRANS, m, n, done, amat, m, x, INCX, dzero, y, INCY);
}

inline static void gemv(int n,
                        int m,
                        const std::complex<double>* restrict amat,
                        const std::complex<double>* restrict x,
                        std::complex<double>* restrict y)
{
  zgemv(NOTRANS, m, n, zone, amat, m, x, INCX, zzero, y, INCY);
}

inline static void gemv(char trans_in,
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
  dgemv(trans_in, n, m, alpha, amat, lda, x, incx, beta, y, incy);
}

inline static void gemv(char trans_in,
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
  sgemv(trans_in, n, m, alpha, amat, lda, x, incx, beta, y, incy);
}

inline static void gemv(char trans_in,
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
  zgemv(trans_in, n, m, alpha, amat, lda, x, incx, beta, y, incy);
}

inline static void gemv(char trans_in,
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
  cgemv(trans_in, n, m, alpha, amat, lda, x, incx, beta, y, incy);
}

#if defined(HAVE_MKL)
inline static void gemv(char trans_in,
                        int n,
                        int m,
                        const double& alpha,
                        const double* restrict amat,
                        int lda,
                        const std::complex<double>* restrict x,
                        int incx,
                        const double& beta,
                        std::complex<double>* y,
                        int incy)
{
  dzgemv(trans_in, n, m, std::complex<double>(alpha), amat, lda, x, incx, std::complex<double>(beta), y, incy);
}

inline static void gemv(char trans_in,
                        int n,
                        int m,
                        const float& alpha,
                        const float* restrict amat,
                        int lda,
                        const std::complex<float>* restrict x,
                        int incx,
                        const float& beta,
                        std::complex<float>* y,
                        int incy)
{
  scgemv(trans_in, n, m, std::complex<float>(alpha), amat, lda, x, incx, std::complex<float>(beta), y, incy);
}
#else
inline static void gemv(char trans_in,
                        int n,
                        int m,
                        const double& alpha,
                        const double* restrict amat,
                        int lda,
                        const std::complex<double>* restrict x,
                        int incx,
                        const double& beta,
                        std::complex<double>* y,
                        int incy)
{
  // A * x  --(Fortran)--> A^T * x --> X * A, where X is the interpretation of x as a 2 row matrix
  // A^T * x  --(Fortran)--> A * x --> X * A^T, where X is the interpretation of x as a 2 row matrix
  if (trans_in == 'n' || trans_in == 'N')
    dgemm('N', 'T', 2, n, m, alpha, reinterpret_cast<double const*>(x), 2 * incx, amat, lda, beta,
          reinterpret_cast<double*>(y), 2 * incy);
  else if (trans_in == 't' || trans_in == 'T')
    dgemm('N', 'N', 2, m, n, alpha, reinterpret_cast<double const*>(x), 2 * incx, amat, lda, beta,
          reinterpret_cast<double*>(y), 2 * incy);
  else
  {
    print_stacktrace throw std::runtime_error("Error: Incorrect trans_in. \n");
  }
}

inline static void gemv(char trans_in,
                        int n,
                        int m,
                        const float& alpha,
                        const float* restrict amat,
                        int lda,
                        const std::complex<float>* restrict x,
                        int incx,
                        const float& beta,
                        std::complex<float>* y,
                        int incy)
{
  // A * x  --(Fortran)--> A^T * x --> X * A, where X is the interpretation of x as a 2 row matrix
  // A^T * x  --(Fortran)--> A * x --> X * A^T, where X is the interpretation of x as a 2 row matrix
  if (trans_in == 'n' || trans_in == 'N')
    sgemm('N', 'T', 2, n, m, alpha, reinterpret_cast<float const*>(x), 2 * incx, amat, lda, beta,
          reinterpret_cast<float*>(y), 2 * incy);
  else if (trans_in == 't' || trans_in == 'T')
    sgemm('N', 'N', 2, m, n, alpha, reinterpret_cast<float const*>(x), 2 * incx, amat, lda, beta,
          reinterpret_cast<float*>(y), 2 * incy);
  else
  {
    print_stacktrace throw std::runtime_error("Error: Incorrect trans_in. \n");
  }
}
#endif

inline static void gemm(char Atrans,
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
  dgemm(Atrans, Btrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

inline static void gemm(char Atrans,
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
  sgemm(Atrans, Btrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

inline static void gemm(char Atrans,
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
  zgemm(Atrans, Btrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

inline static void gemm(char Atrans,
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
  cgemm(Atrans, Btrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

inline static void gemm(char Atrans,
                        char Btrans,
                        int M,
                        int N,
                        int K,
                        double alpha,
                        const std::complex<double>* A,
                        int lda,
                        const double* restrict B,
                        int ldb,
                        double beta,
                        std::complex<double>* restrict C,
                        int ldc)
{
  assert(Atrans == 'n' || Atrans == 'N');
  dgemm(Atrans, Btrans, 2 * M, N, K, alpha, reinterpret_cast<double const*>(A), 2 * lda, B, ldb, beta,
        reinterpret_cast<double*>(C), 2 * ldc);
}

inline static void gemm(char Atrans,
                        char Btrans,
                        int M,
                        int N,
                        int K,
                        float alpha,
                        const std::complex<float>* A,
                        int lda,
                        const float* restrict B,
                        int ldb,
                        float beta,
                        std::complex<float>* restrict C,
                        int ldc)
{
  assert(Atrans == 'n' || Atrans == 'N');
  sgemm(Atrans, Btrans, 2 * M, N, K, alpha, reinterpret_cast<float const*>(A), 2 * lda, B, ldb, beta,
        reinterpret_cast<float*>(C), 2 * ldc);
}


template<typename T>
inline static T dot(int n, const T* restrict a, const T* restrict b)
{
  T res = T(0);
  for (int i = 0; i < n; ++i)
    res += a[i] * b[i];
  return res;
}

template<typename T>
inline static std::complex<T> dot(int n, const std::complex<T>* restrict a, const T* restrict b)
{
  std::complex<T> res = T(0);
  for (int i = 0; i < n; ++i)
    res += a[i] * b[i];
  return res;
}

template<typename T>
inline static std::complex<T> dot(int n, const std::complex<T>* restrict a, const std::complex<T>* restrict b)
{
  std::complex<T> res = 0.0;
  for (int i = 0; i < n; ++i)
    res += a[i] * b[i];
  return res;
}

template<typename T>
inline static std::complex<T> dot(int n, const T* restrict a, const std::complex<T>* restrict b)
{
  std::complex<T> res = 0.0;
  for (int i = 0; i < n; ++i)
    res += a[i] * b[i];
  return res;
}

template<typename T>
inline static T dot(int n, const T* restrict a, int incx, const T* restrict b, int incy)
{
  T res = T(0);
  for (int i = 0, ia = 0, ib = 0; i < n; ++i, ia += incx, ib += incy)
    res += a[ia] * b[ib];
  return res;
}

template<typename T>
inline static std::complex<T> dot(int n, const std::complex<T>* restrict a, int incx, const T* restrict b, int incy)
{
  std::complex<T> res = T(0);
  for (int i = 0, ia = 0, ib = 0; i < n; ++i, ia += incx, ib += incy)
    res += a[ia] * b[ib];
  return res;
}

template<typename T>
inline static std::complex<T> dot(int n, const T* restrict a, int incx, const std::complex<T>* restrict b, int incy)
{
  std::complex<T> res = T(0);
  for (int i = 0, ia = 0, ib = 0; i < n; ++i, ia += incx, ib += incy)
    res += a[ia] * b[ib];
  return res;
}

template<typename T>
std::complex<T> dot(int n, const std::complex<T>* a, int incx, const std::complex<T>* b, int incy)
{
  std::complex<T> res = T(0);
  for (int i = 0, ia = 0, ib = 0; i < n; ++i, ia += incx, ib += incy)
    res += a[ia] * b[ib];
  return res;
}

template<typename T>
inline static void copy(int n, const T* restrict a, T* restrict b)
{
  memcpy(b, a, sizeof(T) * n);
}

/** copy using memcpy(target,source,size)
   * @param target starting address of the targe
   * @param source starting address of the source
   * @param number of elements to copy
   */
template<typename T>
inline static void copy(T* restrict target, const T* restrict source, int n)
{
  memcpy(target, source, sizeof(T) * n);
}

template<typename T>
inline static void copy(int n, const std::complex<T>* restrict a, T* restrict b)
{
  for (int i = 0; i < n; ++i)
    b[i] = a[i].real();
}

template<typename T>
inline static void copy(int n, const T* restrict a, std::complex<T>* restrict b)
{
  for (int i = 0; i < n; ++i)
    b[i] = a[i];
}

template<typename T>
inline static void copy(int n, const T* restrict x, int incx, T* restrict y, int incy)
{
  const int xmax = incx * n;
  for (int ic = 0, jc = 0; ic < xmax; ic += incx, jc += incy)
    y[jc] = x[ic];
}

template<typename T>
inline static void copy2D(int N, int M, T const* src, int lda, T* dst, int ldb)
{
  for (int i = 0; i < N; ++i, src += lda, dst += ldb)
    copy(M, src, 1, dst, 1);
}

inline static void ger(int m,
                       int n,
                       double alpha,
                       const double* x,
                       int incx,
                       const double* y,
                       int incy,
                       double* a,
                       int lda)
{
  dger(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

inline static void ger(int m, int n, float alpha, const float* x, int incx, const float* y, int incy, float* a, int lda)
{
  sger(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

inline static void ger(int m,
                       int n,
                       const std::complex<double>& alpha,
                       const std::complex<double>* x,
                       int incx,
                       const std::complex<double>* y,
                       int incy,
                       std::complex<double>* a,
                       int lda)
{
  zgeru(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

inline static void ger(int m,
                       int n,
                       const std::complex<float>& alpha,
                       const std::complex<float>* x,
                       int incx,
                       const std::complex<float>* y,
                       int incy,
                       std::complex<float>* a,
                       int lda)
{
  cgeru(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

// Blas extensions

// C = alpha*op(A) + beta*op(B)
// unoptimized implementation
// assumes fortran ordering
template<typename T>
inline static void geam(char Atrans,
                        char Btrans,
                        int m,
                        int n,
                        T const alpha,
                        T const* A,
                        int lda,
                        T const beta,
                        T const* B,
                        int ldb,
                        T* C,
                        int ldc)
{
  if (n == 0 || m == 0)
    return;
  assert(ldc >= m);
  // Special cases
  // S1. set C to zero
  if (alpha == T(0) && beta == T(0))
  {
    T zero(0);
    for (int j = 0; j < n; j++)
      for (int i = 0; i < m; i++)
        C[i + j * ldc] = zero;
    return;
  }
  // S2.
  if (alpha == T(1) && beta == T(0))
  {
    if (std::distance<T const*>(A, C) > 0)
    {
      if (Atrans == 'N' || Atrans == 'n')
        assert(std::distance<T const*>(A, C) >= n * lda);
      else
        assert(std::distance<T const*>(A, C) >= m * lda);
    }
    else
    {
      assert(std::distance<T const*>(C, A) >= n * ldc);
    }
    if (Atrans == 'N' || Atrans == 'n')
    {
      assert(lda >= m);
      for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
          C[i + j * ldc] = A[i + j * lda];
    }
    else if (Atrans == 'T' || Atrans == 't')
    {
      assert(lda >= n);
      for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
          C[i + j * ldc] = A[j + i * lda];
    }
    else if (Atrans == 'C' || Atrans == 'c')
    {
      assert(lda >= n);
      for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
          C[i + j * ldc] = ma::conj(A[j + i * lda]);
    }
    else
    {
      throw std::runtime_error("Error: Invalid Atrans in geam.");
    }
    return;
  }
  if (alpha == T(0) && beta == T(1))
  {
    if (std::distance<T const*>(B, C) > 0)
    {
      if (Btrans == 'N' || Btrans == 'n')
        assert(std::distance<T const*>(B, C) >= n * ldb);
      else
        assert(std::distance<T const*>(B, C) >= m * ldb);
    }
    else
    {
      assert(std::distance<T const*>(C, B) >= n * ldc);
    }
    if (Btrans == 'N' || Btrans == 'n')
    {
      assert(ldb >= m);
      for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
          C[i + j * ldc] = B[i + j * ldb];
    }
    else if (Btrans == 'T' || Btrans == 't')
    {
      assert(ldb >= n);
      for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
          C[i + j * ldc] = B[j + i * ldb];
    }
    else if (Btrans == 'C' || Btrans == 'c')
    {
      assert(ldb >= n);
      for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
          C[i + j * ldc] = ma::conj(B[j + i * ldb]);
    }
    else
    {
      throw std::runtime_error("Error: Invalid Btrans in geam.");
    }
    return;
  }
  // Check for in-place modes
  // I1. A==C, lda==ldc, Atrans=='N'
  if (std::distance<T const*>(A, C) == 0)
  {
    if ((lda != ldc) || ((Atrans != 'N') && (Atrans != 'n')))
      throw std::runtime_error("Error: In-place mode requires Op(A)='N' and lda=ldc.");
    if (std::distance<T const*>(B, C) > 0)
    {
      if (Btrans == 'N' || Btrans == 'n')
        assert(std::distance<T const*>(B, C) >= n * ldb);
      else
        assert(std::distance<T const*>(B, C) >= m * ldb);
    }
    else
    {
      assert(std::distance<T const*>(C, B) >= n * ldc);
    }
    if (Btrans == 'N' || Btrans == 'n')
    {
      assert(ldb >= m);
      for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
          C[i + j * ldc] = alpha * C[i + j * ldc] + beta * B[i + j * ldb];
    }
    else if (Btrans == 'T' || Btrans == 't')
    {
      assert(ldb >= n);
      for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
          C[i + j * ldc] = alpha * C[i + j * ldc] + beta * B[j + i * ldb];
    }
    else if (Btrans == 'C' || Btrans == 'c')
    {
      assert(ldb >= n);
      for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
          C[i + j * ldc] = alpha * C[i + j * ldc] + beta * ma::conj(B[j + i * ldb]);
    }
    else
    {
      throw std::runtime_error("Error: Invalid Btrans in geam.");
    }
    return;
  }
  // I2.  B==C, ldb==ldc, Btrans=='N'
  if (std::distance<T const*>(B, C) == 0)
  {
    if ((ldb != ldc) || ((Btrans != 'N') && (Btrans != 'n')))
      throw std::runtime_error("Error: In-place mode requires Op(B)='N' and ldb=ldc.");
    if (std::distance<T const*>(A, C) > 0)
    {
      if (Atrans == 'N' || Atrans == 'n')
        assert(std::distance<T const*>(A, C) >= n * lda);
      else
        assert(std::distance<T const*>(A, C) >= m * lda);
    }
    else
    {
      assert(std::distance<T const*>(C, A) >= n * ldc);
    }
    if (Atrans == 'N' || Atrans == 'n')
    {
      assert(lda >= m);
      for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
          C[i + j * ldc] = beta * C[i + j * ldc] + alpha * A[i + j * lda];
    }
    else if (Atrans == 'T' || Atrans == 't')
    {
      assert(lda >= n);
      for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
          C[i + j * ldc] = beta * C[i + j * ldc] + alpha * A[j + i * lda];
    }
    else if (Atrans == 'C' || Atrans == 'c')
    {
      assert(lda >= n);
      for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
          C[i + j * ldc] = beta * C[i + j * ldc] + alpha * ma::conj(A[j + i * lda]);
    }
    else
    {
      throw std::runtime_error("Error: Invalid Atrans in geam.");
    }
    return;
  }
  // check that C does not overlap A or B
  if (std::distance<T const*>(A, C) > 0)
  {
    if (Atrans == 'N' || Atrans == 'n')
      assert(std::distance<T const*>(A, C) >= n * lda);
    else
      assert(std::distance<T const*>(A, C) >= m * lda);
  }
  else
  {
    assert(std::distance<T const*>(C, A) >= n * ldc);
  }
  if (std::distance<T const*>(B, C) > 0)
  {
    if (Btrans == 'N' || Btrans == 'n')
      assert(std::distance<T const*>(B, C) >= n * ldb);
    else
      assert(std::distance<T const*>(B, C) >= m * ldb);
  }
  else
  {
    assert(std::distance<T const*>(C, B) >= n * ldc);
  }
  // Generic case
  if (Atrans == 'N' || Atrans == 'n')
  {
    assert(lda >= m);
    if (Btrans == 'N' || Btrans == 'n')
    {
      assert(ldb >= m);
      for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
          C[i + j * ldc] = alpha * A[i + j * lda] + beta * B[i + j * ldb];
    }
    else if (Btrans == 'T' || Btrans == 't')
    {
      assert(ldb >= n);
      for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
          C[i + j * ldc] = alpha * A[i + j * lda] + beta * B[j + i * ldb];
    }
    else if (Btrans == 'C' || Btrans == 'c')
    {
      assert(ldb >= n);
      for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
          C[i + j * ldc] = alpha * A[i + j * lda] + beta * ma::conj(B[j + i * ldb]);
    }
    else
    {
      throw std::runtime_error("Error: Invalid Btrans in geam.");
    }
  }
  else if (Atrans == 'T' || Atrans == 't')
  {
    assert(lda >= n);
    if (Btrans == 'N' || Btrans == 'n')
    {
      assert(ldb >= m);
      for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
          C[i + j * ldc] = alpha * A[j + i * lda] + beta * B[i + j * ldb];
    }
    else if (Btrans == 'T' || Btrans == 't')
    {
      assert(ldb >= n);
      for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
          C[i + j * ldc] = alpha * A[j + i * lda] + beta * B[j + i * ldb];
    }
    else if (Btrans == 'C' || Btrans == 'c')
    {
      assert(ldb >= n);
      for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
          C[i + j * ldc] = alpha * A[j + i * lda] + beta * ma::conj(B[j + i * ldb]);
    }
    else
    {
      throw std::runtime_error("Error: Invalid Btrans in geam.");
    }
  }
  else if (Atrans == 'C' || Atrans == 'c')
  {
    assert(lda >= n);
    if (Btrans == 'N' || Btrans == 'n')
    {
      assert(ldb >= m);
      for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
          C[i + j * ldc] = alpha * ma::conj(A[j + i * lda]) + beta * B[i + j * ldb];
    }
    else if (Btrans == 'T' || Btrans == 't')
    {
      assert(ldb >= n);
      for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
          C[i + j * ldc] = alpha * ma::conj(A[j + i * lda]) + beta * B[j + i * ldb];
    }
    else if (Btrans == 'C' || Btrans == 'c')
    {
      assert(ldb >= n);
      for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
          C[i + j * ldc] = alpha * ma::conj(A[j + i * lda]) + beta * ma::conj(B[j + i * ldb]);
    }
    else
    {
      throw std::runtime_error("Error: Invalid Btrans in geam.");
    }
  }
  else
  {
    throw std::runtime_error("Error: Invalid Atrans in geam.");
  }
}

template<typename T>
inline static void zero_complex_part(int n, T* x)
{}

template<typename T>
inline static void zero_complex_part(int n, std::complex<T>* x)
{
  for (int i = 0; i < n; ++i, ++x)
    *x = std::complex<T>(real(*x), 0.0);
}

template<typename T>
inline static void set1D(int n, T const alpha, T* x, int incx)
{
  for (int i = 0; i < n; i++, x += incx)
    *x = alpha;
}

// y = beta*y + alpha * dot(a,b)
// Move to kernels
template<typename T, typename Q>
inline static void adotpby(int n,
                           T const alpha,
                           const T* restrict a,
                           int incx,
                           const T* restrict b,
                           int incy,
                           Q const beta,
                           Q* result)
{
  T res = T(0);
  for (int i = 0, ia = 0, ib = 0; i < n; ++i, ia += incx, ib += incy)
    res += a[ia] * b[ib];
  *result = beta * (*result) + static_cast<Q>(alpha * res);
}

// y[n*inc] = beta*y[n*inc] + alpha * dot(a[n*lda],b[n*lda])
template<typename T, typename Q>
inline static void strided_adotpby(int nb,
                                   int n,
                                   T const alpha,
                                   const T* restrict a,
                                   int lda,
                                   const T* restrict b,
                                   int ldb,
                                   Q const beta,
                                   Q* result,
                                   int inc)
{
  for (int k = 0; k < nb; k++)
    adotpby(n, alpha, a + k * lda, 1, b + k * ldb, 1, beta, result + k * inc);
}

/*
  template<typename T, typename Q>
  inline static
  void adotpby(int n, Q const alpha, const std::complex<T>* restrict a, int incx, const T* restrict b, int incy, Q const beta, std::complex<T>* result)
  {
    std::complex<T> res=T(0);
    for(int i=0, ia=0, ib=0; i<n; ++i, ia+=incx, ib+=incy)
      res += a[ia]*b[ib];
    *result = beta*(*result) + alpha*res;
  }

  template<typename T, typename Q>
  inline static
  void adotpby(int n, Q const alpha, const T* restrict a, int incx, const std::complex<T>* restrict b, int incy, Q const beta, std::complex<T>* result)
  {
    std::complex<T> res=T(0);
    for(int i=0, ia=0, ib=0; i<n; ++i, ia+=incx, ib+=incy)
      res += a[ia]*b[ib];
    *result = beta*(*result) + alpha*res;
  }

  template<typename T, typename Q>
  inline static
  void adotpby(int n, Q const alpha, const std::complex<T>* restrict a, int incx, const std::complex<T>* restrict b, int incy, Q const beta, std::complex<T>* result)
  {
    std::complex<T> res=T(0);
    for(int i=0, ia=0, ib=0; i<n; ++i, ia+=incx, ib+=incy)
      res += a[ia]*b[ib];
    *result = beta*(*result) + alpha*res;
  }
*/

template<typename T>
inline static void axty(int n, T const alpha, T const* x, int incx, T* y, int incy)
{
  for (int i = 0; i < n; i++)
    y[i * incy] *= alpha * x[i * incx];
}

// implements z[i][j] = beta * z[i][j] + alpha * conj(y[i][j]) * x[i]
template<typename T>
inline static void acAxpbB(int m,
                           int n,
                           T const alpha,
                           T const* A,
                           int lda,
                           T const* x,
                           int incx,
                           T const beta,
                           T* B,
                           int ldb)
{
  for (int j = 0; j < n; j++)
  {
    auto Bj = B + j * ldb;
    auto Aj = A + j * lda;
    for (int i = 0; i < m; i++)
      Bj[i] = beta * Bj[i] + alpha * ma::conj(Aj[i]) * x[i * incx];
  }
}


// y[i] = y[i] + alpha*A[i][i]
template<typename T>
inline static void adiagApy(int n, T const alpha, T const* A, int lda, T* y, int incy)
{
  for (int i = 0; i < n; i++)
    y[i * incy] += alpha * A[i * lda + i];
}

template<typename T>
inline static T sum(int n, T const* x, int incx)
{
  T res(0);
  for (int i = 0; i < n; i++)
    res += x[i * incx];
  return res;
}

// assume Fortran ordering like all blas calls
template<typename T>
inline static T sum(int m, int n, T const* A, int lda)
{
  T res(0);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      res += A[i * lda + j];
  return res;
}

// assume Fortran ordering like all blas calls
template<typename T>
void set_identity(int m, int n, T* A, int lda)
{
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      A[i * lda + j] = ((i == j) ? 1 : 0);
}

template<typename T>
void set_identity_strided(int nbatch, int stride, int m, int n, T* A, int lda)
{
  for (int b = 0, b0 = 0; b < nbatch; b++, b0 += stride)
    set_identity(m, n, A + b0, lda);
}

#ifdef HAVE_MKL

inline static void gemm_batch(const CBLAS_LAYOUT Layout,
                              const CBLAS_TRANSPOSE* transa_array,
                              const CBLAS_TRANSPOSE* transb_array,
                              const int* m_array,
                              const int* n_array,
                              const int* k_array,
                              const float* alpha_array,
                              const void** a_array,
                              const int* lda_array,
                              const void** b_array,
                              const int* ldb_array,
                              const void* beta_array,
                              void** c_array,
                              const int* ldc_array,
                              const int group_count,
                              const int* group_size)
{
  cblas_sgemm_batch(Layout, transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array,
                    b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size);
}

inline static void gemm_batch(const CBLAS_LAYOUT Layout,
                              const CBLAS_TRANSPOSE* transa_array,
                              const CBLAS_TRANSPOSE* transb_array,
                              const int* m_array,
                              const int* n_array,
                              const int* k_array,
                              const double* alpha_array,
                              const void** a_array,
                              const int* lda_array,
                              const void** b_array,
                              const int* ldb_array,
                              const void* beta_array,
                              void** c_array,
                              const int* ldc_array,
                              const int group_count,
                              const int* group_size)
{
  cblas_dgemm_batch(Layout, transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array,
                    b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size);
}

inline static void gemm_batch(const CBLAS_LAYOUT Layout,
                              const CBLAS_TRANSPOSE* transa_array,
                              const CBLAS_TRANSPOSE* transb_array,
                              const int* m_array,
                              const int* n_array,
                              const int* k_array,
                              const std::complex<float>* alpha_array,
                              const void** a_array,
                              const int* lda_array,
                              const void** b_array,
                              const int* ldb_array,
                              const void* beta_array,
                              void** c_array,
                              const int* ldc_array,
                              const int group_count,
                              const int* group_size)
{
  cblas_cgemm_batch(Layout, transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array,
                    b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size);
}

inline static void gemm_batch(const CBLAS_LAYOUT Layout,
                              const CBLAS_TRANSPOSE* transa_array,
                              const CBLAS_TRANSPOSE* transb_array,
                              const int* m_array,
                              const int* n_array,
                              const int* k_array,
                              const std::complex<double>* alpha_array,
                              const void** a_array,
                              const int* lda_array,
                              const void** b_array,
                              const int* ldb_array,
                              const void* beta_array,
                              void** c_array,
                              const int* ldc_array,
                              const int group_count,
                              const int* group_size)
{
  cblas_zgemm_batch(Layout, transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array,
                    b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size);
}

#endif

template<typename T>
inline static void gemmStridedBatched(char Atrans,
                                      char Btrans,
                                      int M,
                                      int N,
                                      int K,
                                      T alpha,
                                      const T* A,
                                      int lda,
                                      int strideA,
                                      const T* restrict B,
                                      int ldb,
                                      int strideB,
                                      T beta,
                                      T* restrict C,
                                      int ldc,
                                      int strideC,
                                      int batchSize)
{
#ifdef HAVE_MKL
  // MKL has batched gemm, but with pointer interface. Translate here
  std::vector<const void*> Aptrs(batchSize);
  std::vector<const void*> Bptrs(batchSize);
  std::vector<void*> Cptrs(batchSize);

  for (int i = 0; i < batchSize; i++)
  {
    Aptrs[i] = static_cast<const void*>(A + i * strideA);
    Bptrs[i] = static_cast<const void*>(B + i * strideB);
    Cptrs[i] = static_cast<void*>(C + i * strideC);
  }
  CBLAS_TRANSPOSE opA(cblas_operation(Atrans));
  CBLAS_TRANSPOSE opB(cblas_operation(Btrans));

  // Expect arrays of size group_count.
  // This is 1 in strided case, so passing pointers to everything
  gemm_batch(CblasColMajor, &opA, &opB, &M, &N, &K, &alpha, Aptrs.data(), &lda, Bptrs.data(), &ldb, &beta, Cptrs.data(),
             &ldc, 1, &batchSize);
#else
  // No batched gemm, :-( gemm loop
  for (int i = 0; i < batchSize; i++)
    gemm(Atrans, Btrans, M, N, K, alpha, A + i * strideA, lda, B + i * strideB, ldb, beta, C + i * strideC, ldc);
#endif
}

template<typename T,
         typename Q1,
         typename Q2,
         typename = typename std::enable_if_t<std::is_same<typename std::decay<Q1>::type, T>::value>,
         typename = typename std::enable_if_t<std::is_same<typename std::decay<Q2>::type, T>::value>>
inline static void gemmBatched(char Atrans,
                               char Btrans,
                               int M,
                               int N,
                               int K,
                               T alpha,
                               Q1** A,
                               int lda,
                               Q2** B,
                               int ldb,
                               T beta,
                               T** C,
                               int ldc,
                               int batchSize)
{
#ifdef HAVE_MKL
  CBLAS_TRANSPOSE opA(cblas_operation(Atrans));
  CBLAS_TRANSPOSE opB(cblas_operation(Btrans));
  gemm_batch(CblasColMajor, &opA, &opB, &M, &N, &K, &alpha, (const void**)A, &lda, (const void**)B, &ldb, &beta,
             (void**)C, &ldc, 1, &batchSize);
#else
  // No batched gemm, :-( gemm loop
  for (int i = 0; i < batchSize; i++)
    gemm(Atrans, Btrans, M, N, K, alpha, A[i], lda, B[i], ldb, beta, C[i], ldc);
#endif
}

//  template<typename T, typename Q>
template<typename T,
         typename Q1,
         typename Q2,
         typename T2,
         typename = typename std::enable_if_t<std::is_same<typename std::decay<Q1>::type, T2>::value>,
         typename = typename std::enable_if_t<std::is_same<typename std::decay<Q2>::type, T>::value>,
         typename = typename std::enable_if_t<std::is_same<std::complex<T>, T2>::value>>
inline static void gemmBatched(char Atrans,
                               char Btrans,
                               int M,
                               int N,
                               int K,
                               T alpha,
                               Q1** A,
                               int lda,
                               Q2** B,
                               int ldb,
                               T beta,
                               T2** C,
                               int ldc,
                               int batchSize)
{
#ifdef HAVE_MKL
  assert(Atrans == 'n' || Atrans == 'N');
  CBLAS_TRANSPOSE opA(cblas_operation(Atrans));
  CBLAS_TRANSPOSE opB(cblas_operation(Btrans));
  int M_(2 * M);
  int lda_(2 * lda);
  int ldc_(2 * ldc);
  gemm_batch(CblasColMajor, &opA, &opB, &M_, &N, &K, &alpha, (const void**)A, &lda_, (const void**)B, &ldb, &beta,
             (void**)C, &ldc_, 1, &batchSize);
#else
  // No batched gemm, :-( gemm loop
  for (int i = 0; i < batchSize; i++)
    gemm(Atrans, Btrans, M, N, K, alpha, A[i], lda, B[i], ldb, beta, C[i], ldc);
#endif
}

template<typename T>
inline static void axpyBatched(int n, T* x, T* const* a, int inca, T** b, int incb, int batchSize)
{
  for (int i = 0; i < batchSize; i++)
    axpy(n, x[i], a[i], inca, b[i], incb);
}

template<typename T>
inline static void sumGwBatched(int n, T* x, T* const* a, int inca, T** b, int incb, int b0, int nw, int batchSize)
{
  // since this is not parallel, no need to coordinate over iw
  for (int i = 0; i < batchSize; i++)
    axpy(n, x[i], a[i], inca, b[i], incb);
}

// A[k][i] = B[k][i][i]
template<typename T>
inline static void get_diagonal_strided(int nk, int ni, T const* B, int ldb, int stride, T* A, int lda)
{
  for (int k = 0; k < nk; ++k)
  {
    auto Bk(B + k * stride);
    auto Ak(A + k * lda);
    for (int i = 0; i < ni; ++i)
      Ak[i] += Bk[i * ldb + i];
  }
}

} // namespace ma

/*
#include "mpi3/shared_window.hpp"
// temporary
namespace boost {
namespace mpi3 {
namespace intranode {
template<typename T, typename Q>
auto dot(int const n, array_ptr<T> x, int const incx, array_ptr<Q> y, int const incy)
{
  using ma::dot;
  return dot(n,to_address(x),incx,to_address(y),incy);
}
} // intranode
} // mpi3
} // boost
*/
#endif
