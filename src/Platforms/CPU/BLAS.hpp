//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_NUMERIC_BLAS_H
#define QMCPLUSPLUS_NUMERIC_BLAS_H

#include <complex>
#include <cstring>

/** Interfaces to blas library
 *
 *   static data members to facilitate /Fortran blas interface
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
namespace BLAS
{
  constexpr int INCX     = 1;
  constexpr int INCY     = 1;
  constexpr char UPLO    = 'L';
  constexpr char TRANS   = 'T';
  constexpr char NOTRANS = 'N';

  constexpr float sone                 = 1.0e0;
  constexpr float szero                = 0.0e0;
  constexpr double done                = 1.0e0;
  constexpr double dzero               = 0.0e0;
  constexpr std::complex<float> cone   = 1.0e0;
  constexpr std::complex<float> czero  = 0.0e0;
  constexpr std::complex<double> zone  = 1.0e0;
  constexpr std::complex<double> zzero = 0.0e0;

  void axpy(int n, double x, const double* a, double* b);

  void axpy(int n, double x, const double* a, int incx, double* b, int incy);

  void axpy(int n, const double* a, double* b);

  void axpy(int n, float x, const float* a, int incx, float* b, int incy);

  void axpy(int n, float x, const float* a, float* b);

  void axpy(int n,
            const std::complex<float> x,
            const std::complex<float>* a,
            int incx,
            std::complex<float>* b,
            int incy);

  void axpy(int n,
            const std::complex<double> x,
            const std::complex<double>* a,
            int incx,
            std::complex<double>* b,
            int incy);

  float norm2(int n, const float* a, int incx = 1);

  float norm2(int n, const std::complex<float>* a, int incx = 1);

  double norm2(int n, const double* a, int incx = 1);

  double norm2(int n, const std::complex<double>* a, int incx = 1);

  void scal(int n, float alpha, float* x, int incx = 1);

  void scal(int n, std::complex<float> alpha, std::complex<float>* x, int incx = 1);

  void scal(int n, double alpha, double* x, int incx = 1);

  void scal(int n, std::complex<double> alpha, std::complex<double>* x, int incx = 1);

  void scal(int n, double alpha, std::complex<double>* x, int incx = 1);

  void scal(int n, float alpha, std::complex<float>* x, int incx = 1);

  // amat is [n][m] in C
  void gemv(int n, int m, const double* restrict amat, const double* restrict x, double* restrict y);

  void gemv(int n, int m, const float* restrict amat, const float* restrict x, float* restrict y);

  void gemv(int n,
            int m,
            const std::complex<double>* restrict amat,
            const std::complex<double>* restrict x,
            std::complex<double>* restrict y);

  void gemv(int n,
            int m,
            const std::complex<float>* restrict amat,
            const std::complex<float>* restrict x,
            std::complex<float>* restrict y);

  // amat is [n][m] in C
  void gemv_trans(int n, int m, const double* restrict amat, const double* restrict x, double* restrict y);

  void gemv_trans(int n, int m, const float* restrict amat, const float* restrict x, float* restrict y);

  void gemv_trans(int n,
                  int m,
                  const std::complex<double>* restrict amat,
                  const std::complex<double>* restrict x,
                  std::complex<double>* restrict y);

  void gemv_trans(int n,
                  int m,
                  const std::complex<float>* restrict amat,
                  const std::complex<float>* restrict x,
                  std::complex<float>* restrict y);

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
            int incy);

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
            int incy);

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
            int incy);

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
            int incy);

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
            int ldc);

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
            int ldc);

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
            int ldc);

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
            int ldc);

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
  inline static std::complex<T> dot(int n,
                                    const std::complex<T>* restrict a,
                                    int incx,
                                    const std::complex<T>* restrict b,
                                    int incy)
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

  /*
    inline static
    void copy(int n, const double* restrict a, double* restrict b) {
      dcopy(n,a,1,b,1);
    }

    inline static
    void copy(int n, const std::complex<double>* restrict a, std::complex<double>* restrict b) {
      zcopy(n,a,1,b,1);
    }

    inline static
    void copy(int n, const double* restrict a, int ia, double* restrict b, int ib) {
      dcopy(n,a,ia,b,ib);
    }

    inline static
    void copy(int n, const std::complex<double>* restrict a, int ia, std::complex<double>* restrict b, int ib) {
      zcopy(n,a,ia,b,ib);
    }
  */

  void ger(int m, int n, double alpha, const double* x, int incx, const double* y, int incy, double* a, int lda);

  void ger(int m, int n, float alpha, const float* x, int incx, const float* y, int incy, float* a, int lda);

  void ger(int m,
           int n,
           const std::complex<double>& alpha,
           const std::complex<double>* x,
           int incx,
           const std::complex<double>* y,
           int incy,
           std::complex<double>* a,
           int lda);

  void ger(int m,
           int n,
           const std::complex<float>& alpha,
           const std::complex<float>* x,
           int incx,
           const std::complex<float>* y,
           int incy,
           std::complex<float>* a,
           int lda);
}; // namespace BLAS

struct LAPACK
{
  static void heev(char& jobz,
                   char& uplo,
                   int& n,
                   std::complex<float>* a,
                   int& lda,
                   float* w,
                   std::complex<float>* work,
                   int& lwork,
                   float* rwork,
                   int& info);

  static void heev(char& jobz,
                   char& uplo,
                   int& n,
                   std::complex<double>* a,
                   int& lda,
                   double* w,
                   std::complex<double>* work,
                   int& lwork,
                   double* rwork,
                   int& info);

  static void gesvd(const char& jobu,
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
                    int& info);

  static void gesvd(const char& jobu,
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
                    int& info);

  static void gesvd(const char& jobu,
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
                    int& info);

  static void gesvd(const char& jobu,
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
                    int& info);

  static void geev(char* jobvl,
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
                   int* info);

  static void geev(char* jobvl,
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
                   int* info);

  static void geev(char* jobvl,
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
                   int* info);

  static void geev(char* jobvl,
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
                   int* info);

  static void ggev(char* jobvl,
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
                   int* info);

  static void ggev(char* jobvl,
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
                   int* info);

  static void hevr(char& JOBZ,
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
                   int& INFO);

  static void hevr(char& JOBZ,
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
                   int& INFO);

  static void hevr(char& JOBZ,
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
                   int& INFO);

  static void hevr(char& JOBZ,
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
                   int& INFO);

  static void getrf(const int& n, const int& m, double* a, const int& n0, int* piv, int& st);

  static void getrf(const int& n, const int& m, float* a, const int& n0, int* piv, int& st);

  static void getrf(const int& n, const int& m, std::complex<double>* a, const int& n0, int* piv, int& st);

  static void getrf(const int& n, const int& m, std::complex<float>* a, const int& n0, int* piv, int& st);

  static void getri(int n,
                    float* restrict a,
                    int n0,
                    int const* restrict piv,
                    float* restrict work,
                    int const& n1,
                    int& status);

  static void getri(int n,
                    double* restrict a,
                    int n0,
                    int const* restrict piv,
                    double* restrict work,
                    int const& n1,
                    int& status);

  static void getri(int n,
                    std::complex<float>* restrict a,
                    int n0,
                    int const* restrict piv,
                    std::complex<float>* restrict work,
                    int const& n1,
                    int& status);

  static void getri(int n,
                    std::complex<double>* restrict a,
                    int n0,
                    int const* restrict piv,
                    std::complex<double>* restrict work,
                    int const& n1,
                    int& status);

  static void geqrf(int M,
                    int N,
                    std::complex<double>* A,
                    const int LDA,
                    std::complex<double>* TAU,
                    std::complex<double>* WORK,
                    int LWORK,
                    int& INFO);

  static void geqrf(int M, int N, double* A, const int LDA, double* TAU, double* WORK, int LWORK, int& INFO);

  static void geqrf(int M,
                    int N,
                    std::complex<float>* A,
                    const int LDA,
                    std::complex<float>* TAU,
                    std::complex<float>* WORK,
                    int LWORK,
                    int& INFO);

  static void geqrf(int M, int N, float* A, const int LDA, float* TAU, float* WORK, int LWORK, int& INFO);

  static void gelqf(int M,
                    int N,
                    std::complex<double>* A,
                    const int LDA,
                    std::complex<double>* TAU,
                    std::complex<double>* WORK,
                    int LWORK,
                    int& INFO);

  static void gelqf(int M, int N, double* A, const int LDA, double* TAU, double* WORK, int LWORK, int& INFO);

  static void gelqf(int M,
                    int N,
                    std::complex<float>* A,
                    const int LDA,
                    std::complex<float>* TAU,
                    std::complex<float>* WORK,
                    int LWORK,
                    int& INFO);

  static void gelqf(int M, int N, float* A, const int LDA, float* TAU, float* WORK, int LWORK, int& INFO);

  static void gqr(int M,
                  int N,
                  int K,
                  std::complex<double>* A,
                  const int LDA,
                  std::complex<double>* TAU,
                  std::complex<double>* WORK,
                  int LWORK,
                  int& INFO);

  static void gqr(int M, int N, int K, double* A, const int LDA, double* TAU, double* WORK, int LWORK, int& INFO);

  static void gqr(int M,
                  int N,
                  int K,
                  std::complex<float>* A,
                  const int LDA,
                  std::complex<float>* TAU,
                  std::complex<float>* WORK,
                  int LWORK,
                  int& INFO);

  static void gqr(int M, int N, int K, float* A, const int LDA, float* TAU, float* WORK, int LWORK, int& INFO);

  static void glq(int M,
                  int N,
                  int K,
                  std::complex<double>* A,
                  const int LDA,
                  std::complex<double>* TAU,
                  std::complex<double>* WORK,
                  int LWORK,
                  int& INFO);

  static void glq(int M, int N, int K, double* A, const int LDA, double* TAU, double* WORK, int LWORK, int& INFO);

  static void glq(int M,
                  int N,
                  int K,
                  std::complex<float>* A,
                  const int LDA,
                  std::complex<float>* TAU,
                  std::complex<float>* WORK,
                  int LWORK,
                  int& INFO);

  static void glq(int M, int N, int K, float* A, const int LDA, float* TAU, float* WORK, int const LWORK, int& INFO);

  static void potrf(const char& UPLO, const int& N, float* A, const int& LDA, int& INFO);

  static void potrf(const char& UPLO, const int& N, double* A, const int& LDA, int& INFO);

  static void potrf(const char& UPLO, const int& N, std::complex<float>* A, const int& LDA, int& INFO);

  static void potrf(const char& UPLO, const int& N, std::complex<double>* A, const int& LDA, int& INFO);
};


#endif // OHMMS_BLAS_H
