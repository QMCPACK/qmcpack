//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "syclBLAS.hpp"
#include "oneapi/mkl/blas.hpp"

namespace qmcplusplus
{
namespace syclBLAS
{
inline oneapi::mkl::transpose convertTransEnum(char trans)
{
  return trans == 'T' ? oneapi::mkl::transpose::trans : oneapi::mkl::transpose::nontrans;
}

template<typename T>
sycl::event gemv(sycl::queue& handle,
                 const char trans,
                 const int m,
                 const int n,
                 const T alpha,
                 const T* const A,
                 const int lda,
                 const T* const x,
                 const int incx,
                 const T beta,
                 T* const y,
                 const int incy,
                 const std::vector<sycl::event>& events)
{
  return oneapi::mkl::blas::gemv(handle, convertTransEnum(trans), m, n, alpha, A, lda, x, incx, beta, y, incy, events);
}

template sycl::event gemv(sycl::queue& handle,
                          const char trans,
                          const int m,
                          const int n,
                          const double alpha,
                          const double* const A,
                          const int lda,
                          const double* const x,
                          const int incx,
                          const double beta,
                          double* const y,
                          const int incy,
                          const std::vector<sycl::event>& events);

template sycl::event gemv(sycl::queue& handle,
                          const char trans,
                          const int m,
                          const int n,
                          const float alpha,
                          const float* const A,
                          const int lda,
                          const float* const x,
                          const int incx,
                          const float beta,
                          float* const y,
                          const int incy,
                          const std::vector<sycl::event>& events);

template sycl::event gemv(sycl::queue& handle,
                          const char trans,
                          const int m,
                          const int n,
                          const std::complex<double> alpha,
                          const std::complex<double>* const A,
                          const int lda,
                          const std::complex<double>* const x,
                          const int incx,
                          const std::complex<double> beta,
                          std::complex<double>* const y,
                          const int incy,
                          const std::vector<sycl::event>& events);

template sycl::event gemv(sycl::queue& handle,
                          const char trans,
                          const int m,
                          const int n,
                          const std::complex<float> alpha,
                          const std::complex<float>* const A,
                          const int lda,
                          const std::complex<float>* const x,
                          const int incx,
                          const std::complex<float> beta,
                          std::complex<float>* const y,
                          const int incy,
                          const std::vector<sycl::event>& events);

template<typename T>
sycl::event gemm(sycl::queue& handle,
                 const char tA,
                 const char tB,
                 const int m,
                 const int n,
                 const int k,
                 const T alpha,
                 const T* A,
                 const int lda,
                 const T* B,
                 const int ldb,
                 const T beta,
                 T* C,
                 const int ldc,
                 const std::vector<sycl::event>& events)
{
   return oneapi::mkl::blas::gemm(handle, convertTransEnum(tA), convertTransEnum(tB), m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, events);
}



template sycl::event gemm(sycl::queue& handle,
                          const char tA,
                          const char tB,
                          const int m,
                          const int n,
                          const int k,
                          const float alpha,
                          const float* const A,
                          const int lda,
                          const float* const B,
                          const int ldb,
                          const float beta,
                          float* const C,
                          const int ldc,
                          const std::vector<sycl::event>& events);

template sycl::event gemm(sycl::queue& handle,
                          const char tA,
                          const char tB,
                          const int m,
                          const int n,
                          const int k,
                          const double alpha,
                          const double* const A,
                          const int lda,
                          const double* const B,
                          const int ldb,
                          const double beta,
                          double* const C,
                          const int ldc,
                          const std::vector<sycl::event>& events);

template sycl::event gemm(sycl::queue& handle,
                          const char tA,
                          const char tB,
                          const int m,
                          const int n,
                          const int k,
                          const std::complex<float> alpha,
                          const std::complex<float>* const A,
                          const int lda,
                          const std::complex<float>* const B,
                          const int ldb,
                          const std::complex<float> beta,
                          std::complex<float>* const C,
                          const int ldc,
                          const std::vector<sycl::event>& events);

template sycl::event gemm(sycl::queue& handle,
                          const char tA,
                          const char tB,
                          const int m,
                          const int n,
                          const int k,
                          const std::complex<double> alpha,
                          const std::complex<double>* const A,
                          const int lda,
                          const std::complex<double>* const B,
                          const int ldb,
                          const std::complex<double> beta,
                          std::complex<double>* const C,
                          const int ldc,
                          const std::vector<sycl::event>& events);

} // namespace syclBLAS

} // namespace qmcplusplus
