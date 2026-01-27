//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_CUDA_ACCELBLAS_CUDA_H
#define QMCPLUSPLUS_CUDA_ACCELBLAS_CUDA_H

#include "Common/AccelBLASHandle.hpp"
#include "CUDA/CUDAruntime.hpp"
#include "CUDA/QueueCUDA.hpp"
#include "CUDA/cuBLAS.hpp"
#include "CUDA/cuBLAS_missing_functions.hpp"

#ifndef QMC_CUDA2HIP
#define castNativeType castCUDAType
#else
#define castNativeType casthipblasType
#endif

namespace qmcplusplus
{
namespace compute
{
template<>
class BLASHandle<PlatformKind::CUDA>
{
public:
  // cuda stream, not owned, reference-only
  const cudaStream_t h_stream;
  // cublas handle
  cublasHandle_t h_cublas;

  BLASHandle(Queue<PlatformKind::CUDA>& queue) : h_stream(queue.getNative())
  {
    cublasErrorCheck(cublasCreate(&h_cublas), "cublasCreate failed!");
    cublasErrorCheck(cublasSetStream(h_cublas, h_stream), "cublasSetStream failed!");
  }

  ~BLASHandle() { cublasErrorCheck(cublasDestroy(h_cublas), "cublasDestroy failed!"); }
};

namespace BLAS
{
inline void gemm(BLASHandle<PlatformKind::CUDA>& handle,
                 const char transa,
                 const char transb,
                 int m,
                 int n,
                 int k,
                 const float& alpha,
                 const float* A,
                 int lda,
                 const float* B,
                 int ldb,
                 const float& beta,
                 float* C,
                 int ldc)
{
  cublasErrorCheck(cublasSgemm(handle.h_cublas, cuBLAS::convertOperation(transa), cuBLAS::convertOperation(transb), m,
                               n, k, &alpha, A, lda, B, ldb, &beta, C, ldc),
                   "cublasSgemm failed!");
}

inline void gemm(BLASHandle<PlatformKind::CUDA>& handle,
                 const char transa,
                 const char transb,
                 int m,
                 int n,
                 int k,
                 const double& alpha,
                 const double* A,
                 int lda,
                 const double* B,
                 int ldb,
                 const double& beta,
                 double* C,
                 int ldc)
{
  cublasErrorCheck(cublasDgemm(handle.h_cublas, cuBLAS::convertOperation(transa), cuBLAS::convertOperation(transb), m,
                               n, k, &alpha, A, lda, B, ldb, &beta, C, ldc),
                   "cublasDgemm failed!");
}

inline void gemm(BLASHandle<PlatformKind::CUDA>& handle,
                 const char transa,
                 const char transb,
                 int m,
                 int n,
                 int k,
                 const std::complex<float>& alpha,
                 const std::complex<float>* A,
                 int lda,
                 const std::complex<float>* B,
                 int ldb,
                 const std::complex<float>& beta,
                 std::complex<float>* C,
                 int ldc)
{
  const cuComplex alpha_cu = make_cuComplex(alpha.real(), alpha.imag());
  const cuComplex beta_cu  = make_cuComplex(beta.real(), beta.imag());
  cublasErrorCheck(cublasCgemm(handle.h_cublas, cuBLAS::convertOperation(transa), cuBLAS::convertOperation(transb), m,
                               n, k, &alpha_cu, castNativeType(A), lda, castNativeType(B), ldb, &beta_cu,
                               castNativeType(C), ldc),
                   "cublasCgemm failed!");
}

inline void gemm(BLASHandle<PlatformKind::CUDA>& handle,
                 const char transa,
                 const char transb,
                 int m,
                 int n,
                 int k,
                 const std::complex<double>& alpha,
                 const std::complex<double>* A,
                 int lda,
                 const std::complex<double>* B,
                 int ldb,
                 const std::complex<double>& beta,
                 std::complex<double>* C,
                 int ldc)
{
  const cuDoubleComplex alpha_cu = make_cuDoubleComplex(alpha.real(), alpha.imag());
  const cuDoubleComplex beta_cu  = make_cuDoubleComplex(beta.real(), beta.imag());
  cublasErrorCheck(cublasZgemm(handle.h_cublas, cuBLAS::convertOperation(transa), cuBLAS::convertOperation(transb), m,
                               n, k, &alpha_cu, castNativeType(A), lda, castNativeType(B), ldb, &beta_cu,
                               castNativeType(C), ldc),
                   "cublasZgemm failed!");
}

inline void gemv(BLASHandle<PlatformKind::CUDA>& handle,
                 const char trans,
                 const int m,
                 const int n,
                 const float& alpha,
                 const float* const A,
                 const int lda,
                 const float* const x,
                 const int incx,
                 const float& beta,
                 float* const y,
                 const int incy)
{
  cublasErrorCheck(cublasSgemv(handle.h_cublas, cuBLAS::convertOperation(trans), m, n, &alpha, A, lda, x, incx, &beta,
                               y, incy),
                   "cublasSgemv failed!");
}

inline void gemv(BLASHandle<PlatformKind::CUDA>& handle,
                 const char trans,
                 const int m,
                 const int n,
                 const double& alpha,
                 const double* const A,
                 const int lda,
                 const double* const x,
                 const int incx,
                 const double& beta,
                 double* const y,
                 const int incy)
{
  cublasErrorCheck(cublasDgemv(handle.h_cublas, cuBLAS::convertOperation(trans), m, n, &alpha, A, lda, x, incx, &beta,
                               y, incy),
                   "cublasDgemv failed!");
}

inline void gemv(BLASHandle<PlatformKind::CUDA>& handle,
                 const char trans,
                 const int m,
                 const int n,
                 const std::complex<float>& alpha,
                 const std::complex<float>* A,
                 const int lda,
                 const std::complex<float>* x,
                 const int incx,
                 const std::complex<float>& beta,
                 std::complex<float>* y,
                 const int incy)
{
  const cuComplex alpha_cu = make_cuComplex(alpha.real(), alpha.imag());
  const cuComplex beta_cu  = make_cuComplex(beta.real(), beta.imag());
  cublasErrorCheck(cublasCgemv(handle.h_cublas, cuBLAS::convertOperation(trans), m, n, &alpha_cu, castNativeType(A),
                               lda, castNativeType(x), incx, &beta_cu, castNativeType(y), incy),
                   "cublasCgemv failed!");
}

inline void gemv(BLASHandle<PlatformKind::CUDA>& handle,
                 const char trans,
                 const int m,
                 const int n,
                 const std::complex<double>& alpha,
                 const std::complex<double>* A,
                 const int lda,
                 const std::complex<double>* x,
                 const int incx,
                 const std::complex<double>& beta,
                 std::complex<double>* y,
                 const int incy)
{
  const cuDoubleComplex alpha_cu = make_cuDoubleComplex(alpha.real(), alpha.imag());
  const cuDoubleComplex beta_cu  = make_cuDoubleComplex(beta.real(), beta.imag());
  cublasErrorCheck(cublasZgemv(handle.h_cublas, cuBLAS::convertOperation(trans), m, n, &alpha_cu, castNativeType(A),
                               lda, castNativeType(x), incx, &beta_cu, castNativeType(y), incy),
                   "cublasZgemv failed!");
}

template<typename T>
inline void gemv_batched(BLASHandle<PlatformKind::CUDA>& handle,
                         const char trans,
                         const int m,
                         const int n,
                         const T* alpha,
                         const T* const A[],
                         const int lda,
                         const T* const x[],
                         const int incx,
                         const T* beta,
                         T* const y[],
                         const int incy,
                         const int batch_count)
{
  cudaErrorCheck(cuBLAS_MFs::gemv_batched(handle.h_stream, trans, m, n, alpha, A, lda, x, incx, beta, y, incy,
                                          batch_count),
                 "cuBLAS_MFs::gemv_batched failed!");
}

inline void ger(BLASHandle<PlatformKind::CUDA>& handle,
                const int m,
                const int n,
                const float& alpha,
                const float* const x,
                const int incx,
                const float* const y,
                const int incy,
                float* const A,
                const int lda)
{
  cublasErrorCheck(cublasSger(handle.h_cublas, m, n, &alpha, x, incx, y, incy, A, lda), "cublasSger failed!");
}

inline void ger(BLASHandle<PlatformKind::CUDA>& handle,
                const int m,
                const int n,
                const double& alpha,
                const double* const x,
                const int incx,
                const double* const y,
                const int incy,
                double* const A,
                const int lda)
{
  cublasErrorCheck(cublasDger(handle.h_cublas, m, n, &alpha, x, incx, y, incy, A, lda), "cublasDger failed!");
}

inline void ger(BLASHandle<PlatformKind::CUDA>& handle,
                const int m,
                const int n,
                const std::complex<float>& alpha,
                const std::complex<float>* x,
                const int incx,
                const std::complex<float>* y,
                const int incy,
                std::complex<float>* A,
                const int lda)
{
  const cuComplex alpha_cu = make_cuComplex(alpha.real(), alpha.imag());
  cublasErrorCheck(cublasCgeru(handle.h_cublas, m, n, &alpha_cu, castNativeType(x), incx, castNativeType(y), incy,
                               castNativeType(A), lda),
                   "cublasCger failed!");
}

inline void ger(BLASHandle<PlatformKind::CUDA>& handle,
                const int m,
                const int n,
                const std::complex<double>& alpha,
                const std::complex<double>* x,
                const int incx,
                const std::complex<double>* y,
                const int incy,
                std::complex<double>* A,
                const int lda)
{
  const cuDoubleComplex alpha_cu = make_cuDoubleComplex(alpha.real(), alpha.imag());
  cublasErrorCheck(cublasZgeru(handle.h_cublas, m, n, &alpha_cu, castNativeType(x), incx, castNativeType(y), incy,
                               castNativeType(A), lda),
                   "cublasZger failed!");
}

template<typename T>
inline void ger_batched(BLASHandle<PlatformKind::CUDA>& handle,
                        const int m,
                        const int n,
                        const T* alpha,
                        const T* const x[],
                        const int incx,
                        const T* const y[],
                        const int incy,
                        T* const A[],
                        const int lda,
                        const int batch_count)
{
  cudaErrorCheck(cuBLAS_MFs::ger_batched(handle.h_stream, m, n, alpha, x, incx, y, incy, A, lda, batch_count),
                 "cuBLAS_MFs::ger_batched failed!");
}

template<typename T>
inline void copy_batched(BLASHandle<PlatformKind::CUDA>& handle,
                         const int n,
                         const T* const in[],
                         const int incx,
                         T* const out[],
                         const int incy,
                         const int batch_count)
{
  cudaErrorCheck(cuBLAS_MFs::copy_batched(handle.h_stream, n, in, incx, out, incy, batch_count),
                 "cuBLAS_MFs::copy_batched failed!");
}

inline void gemm_batched(BLASHandle<PlatformKind::CUDA>& handle,
                         const char transa,
                         const char transb,
                         int m,
                         int n,
                         int k,
                         const float& alpha,
                         const float* const A[],
                         int lda,
                         const float* const B[],
                         int ldb,
                         const float& beta,
                         float* const C[],
                         int ldc,
                         int batchCount)
{
  cublasErrorCheck(cublasSgemmBatched(handle.h_cublas, cuBLAS::convertOperation(transa),
                                      cuBLAS::convertOperation(transb), m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc,
                                      batchCount),
                   "cublasSgemmBatched failed!");
}

inline void gemm_batched(BLASHandle<PlatformKind::CUDA>& handle,
                         const char transa,
                         const char transb,
                         int m,
                         int n,
                         int k,
                         const std::complex<float>& alpha,
                         const std::complex<float>* const A[],
                         int lda,
                         const std::complex<float>* const B[],
                         int ldb,
                         const std::complex<float>& beta,
                         std::complex<float>* const C[],
                         int ldc,
                         int batchCount)
{
  // This is necessary to not break the complex CUDA type mapping semantics while
  // dealing with the const cuComplex * A[] style API of cuBLAS
  // C++ makes you jump through some hoops to remove the bottom const on a double pointer.
  // see typetraits/type_manipulation.hpp
  auto non_const_A = const_cast<BottomConstRemoved<decltype(A)>::type>(A);
  auto non_const_B = const_cast<BottomConstRemoved<decltype(B)>::type>(B);
  auto non_const_C = const_cast<BottomConstRemoved<decltype(C)>::type>(C);

  const cuComplex alpha_cu = make_cuComplex(alpha.real(), alpha.imag());
  const cuComplex beta_cu  = make_cuComplex(beta.real(), beta.imag());

  cublasErrorCheck(cublasCgemmBatched(handle.h_cublas, cuBLAS::convertOperation(transa),
                                      cuBLAS::convertOperation(transb), m, n, k, &alpha_cu,
                                      castNativeType(non_const_A), lda, castNativeType(non_const_B), ldb,
                                      &beta_cu, castNativeType(non_const_C), ldc, batchCount),
                   "cublasCgemmBatched failed!");
}

inline void gemm_batched(BLASHandle<PlatformKind::CUDA>& handle,
                         const char transa,
                         const char transb,
                         int m,
                         int n,
                         int k,
                         const double& alpha,
                         const double* const A[],
                         int lda,
                         const double* const B[],
                         int ldb,
                         const double& beta,
                         double* const C[],
                         int ldc,
                         int batchCount)
{
  cublasErrorCheck(cublasDgemmBatched(handle.h_cublas, cuBLAS::convertOperation(transa),
                                      cuBLAS::convertOperation(transb), m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc,
                                      batchCount),
                   "cublasDgemmBatched failed!");
}

inline void gemm_batched(BLASHandle<PlatformKind::CUDA>& handle,
                         const char transa,
                         const char transb,
                         int m,
                         int n,
                         int k,
                         const std::complex<double>& alpha,
                         const std::complex<double>* const A[],
                         int lda,
                         const std::complex<double>* const B[],
                         int ldb,
                         const std::complex<double>& beta,
                         std::complex<double>* const C[],
                         int ldc,
                         int batchCount)
{
  auto non_const_A = const_cast<BottomConstRemoved<decltype(A)>::type>(A);
  auto non_const_B = const_cast<BottomConstRemoved<decltype(B)>::type>(B);
  auto non_const_C = const_cast<BottomConstRemoved<decltype(C)>::type>(C);

  const cuDoubleComplex alpha_cu = make_cuDoubleComplex(alpha.real(), alpha.imag());
  const cuDoubleComplex beta_cu  = make_cuDoubleComplex(beta.real(), beta.imag());

  cublasErrorCheck(cublasZgemmBatched(handle.h_cublas, cuBLAS::convertOperation(transa),
                                      cuBLAS::convertOperation(transb), m, n, k, &alpha_cu,
                                      castNativeType(non_const_A), lda, castNativeType(non_const_B), ldb,
                                      &beta_cu, castNativeType(non_const_C), ldc, batchCount),
                   "cublasZgemmBatched failed!");
}

} // namespace BLAS
} // namespace compute
} // namespace qmcplusplus
#undef castNativeType
#endif
