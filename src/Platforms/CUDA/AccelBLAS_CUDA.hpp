//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
// File refactored from: MatrixDelayedUpdateCUDA.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_CUDA_ACCELBLAS_CUDA_H
#define QMCPLUSPLUS_CUDA_ACCELBLAS_CUDA_H

#include "AccelBLAS.hpp"
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
  cudaStream_t h_stream;
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
template<>
inline void gemm<PlatformKind::CUDA, float>(BLASHandle<PlatformKind::CUDA>& handle,
                                            const char transa,
                                            const char transb,
                                            int m,
                                            int n,
                                            int k,
                                            const float* alpha,
                                            const float* A,
                                            int lda,
                                            const float* B,
                                            int ldb,
                                            const float* beta,
                                            float* C,
                                            int ldc)
{
  cublasErrorCheck(cublasSgemm(handle.h_cublas, cuBLAS::convertOperation(transa), cuBLAS::convertOperation(transb), m,
                               n, k, alpha, A, lda, B, ldb, beta, C, ldc),
                   "cublasSgemm failed!");
}

template<>
inline void gemm<PlatformKind::CUDA, double>(BLASHandle<PlatformKind::CUDA>& handle,
                                             const char transa,
                                             const char transb,
                                             int m,
                                             int n,
                                             int k,
                                             const double* alpha,
                                             const double* A,
                                             int lda,
                                             const double* B,
                                             int ldb,
                                             const double* beta,
                                             double* C,
                                             int ldc)
{
  cublasErrorCheck(cublasDgemm(handle.h_cublas, cuBLAS::convertOperation(transa), cuBLAS::convertOperation(transb), m,
                               n, k, alpha, A, lda, B, ldb, beta, C, ldc),
                   "cublasDgemm failed!");
}

template<>
inline void gemm<PlatformKind::CUDA, std::complex<float>>(BLASHandle<PlatformKind::CUDA>& handle,
                                                          const char transa,
                                                          const char transb,
                                                          int m,
                                                          int n,
                                                          int k,
                                                          const std::complex<float>* alpha,
                                                          const std::complex<float>* A,
                                                          int lda,
                                                          const std::complex<float>* B,
                                                          int ldb,
                                                          const std::complex<float>* beta,
                                                          std::complex<float>* C,
                                                          int ldc)
{
  cublasErrorCheck(cublasCgemm(handle.h_cublas, cuBLAS::convertOperation(transa), cuBLAS::convertOperation(transb), m,
                               n, k, castNativeType(alpha), castNativeType(A), lda, castNativeType(B), ldb,
                               castNativeType(beta), castNativeType(C), ldc),
                   "cublasCgemm failed!");
}

template<>
inline void gemm<PlatformKind::CUDA, std::complex<double>>(BLASHandle<PlatformKind::CUDA>& handle,
                                                           const char transa,
                                                           const char transb,
                                                           int m,
                                                           int n,
                                                           int k,
                                                           const std::complex<double>* alpha,
                                                           const std::complex<double>* A,
                                                           int lda,
                                                           const std::complex<double>* B,
                                                           int ldb,
                                                           const std::complex<double>* beta,
                                                           std::complex<double>* C,
                                                           int ldc)
{
  cublasErrorCheck(cublasZgemm(handle.h_cublas, cuBLAS::convertOperation(transa), cuBLAS::convertOperation(transb), m,
                               n, k, castNativeType(alpha), castNativeType(A), lda, castNativeType(B), ldb,
                               castNativeType(beta), castNativeType(C), ldc),
                   "cublasZgemm failed!");
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

template<>
inline void gemm_batched<PlatformKind::CUDA, float>(BLASHandle<PlatformKind::CUDA>& handle,
                                                    const char transa,
                                                    const char transb,
                                                    int m,
                                                    int n,
                                                    int k,
                                                    const float* alpha,
                                                    const float* const A[],
                                                    int lda,
                                                    const float* const B[],
                                                    int ldb,
                                                    const float* beta,
                                                    float* const C[],
                                                    int ldc,
                                                    int batchCount)
{
  cublasErrorCheck(cublasSgemmBatched(handle.h_cublas, cuBLAS::convertOperation(transa),
                                      cuBLAS::convertOperation(transb), m, n, k, alpha, A, lda, B, ldb, beta, C, ldc,
                                      batchCount),
                   "cublasSgemmBatched failed!");
}

template<>
inline void gemm_batched<PlatformKind::CUDA, std::complex<float>>(BLASHandle<PlatformKind::CUDA>& handle,
                                                                  const char transa,
                                                                  const char transb,
                                                                  int m,
                                                                  int n,
                                                                  int k,
                                                                  const std::complex<float>* alpha,
                                                                  const std::complex<float>* const A[],
                                                                  int lda,
                                                                  const std::complex<float>* const B[],
                                                                  int ldb,
                                                                  const std::complex<float>* beta,
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

  cublasErrorCheck(cublasCgemmBatched(handle.h_cublas, cuBLAS::convertOperation(transa),
                                      cuBLAS::convertOperation(transb), m, n, k, castNativeType(alpha),
                                      castNativeType(non_const_A), lda, castNativeType(non_const_B), ldb,
                                      castNativeType(beta), castNativeType(non_const_C), ldc, batchCount),
                   "cublasCgemmBatched failed!");
}

template<>
inline void gemm_batched<PlatformKind::CUDA, double>(BLASHandle<PlatformKind::CUDA>& handle,
                                                     const char transa,
                                                     const char transb,
                                                     int m,
                                                     int n,
                                                     int k,
                                                     const double* alpha,
                                                     const double* const A[],
                                                     int lda,
                                                     const double* const B[],
                                                     int ldb,
                                                     const double* beta,
                                                     double* const C[],
                                                     int ldc,
                                                     int batchCount)
{
  cublasErrorCheck(cublasDgemmBatched(handle.h_cublas, cuBLAS::convertOperation(transa),
                                      cuBLAS::convertOperation(transb), m, n, k, alpha, A, lda, B, ldb, beta, C, ldc,
                                      batchCount),
                   "cublasDgemmBatched failed!");
}

template<>
inline void gemm_batched<PlatformKind::CUDA, std::complex<double>>(BLASHandle<PlatformKind::CUDA>& handle,
                                                                   const char transa,
                                                                   const char transb,
                                                                   int m,
                                                                   int n,
                                                                   int k,
                                                                   const std::complex<double>* alpha,
                                                                   const std::complex<double>* const A[],
                                                                   int lda,
                                                                   const std::complex<double>* const B[],
                                                                   int ldb,
                                                                   const std::complex<double>* beta,
                                                                   std::complex<double>* const C[],
                                                                   int ldc,
                                                                   int batchCount)
{
  auto non_const_A = const_cast<BottomConstRemoved<decltype(A)>::type>(A);
  auto non_const_B = const_cast<BottomConstRemoved<decltype(B)>::type>(B);
  auto non_const_C = const_cast<BottomConstRemoved<decltype(C)>::type>(C);

  cublasErrorCheck(cublasZgemmBatched(handle.h_cublas, cuBLAS::convertOperation(transa),
                                      cuBLAS::convertOperation(transb), m, n, k, castNativeType(alpha),
                                      castNativeType(non_const_A), lda, castNativeType(non_const_B), ldb,
                                      castNativeType(beta), castNativeType(non_const_C), ldc, batchCount),
                   "cublasZgemmBatched failed!");
}

} // namespace BLAS
} // namespace compute
} // namespace qmcplusplus
#undef castNativeType
#endif
