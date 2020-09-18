//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_BLAS_CUDA_HPP
#define AFQMC_BLAS_CUDA_HPP

#include <cassert>
#include <vector>
#include "AFQMC/Memory/CUDA/cuda_gpu_pointer.hpp"
#include "AFQMC/Numerics/detail/CUDA/cublas_wrapper.hpp"
#include "AFQMC/Numerics/detail/CUDA/cublasXt_wrapper.hpp"
// hand coded kernels for blas extensions
#include "AFQMC/Numerics/detail/CUDA/Kernels/adotpby.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/axty.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/sum.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/adiagApy.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/acAxpbB.cuh"

// Currently available:
// Lvl-1: dot, axpy, scal
// Lvl-2: gemv
// Lvl-3: gemm

namespace qmc_cuda
{
// copy Specializations
template<class ptr, typename = typename std::enable_if_t<(ptr::memory_type != CPU_OUTOFCARS_POINTER_TYPE)>>
inline static void copy(int n, ptr x, int incx, ptr y, int incy)
{
  if (CUBLAS_STATUS_SUCCESS !=
      cublas::cublas_copy(*x.handles.cublas_handle, n, to_address(x), incx, to_address(y), incy))
    throw std::runtime_error("Error: cublas_copy returned error code.");
}

template<class ptr,
         typename = typename std::enable_if_t<(ptr::memory_type == CPU_OUTOFCARS_POINTER_TYPE)>,
         typename = void>
inline static void copy(int n, ptr x, int incx, ptr y, int incy)
{
  using ma::copy;
  return copy(n, to_address(x), incx, to_address(y), incy);
}

// scal Specializations
template<class T, class ptr, typename = typename std::enable_if_t<(ptr::memory_type != CPU_OUTOFCARS_POINTER_TYPE)>>
inline static void scal(int n, T alpha, ptr x, int incx)
{
  if (CUBLAS_STATUS_SUCCESS != cublas::cublas_scal(*x.handles.cublas_handle, n, alpha, to_address(x), incx))
    throw std::runtime_error("Error: cublas_scal returned error code.");
}

template<class T,
         class ptr,
         typename = typename std::enable_if_t<(ptr::memory_type == CPU_OUTOFCARS_POINTER_TYPE)>,
         typename = void>
inline static void scal(int n, T alpha, ptr x, int incx)
{
  using ma::scal;
  return scal(n, alpha, to_address(x), incx);
}

// dot Specializations
template<class ptrA,
         class ptrB,
         typename = typename std::enable_if_t<(ptrA::memory_type != CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrB::memory_type != CPU_OUTOFCARS_POINTER_TYPE)>>
inline static auto dot(int const n, ptrA const& x, int const incx, ptrB const& y, int const incy)
{
  return cublas::cublas_dot(*x.handles.cublas_handle, n, to_address(x), incx, to_address(y), incy);
}

template<class ptrA,
         class ptrB,
         typename = typename std::enable_if_t<(ptrA::memory_type == CPU_OUTOFCARS_POINTER_TYPE) or
                                              (ptrB::memory_type == CPU_OUTOFCARS_POINTER_TYPE)>,
         typename = void>
inline static auto dot(int const n, ptrA const& x, int const incx, ptrB const& y, int const incy)
{
  using ma::dot;
  return dot(n, to_address(x), incx, to_address(y), incy);
}

// axpy Specializations
template<typename T,
         class ptrA,
         class ptrB,
         typename = typename std::enable_if_t<(ptrA::memory_type != CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrB::memory_type != CPU_OUTOFCARS_POINTER_TYPE)>>
inline static void axpy(int n, T const a, ptrA const& x, int incx, ptrB&& y, int incy)
{
  if (CUBLAS_STATUS_SUCCESS !=
      cublas::cublas_axpy(*x.handles.cublas_handle, n, a, to_address(x), incx, to_address(y), incy))
    throw std::runtime_error("Error: cublas_axpy returned error code.");
}

template<typename T,
         class ptrA,
         class ptrB,
         typename = typename std::enable_if_t<(ptrA::memory_type == CPU_OUTOFCARS_POINTER_TYPE) or
                                              (ptrB::memory_type == CPU_OUTOFCARS_POINTER_TYPE)>,
         typename = void>
inline static void axpy(int n, T const a, ptrA const& x, int incx, ptrB&& y, int incy)
{
  using ma::axpy;
  axpy(n, a, to_address(x), incx, to_address(y), incy);
}

// GEMV Specializations
template<typename T,
         class ptrA,
         class ptrB,
         class ptrC,
         typename = typename std::enable_if_t<(ptrA::memory_type != CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrB::memory_type != CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrC::memory_type != CPU_OUTOFCARS_POINTER_TYPE)>>
inline static void gemv(char Atrans,
                        int M,
                        int N,
                        T alpha,
                        ptrA const& A,
                        int lda,
                        ptrB const& x,
                        int incx,
                        T beta,
                        ptrC&& y,
                        int incy)
{
  if (CUBLAS_STATUS_SUCCESS !=
      cublas::cublas_gemv(*A.handles.cublas_handle, Atrans, M, N, alpha, to_address(A), lda, to_address(x), incx, beta,
                          to_address(y), incy))
    throw std::runtime_error("Error: cublas_gemv returned error code.");
}

template<typename T,
         class ptrA,
         class ptrB,
         class ptrC,
         typename = typename std::enable_if_t<(ptrA::memory_type == CPU_OUTOFCARS_POINTER_TYPE) or
                                              (ptrB::memory_type == CPU_OUTOFCARS_POINTER_TYPE) or
                                              (ptrC::memory_type == CPU_OUTOFCARS_POINTER_TYPE)>,
         typename = void>
inline static void gemv(char Atrans,
                        int M,
                        int N,
                        T alpha,
                        ptrA const& A,
                        int lda,
                        ptrB const& x,
                        int incx,
                        T beta,
                        ptrC&& y,
                        int incy)
{
  using ma::gemv;
  gemv(Atrans, M, N, alpha, to_address(A), lda, to_address(x), incx, beta, to_address(y), incy);
  /*
    const char Btrans('N');
    const int one(1);
    if(CUBLAS_STATUS_SUCCESS != cublas::cublasXt_gemm(*A.handles.cublasXt_handle,Atrans,Btrans,
                                            M,one,K,alpha,to_address(A),lda,to_address(x),incx,
                                            beta,to_address(y),incy))
      throw std::runtime_error("Error: cublasXt_gemv (gemm) returned error code.");
*/
}

// GEMM Specializations
template<typename T,
         class ptrA,
         class ptrB,
         class ptrC,
         typename = typename std::enable_if_t<(ptrA::memory_type != CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrB::memory_type != CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrC::memory_type != CPU_OUTOFCARS_POINTER_TYPE)>>
inline static void gemm(char Atrans,
                        char Btrans,
                        int M,
                        int N,
                        int K,
                        T alpha,
                        ptrA const& A,
                        int lda,
                        ptrB const& B,
                        int ldb,
                        T beta,
                        ptrC&& C,
                        int ldc)
{
  if (CUBLAS_STATUS_SUCCESS !=
      cublas::cublas_gemm(*A.handles.cublas_handle, Atrans, Btrans, M, N, K, alpha, to_address(A), lda, to_address(B),
                          ldb, beta, to_address(C), ldc))
    throw std::runtime_error("Error: cublas_gemm returned error code.");
}

template<typename T,
         class ptrA,
         class ptrB,
         class ptrC,
         typename = typename std::enable_if_t<(ptrA::memory_type == CPU_OUTOFCARS_POINTER_TYPE) or
                                              (ptrB::memory_type == CPU_OUTOFCARS_POINTER_TYPE) or
                                              (ptrC::memory_type == CPU_OUTOFCARS_POINTER_TYPE)>,
         typename = void>
inline static void gemm(char Atrans,
                        char Btrans,
                        int M,
                        int N,
                        int K,
                        T alpha,
                        ptrA const& A,
                        int lda,
                        ptrB const& B,
                        int ldb,
                        T beta,
                        ptrC&& C,
                        int ldc)
{
  if (CUBLAS_STATUS_SUCCESS !=
      cublas::cublasXt_gemm(*A.handles.cublasXt_handle, Atrans, Btrans, M, N, K, alpha, to_address(A), lda,
                            to_address(B), ldb, beta, to_address(C), ldc))
    throw std::runtime_error("Error: cublasXt_gemm returned error code.");
}

// Blas Extensions
// geam
template<class T,
         class ptrA,
         class ptrB,
         class ptrC,
         typename = typename std::enable_if_t<(ptrA::memory_type != CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrB::memory_type != CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrC::memory_type != CPU_OUTOFCARS_POINTER_TYPE)>>
inline static void geam(char Atrans,
                        char Btrans,
                        int M,
                        int N,
                        T const alpha,
                        ptrA const& A,
                        int lda,
                        T const beta,
                        ptrB const& B,
                        int ldb,
                        ptrC C,
                        int ldc)
{
  if (CUBLAS_STATUS_SUCCESS !=
      cublas::cublas_geam(*A.handles.cublas_handle, Atrans, Btrans, M, N, alpha, to_address(A), lda, beta,
                          to_address(B), ldb, to_address(C), ldc))
    throw std::runtime_error("Error: cublas_geam returned error code.");
}

template<class T,
         class ptrA,
         class ptrB,
         class ptrC,
         typename = typename std::enable_if_t<(ptrA::memory_type == CPU_OUTOFCARS_POINTER_TYPE) or
                                              (ptrB::memory_type == CPU_OUTOFCARS_POINTER_TYPE) or
                                              (ptrC::memory_type == CPU_OUTOFCARS_POINTER_TYPE)>,
         typename = void>
inline static void geam(char Atrans,
                        char Btrans,
                        int M,
                        int N,
                        T const alpha,
                        ptrA const& A,
                        int lda,
                        T const beta,
                        ptrB const& B,
                        int ldb,
                        ptrC C,
                        int ldc)
{
  using ma::geam;
  return geam(Atrans, Btrans, M, N, alpha, to_address(A), lda, beta, to_address(B), ldb, to_address(C), ldc);
}

//template<class T,
template<class ptr,
         typename = typename std::enable_if_t<not(ptr::memory_type == CPU_OUTOFCARS_POINTER_TYPE)>,
         typename = void>
//inline static void set1D(int n, T const alpha, ptr x, int incx)
inline static void set1D(int n, typename ptr::value_type const alpha, ptr x, int incx)
{
  // No set funcion in cuda!!! Avoiding kernels for now
  std::vector<typename ptr::value_type> buff(n, alpha);
  if (CUBLAS_STATUS_SUCCESS !=
      cublasSetVector(n, sizeof(typename ptr::value_type), buff.data(), 1, to_address(x), incx))
    throw std::runtime_error("Error: cublasSetVector returned error code.");
}

template<class T, class ptr, typename = typename std::enable_if_t<(ptr::memory_type == CPU_OUTOFCARS_POINTER_TYPE)>>
inline static void set1D(int n, T const alpha, ptr x, int incx)
{
  auto y = to_address(x);
  for (int i = 0; i < n; i++, y += incx)
    *y = alpha;
}

// dot extension
template<class T,
         class Q,
         class ptrA,
         class ptrB,
         class ptrC,
         typename = typename std::enable_if_t<(ptrA::memory_type != CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrB::memory_type != CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrC::memory_type != CPU_OUTOFCARS_POINTER_TYPE)>>
inline static void adotpby(int const n,
                           T const alpha,
                           ptrA const& x,
                           int const incx,
                           ptrB const& y,
                           int const incy,
                           Q const beta,
                           ptrC result)
{
  kernels::adotpby(n, alpha, to_address(x), incx, to_address(y), incy, beta, to_address(result));
}

template<class T,
         class Q,
         class ptrA,
         class ptrB,
         class ptrC,
         typename = typename std::enable_if_t<(ptrA::memory_type == CPU_OUTOFCARS_POINTER_TYPE) or
                                              (ptrB::memory_type == CPU_OUTOFCARS_POINTER_TYPE) or
                                              (ptrC::memory_type == CPU_OUTOFCARS_POINTER_TYPE)>,
         typename = void>
inline static void adotpby(int const n,
                           T const alpha,
                           ptrA const& x,
                           int const incx,
                           ptrB const& y,
                           int const incy,
                           Q const beta,
                           ptrC result)
{
  using ma::adotpby;
  adotpby(n, alpha, to_address(x), incx, to_address(y), incy, beta, to_address(result));
}


// axty
template<class T,
         class ptrA,
         class ptrB,
         typename = typename std::enable_if_t<(ptrA::memory_type != CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrB::memory_type != CPU_OUTOFCARS_POINTER_TYPE)>>
inline static void axty(int n, T const alpha, ptrA const x, int incx, ptrB y, int incy)
{
  if (incx != 1 || incy != 1)
    throw std::runtime_error("Error: axty with inc != 1 not implemented.");
  kernels::axty(n, alpha, to_address(x), to_address(y));
}

template<class T,
         class ptrA,
         class ptrB,
         typename = typename std::enable_if_t<(ptrA::memory_type == CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrB::memory_type == CPU_OUTOFCARS_POINTER_TYPE)>,
         typename = void>
inline static void axty(int n, T const alpha, ptrA const x, int incx, ptrB y, int incy)
{
  using ma::axty;
  axty(n, alpha, to_address(x), incx, to_address(y), incy);
}

// acAxpbB
template<class T,
         class ptrA,
         class ptrx,
         class ptrB,
         typename = typename std::enable_if_t<(ptrA::memory_type != CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrx::memory_type != CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrB::memory_type != CPU_OUTOFCARS_POINTER_TYPE)>>
inline static void acAxpbB(int m,
                           int n,
                           T const alpha,
                           ptrA const A,
                           int lda,
                           ptrx const x,
                           int incx,
                           T const beta,
                           ptrB B,
                           int ldb)
{
  kernels::acAxpbB(m, n, alpha, to_address(A), lda, to_address(x), incx, beta, to_address(B), ldb);
}

template<class T,
         class ptrA,
         class ptrx,
         class ptrB,
         typename = typename std::enable_if_t<(ptrA::memory_type == CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrx::memory_type == CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrB::memory_type == CPU_OUTOFCARS_POINTER_TYPE)>,
         typename = void>
inline static void acAxpbB(int m,
                           int n,
                           T const alpha,
                           ptrA const A,
                           int lda,
                           ptrx const x,
                           int incx,
                           T const beta,
                           ptrB B,
                           int ldb)
{
  using ma::acAxpbB;
  acAxpbB(m, n, alpha, to_address(A), lda, to_address(x), incx, beta, to_address(B), ldb);
}

// adiagApy
template<class T,
         class ptrA,
         class ptrB,
         typename = typename std::enable_if_t<(ptrA::memory_type != CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrB::memory_type != CPU_OUTOFCARS_POINTER_TYPE)>>
inline static void adiagApy(int n, T const alpha, ptrA const A, int lda, ptrB y, int incy)
{
  kernels::adiagApy(n, alpha, to_address(A), lda, to_address(y), incy);
}

template<class T,
         class ptrA,
         class ptrB,
         typename = typename std::enable_if_t<(ptrA::memory_type == CPU_OUTOFCARS_POINTER_TYPE) or
                                              (ptrB::memory_type == CPU_OUTOFCARS_POINTER_TYPE)>,
         typename = void>
inline static void adiagApy(int n, T const alpha, ptrA const A, int lda, ptrB y, int incy)
{
  using ma::adiagApy;
  adiagApy(n, alpha, to_address(A), lda, to_address(y), incy);
}

template<class ptr, typename = typename std::enable_if_t<(ptr::memory_type != CPU_OUTOFCARS_POINTER_TYPE)>>
inline static auto sum(int n, ptr const x, int incx)
{
  return kernels::sum(n, to_address(x), incx);
}

template<class ptr, typename = typename std::enable_if_t<(ptr::memory_type != CPU_OUTOFCARS_POINTER_TYPE)>>
inline static auto sum(int m, int n, ptr const A, int lda)
{
  return kernels::sum(m, n, to_address(A), lda);
}

template<class ptr,
         typename = typename std::enable_if_t<(ptr::memory_type == CPU_OUTOFCARS_POINTER_TYPE)>,
         typename = void>
inline static auto sum(int n, ptr const x, int incx)
{
  using ma::sum;
  return sum(n, to_address(x), incx);
}

template<class ptr,
         typename = typename std::enable_if_t<(ptr::memory_type == CPU_OUTOFCARS_POINTER_TYPE)>,
         typename = void>
inline static auto sum(int m, int n, ptr const A, int lda)
{
  using ma::sum;
  return sum(m, n, to_address(A), lda);
}

template<class T,
         class ptrA,
         class ptrB,
         class ptrC,
         typename = typename std::enable_if_t<(ptrA::memory_type != CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrB::memory_type != CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrC::memory_type != CPU_OUTOFCARS_POINTER_TYPE)>>
inline static void gemmStridedBatched(char Atrans,
                                      char Btrans,
                                      int M,
                                      int N,
                                      int K,
                                      T const alpha,
                                      ptrA const A,
                                      int lda,
                                      int strideA,
                                      ptrB const B,
                                      int ldb,
                                      int strideB,
                                      T beta,
                                      ptrC C,
                                      int ldc,
                                      int strideC,
                                      int batchSize)
{
  cublas::cublas_gemmStridedBatched(*A.handles.cublas_handle, Atrans, Btrans, M, N, K, alpha, to_address(A), lda,
                                    strideA, to_address(B), ldb, strideB, beta, to_address(C), ldc, strideC, batchSize);
}

template<class T,
         class ptrA,
         class ptrB,
         class ptrC,
         typename = typename std::enable_if_t<(ptrA::memory_type == CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrB::memory_type == CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrC::memory_type == CPU_OUTOFCARS_POINTER_TYPE)>,
         typename = void>
inline static void gemmStridedBatched(char Atrans,
                                      char Btrans,
                                      int M,
                                      int N,
                                      int K,
                                      T const alpha,
                                      ptrA const A,
                                      int lda,
                                      int strideA,
                                      ptrB const B,
                                      int ldb,
                                      int strideB,
                                      T beta,
                                      ptrC C,
                                      int ldc,
                                      int strideC,
                                      int batchSize)
{
  using ma::gemmStridedBatched;
  gemmStridedBatched(Atrans, Btrans, M, N, K, alpha, to_address(A), lda, strideA, to_address(B), ldb, strideB, beta,
                     to_address(C), ldc, strideC, batchSize);
}

template<class T,
         class ptrA,
         class ptrB,
         class ptrC,
         typename = typename std::enable_if_t<(ptrA::memory_type != CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrB::memory_type != CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrC::memory_type != CPU_OUTOFCARS_POINTER_TYPE)>>
inline static void gemmBatched(char Atrans,
                               char Btrans,
                               int M,
                               int N,
                               int K,
                               T const alpha,
                               ptrA const* A,
                               int lda,
                               ptrB const* B,
                               int ldb,
                               T beta,
                               ptrC* C,
                               int ldc,
                               int batchSize)
{
  using Q = typename ptrA::value_type;
  Q **A_d, **B_d, **C_d;
  Q **A_h, **B_h, **C_h;
  A_h = new Q*[batchSize];
  B_h = new Q*[batchSize];
  C_h = new Q*[batchSize];
  for (int i = 0; i < batchSize; i++)
  {
    A_h[i] = to_address(A[i]);
    B_h[i] = to_address(B[i]);
    C_h[i] = to_address(C[i]);
  }
  arch::malloc((void**)&A_d, batchSize * sizeof(*A_h));
  arch::malloc((void**)&B_d, batchSize * sizeof(*B_h));
  arch::malloc((void**)&C_d, batchSize * sizeof(*C_h));
  arch::memcopy(A_d, A_h, batchSize * sizeof(*A_h), arch::memcopyH2D);
  arch::memcopy(B_d, B_h, batchSize * sizeof(*B_h), arch::memcopyH2D);
  arch::memcopy(C_d, C_h, batchSize * sizeof(*C_h), arch::memcopyH2D);
  cublas::cublas_gemmBatched(*(A[0]).handles.cublas_handle, Atrans, Btrans, M, N, K, alpha, A_d, lda, B_d, ldb, beta,
                             C_d, ldc, batchSize);
  arch::free(A_d);
  arch::free(B_d);
  arch::free(C_d);
  delete[] A_h;
  delete[] B_h;
  delete[] C_h;
}

template<class T,
         class ptrA,
         class ptrB,
         class ptrC,
         typename = typename std::enable_if_t<(ptrA::memory_type == CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrB::memory_type == CPU_OUTOFCARS_POINTER_TYPE) and
                                              (ptrC::memory_type == CPU_OUTOFCARS_POINTER_TYPE)>,
         typename = void>
inline static void gemmBatched(char Atrans,
                               char Btrans,
                               int M,
                               int N,
                               int K,
                               T const alpha,
                               ptrA const* A,
                               int lda,
                               ptrB const* B,
                               int ldb,
                               T beta,
                               ptrC* C,
                               int ldc,
                               int batchSize)
{
  using Q = typename ptrA::value_type;
  Q** A_d;
  Q** B_d;
  Q** C_d;
  A_d = new Q*[batchSize];
  B_d = new Q*[batchSize];
  C_d = new Q*[batchSize];
  for (int i = 0; i < batchSize; i++)
  {
    A_d[i] = to_address(A[i]);
    B_d[i] = to_address(B[i]);
    C_d[i] = to_address(C[i]);
  }
  using ma::gemmBatched;
  gemmBatched(Atrans, Btrans, M, N, K, alpha, A_d, lda, B_d, ldb, beta, C_d, ldc, batchSize);
  delete[] A_d;
  delete[] B_d;
  delete[] C_d;
}

} // namespace qmc_cuda

#endif
