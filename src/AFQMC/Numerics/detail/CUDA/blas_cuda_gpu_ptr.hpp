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

#ifndef AFQMC_BLAS_CUDA_GPU_PTR_HPP
#define AFQMC_BLAS_CUDA_GPU_PTR_HPP

#include <type_traits>
#include <cassert>
#include <vector>
//#include "AFQMC/Memory/CUDA/cuda_gpu_pointer.hpp"
#include "AFQMC/Utilities/type_conversion.hpp"
#include "AFQMC/Memory/device_pointers.hpp"
#include "AFQMC/Memory/arch.hpp"
#include "AFQMC/Numerics/detail/CUDA/cublas_wrapper.hpp"
//#include "AFQMC/Numerics/detail/CUDA/cublasXt_wrapper.hpp"
// hand coded kernels for blas extensions
#include "AFQMC/Numerics/detail/CUDA/Kernels/adotpby.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/setIdentity.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/axty.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/sum.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/adiagApy.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/acAxpbB.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/zero_complex_part.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/axpyBatched.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/get_diagonal.cuh"

// Currently available:
// Lvl-1: dot, axpy, scal
// Lvl-2: gemv
// Lvl-3: gemm

namespace device
{
// copy Specializations
template<typename T, typename Q>
inline static void copy(int n, device_pointer<Q> x, int incx, device_pointer<T> y, int incy)
{
  static_assert(std::is_same<typename std::decay<Q>::type, T>::value, "Wrong dispatch.\n");
  if (CUBLAS_STATUS_SUCCESS !=
      cublas::cublas_copy(*x.handles.cublas_handle, n, to_address(x), incx, to_address(y), incy))
    throw std::runtime_error("Error: cublas_copy returned error code.");
}

template<typename T, typename Q>
inline static void copy(int n, T const* x, int incx, device_pointer<Q> y, int incy)
{
  static_assert(std::is_same<typename std::decay<Q>::type, T>::value, "Wrong dispatch.\n");
  arch::memcopy2D(to_address(y), sizeof(Q) * incy, x, sizeof(T) * incx, sizeof(T), n, arch::memcopyH2D,
                  "lapack_cuda_gpu_ptr::copy");
}

template<typename T, typename Q>
inline static void copy(int n, device_pointer<Q> x, int incx, T* y, int incy)
{
  static_assert(std::is_same<typename std::decay<Q>::type, T>::value, "Wrong dispatch.\n");
  assert(sizeof(Q) == sizeof(T));
  arch::memcopy2D(y, sizeof(T) * incy, to_address(x), sizeof(Q) * incx, sizeof(T), n, arch::memcopyD2H,
                  "lapack_cuda_gpu_ptr::copy");
}

// scal Specializations
template<typename T, typename Q>
inline static void scal(int n, Q alpha, device_pointer<T> x, int incx = 1)
{
  static_assert(std::is_convertible<typename std::decay<Q>::type, T>::value, "Wrong dispatch.\n");
  if (CUBLAS_STATUS_SUCCESS != cublas::cublas_scal(*x.handles.cublas_handle, n, T(alpha), to_address(x), incx))
    throw std::runtime_error("Error: cublas_scal returned error code.");
}

// dot Specializations
template<typename T, typename Q>
inline static auto dot(int const n, device_pointer<Q> x, int const incx, device_pointer<T> y, int const incy)
{
  static_assert(std::is_same<typename std::decay<Q>::type, typename std::decay<T>::type>::value, "Wrong dispatch.\n");
  return cublas::cublas_dot(*x.handles.cublas_handle, n, to_address(x), incx, to_address(y), incy);
}

// axpy Specializations
template<typename T, typename Q>
inline static void axpy(int n, T const a, device_pointer<Q> x, int incx, device_pointer<T> y, int incy)
{
  static_assert(std::is_same<typename std::decay<Q>::type, T>::value, "Wrong dispatch.\n");
  if (CUBLAS_STATUS_SUCCESS !=
      cublas::cublas_axpy(*x.handles.cublas_handle, n, a, to_address(x), incx, to_address(y), incy))
    throw std::runtime_error("Error: cublas_axpy returned error code.");
}

// GEMV Specializations
template<typename T, typename T2, typename Q1, typename Q2>
inline static void gemv(char Atrans,
                        int M,
                        int N,
                        T2 alpha,
                        device_pointer<Q1> A,
                        int lda,
                        device_pointer<Q2> x,
                        int incx,
                        T2 beta,
                        device_pointer<T> y,
                        int incy)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T2>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  if (CUBLAS_STATUS_SUCCESS !=
      cublas::cublas_gemv(*A.handles.cublas_handle, Atrans, M, N, alpha, to_address(A), lda, to_address(x), incx, beta,
                          to_address(y), incy))
    throw std::runtime_error("Error: cublas_gemv returned error code.");
}

// GEMM Specializations
// why is this not working with T const????
template<typename T, typename T2, typename Q1, typename Q2>
inline static void gemm(char Atrans,
                        char Btrans,
                        int M,
                        int N,
                        int K,
                        T2 alpha,
                        device_pointer<Q1> A,
                        int lda,
                        device_pointer<Q2> B,
                        int ldb,
                        T2 beta,
                        device_pointer<T> C,
                        int ldc)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T2>::value, "Wrong dispatch.\n");
  if (CUBLAS_STATUS_SUCCESS !=
      cublas::cublas_gemm(*A.handles.cublas_handle, Atrans, Btrans, M, N, K, alpha, to_address(A), lda, to_address(B),
                          ldb, beta, to_address(C), ldc))
    throw std::runtime_error("Error: cublas_gemm returned error code.");
}

// Blas Extensions
// geam
template<typename T, typename Q1, typename Q2>
inline static void geam(char Atrans,
                        char Btrans,
                        int M,
                        int N,
                        T const alpha,
                        device_pointer<Q1> A,
                        int lda,
                        T const beta,
                        device_pointer<Q2> B,
                        int ldb,
                        device_pointer<T> C,
                        int ldc)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  if (CUBLAS_STATUS_SUCCESS !=
      cublas::cublas_geam(*A.handles.cublas_handle, Atrans, Btrans, M, N, alpha, to_address(A), lda, beta,
                          to_address(B), ldb, to_address(C), ldc))
    throw std::runtime_error("Error: cublas_geam returned error code.");
}

template<typename T>
//inline static void set1D(int n, T const alpha, ptr x, int incx)
inline static void set1D(int n, T const alpha, device_pointer<T> x, int incx)
{
  // No set funcion in cuda!!! Avoiding kernels for now
  //std::vector<T> buff(n,alpha);
  //if(CUBLAS_STATUS_SUCCESS != cublasSetVector(n,sizeof(T),buff.data(),1,to_address(x),incx))
  T alpha_(alpha);
  if (CUBLAS_STATUS_SUCCESS != cublasSetVector(n, sizeof(T), std::addressof(alpha), 1, to_address(x), incx))
    throw std::runtime_error("Error: cublasSetVector returned error code.");
}

// dot extension
template<typename T, typename T1, typename T2, typename Q1, typename Q2>
inline static void adotpby(int const n,
                           T1 const alpha,
                           device_pointer<Q1> x,
                           int const incx,
                           device_pointer<Q2> y,
                           int const incy,
                           T2 const beta,
                           T* result)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T1>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T1>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<T2>::type, T>::value, "Wrong dispatch.\n");
  kernels::adotpby(n, alpha, to_address(x), incx, to_address(y), incy, beta, result);
}

// dot extension
template<typename T, typename T1, typename T2, typename Q1, typename Q2>
inline static void strided_adotpby(int nk,
                                   int const n,
                                   T1 const alpha,
                                   device_pointer<Q1> A,
                                   int const lda,
                                   device_pointer<Q2> B,
                                   int const ldb,
                                   T2 const beta,
                                   T* y,
                                   int inc)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T1>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T1>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<T2>::type, T>::value, "Wrong dispatch.\n");
  kernels::strided_adotpby(nk, n, alpha, to_address(A), lda, to_address(B), ldb, beta, y, inc);
}

// axty
template<typename T, typename Q>
inline static void axty(int n, T const alpha, device_pointer<Q> x, int incx, device_pointer<T> y, int incy)
{
  static_assert(std::is_same<typename std::decay<Q>::type, T>::value, "Wrong dispatch.\n");
  if (incx != 1 || incy != 1)
    throw std::runtime_error("Error: axty with inc != 1 not implemented.");
  kernels::axty(n, alpha, to_address(x), to_address(y));
}

// acAxpbB
template<typename T, typename Q1, typename Q2>
inline static void acAxpbB(int m,
                           int n,
                           T const alpha,
                           device_pointer<Q1> A,
                           int lda,
                           device_pointer<Q2> x,
                           int incx,
                           T const beta,
                           device_pointer<T> B,
                           int ldb)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  kernels::acAxpbB(m, n, alpha, to_address(A), lda, to_address(x), incx, beta, to_address(B), ldb);
}

// adiagApy
template<typename T, typename Q1>
inline static void adiagApy(int n, T const alpha, device_pointer<Q1> A, int lda, device_pointer<T> y, int incy)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  kernels::adiagApy(n, alpha, to_address(A), lda, to_address(y), incy);
}

template<typename T>
inline static void zero_complex_part(int n, device_pointer<T> x)
{
  kernels::zero_complex_part(n, to_address(x));
}

template<typename T>
inline static auto sum(int n, device_pointer<T> x, int incx)
{
  return kernels::sum(n, to_address(x), incx);
}

template<typename T>
inline static auto sum(int m, int n, device_pointer<T> A, int lda)
{
  return kernels::sum(m, n, to_address(A), lda);
}

template<typename T>
void set_identity(int m, int n, device_pointer<T> A, int lda)
{
  kernels::set_identity(m, n, to_address(A), lda);
}

template<typename T>
void set_identity_strided(int nbatch, int stride, int m, int n, device_pointer<T> A, int lda)
{
  kernels::set_identity_strided(nbatch, stride, m, n, to_address(A), lda);
}

template<typename T, typename Q1, typename Q2>
inline static void gemmStridedBatched(char Atrans,
                                      char Btrans,
                                      int M,
                                      int N,
                                      int K,
                                      T const alpha,
                                      device_pointer<Q1> A,
                                      int lda,
                                      int strideA,
                                      device_pointer<Q2> B,
                                      int ldb,
                                      int strideB,
                                      T beta,
                                      device_pointer<T> C,
                                      int ldc,
                                      int strideC,
                                      int batchSize)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  cublas::cublas_gemmStridedBatched(*A.handles.cublas_handle, Atrans, Btrans, M, N, K, alpha, to_address(A), lda,
                                    strideA, to_address(B), ldb, strideB, beta, to_address(C), ldc, strideC, batchSize);
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
                               T const alpha,
                               device_pointer<Q1>* A,
                               int lda,
                               device_pointer<Q2>* B,
                               int ldb,
                               T const beta,
                               device_pointer<T>* C,
                               int ldc,
                               int batchSize)
{
  static_assert(std::is_same<typename std::decay<Q1>::type, T>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  // replace with single call to arch::malloc and arch::memcopy
  T **A_d, **B_d, **C_d;
  Q1** A_h;
  Q2** B_h;
  T** C_h;
  A_h = new Q1*[batchSize];
  B_h = new Q2*[batchSize];
  C_h = new T*[batchSize];
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
                               T const alpha,
                               device_pointer<Q1>* A,
                               int lda,
                               device_pointer<Q2>* B,
                               int ldb,
                               T const beta,
                               device_pointer<T2>* C,
                               int ldc,
                               int batchSize)
{
  // check that remove_complex<T2> == T ???
  static_assert(std::is_same<typename std::decay<Q1>::type, T2>::value, "Wrong dispatch.\n");
  static_assert(std::is_same<typename std::decay<Q2>::type, T>::value, "Wrong dispatch.\n");
  assert(Atrans == 'N' || Atrans == 'n');
  // replace with single call to arch::malloc and arch::memcopy
  T2** A_d;
  T** B_d;
  T2** C_d;
  Q1** A_h;
  Q2** B_h;
  T2** C_h;
  A_h = new Q1*[batchSize];
  B_h = new Q2*[batchSize];
  C_h = new T2*[batchSize];
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

template<typename T1, typename T2, typename T3>
inline static void axpyBatched(int n,
                               T1* x,
                               device_pointer<T2>* a,
                               int inca,
                               device_pointer<T3>* b,
                               int incb,
                               int batchSize)
{
  T2 const** a_ = new T2 const*[batchSize];
  T3** b_       = new T3*[batchSize];
  for (int i = 0; i < batchSize; i++)
  {
    a_[i] = to_address(a[i]);
    b_[i] = to_address(b[i]);
  }
  kernels::axpy_batched_gpu(n, x, a_, inca, b_, incb, batchSize);
  delete[] a_;
  delete[] b_;
}

template<typename T1, typename T2, typename T3>
inline static void sumGwBatched(int n,
                                T1* x,
                                device_pointer<T2>* a,
                                int inca,
                                device_pointer<T3>* b,
                                int incb,
                                int b0,
                                int nw,
                                int batchSize)
{
  T2 const** a_ = new T2 const*[batchSize];
  T3** b_       = new T3*[batchSize];
  for (int i = 0; i < batchSize; i++)
  {
    a_[i] = to_address(a[i]);
    b_[i] = to_address(b[i]);
  }
  kernels::sumGw_batched_gpu(n, x, a_, inca, b_, incb, b0, nw, batchSize);
  delete[] a_;
  delete[] b_;
}

template<typename T, typename T2>
inline static void copy2D(int N, int M, device_pointer<T> src, int lda, device_pointer<T2> dst, int ldb)
{
  static_assert(std::is_same<typename std::decay<T>::type, T2>::value, "Wrong dispatch.\n");
  arch::memcopy2D(to_address(dst), sizeof(T2) * ldb, to_address(src), sizeof(T) * lda, M * sizeof(T), N,
                  arch::memcopyD2D, "blas_cuda_gpu_ptr::copy2D");
}

template<typename T, typename T2>
inline static void copy2D(int N, int M, T const* src, int lda, device_pointer<T2> dst, int ldb)
{
  static_assert(std::is_same<typename std::decay<T>::type, T2>::value, "Wrong dispatch.\n");
  arch::memcopy2D(to_address(dst), sizeof(T2) * ldb, src, sizeof(T) * lda, M * sizeof(T), N, arch::memcopyH2D,
                  "blas_cuda_gpu_ptr::copy2D");
}

template<typename T, typename T2>
inline static void copy2D(int N, int M, device_pointer<T> src, int lda, T2* dst, int ldb)
{
  static_assert(std::is_same<typename std::decay<T>::type, T2>::value, "Wrong dispatch.\n");
  arch::memcopy2D(dst, sizeof(T2) * ldb, to_address(src), sizeof(T) * lda, M * sizeof(T), N, arch::memcopyD2H,
                  "blas_cuda_gpu_ptr::copy2D");
}

template<typename T, typename T2>
inline static void get_diagonal_strided(int nk,
                                        int ni,
                                        device_pointer<T> A,
                                        int lda,
                                        int stride,
                                        device_pointer<T2> B,
                                        int ldb)
{
  kernels::get_diagonal_strided(nk, ni, to_address(A), lda, stride, to_address(B), ldb);
}
} // namespace device

#endif
