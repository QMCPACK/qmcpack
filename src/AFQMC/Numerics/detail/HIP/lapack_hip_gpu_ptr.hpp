///////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Fionn Malone, malone14@llnl.gov, LLNL
//
// File created by: Fionn Malone, malone14@llnl.gov, LLNL
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_LAPACK_GPU_HPP
#define AFQMC_LAPACK_GPU_HPP

#include <cassert>
#include "AFQMC/Utilities/type_conversion.hpp"
#include "AFQMC/Memory/custom_pointers.hpp"
#include "AFQMC/Numerics/detail/HIP/hipblas_wrapper.hpp"
#include "AFQMC/Numerics/detail/HIP/rocsolver_wrapper.hpp"
#include "AFQMC/Numerics/detail/HIP/Kernels/setIdentity.hip.h"

namespace device
{
using qmc_hip::rocsolverHandle_t;
using qmc_hip::rocsolverStatus_t;
using qmcplusplus::afqmc::remove_complex;

template<typename T, typename R, typename I>
inline static void hevr(char JOBZ,
                        char RANGE,
                        char UPLO,
                        int N,
                        device_pointer<T> A,
                        int LDA,
                        T VL,
                        T VU,
                        int IL,
                        int IU,
                        T ABSTOL,
                        int& M,
                        device_pointer<T> W,
                        device_pointer<T> Z,
                        int LDZ,
                        device_pointer<I> ISUPPZ,
                        device_pointer<T> WORK,
                        int& LWORK,
                        device_pointer<R> RWORK,
                        int& LRWORK,
                        device_pointer<I> IWORK,
                        int& LIWORK,
                        int& INFO)
{
  throw std::runtime_error("Error: hevr not implemented in gpu.");
}

// getrf_bufferSize
template<typename T>
inline static void getrf_bufferSize(const int n, const int m, device_pointer<T> a, int lda, int& lwork)
{
  rocsolver::rocsolver_getrf_bufferSize(*a.handles.rocsolver_handle_, n, m, to_address(a), lda, &lwork);
}

template<typename T, typename R, typename I>
inline static void getrf(const int n,
                         const int m,
                         device_pointer<T>&& a,
                         int lda,
                         device_pointer<I>&& piv,
                         int& st,
                         device_pointer<R> work)
{
  rocsolverStatus_t status = rocsolver::rocsolver_getrf(*a.handles.rocsolver_handle_, n, m, to_address(a), lda,
                                                        to_address(work), to_address(piv), to_address(piv) + n);
  arch::memcopy(&st, to_address(piv) + n, sizeof(int), arch::memcopyD2H);
  if (rocblas_status_success != status)
  {
    std::cerr << " hipblas_getrf status, info: " << status << " " << st << std::endl;
    std::cerr.flush();
    throw std::runtime_error("Error: hipblas_getrf returned error code.");
  }
}

// getrfBatched
template<typename T, typename I>
inline static void getrfBatched(const int n,
                                device_pointer<T>* a,
                                int lda,
                                device_pointer<I> piv,
                                device_pointer<I> info,
                                int batchSize)
{
  T** A_d;
  T** A_h;
  A_h = new T*[batchSize];
  for (int i = 0; i < batchSize; i++)
    A_h[i] = to_address(a[i]);
  arch::malloc((void**)&A_d, batchSize * sizeof(*A_h));
  arch::memcopy(A_d, A_h, batchSize * sizeof(*A_h), arch::memcopyH2D);
  hipblasStatus_t status = hipblas::hipblas_getrfBatched(*(a[0]).handles.hipblas_handle, n, A_d, lda, to_address(piv),
                                                         to_address(info), batchSize);
  if (HIPBLAS_STATUS_SUCCESS != status)
    throw std::runtime_error("Error: hipblas_getrf returned error code.");
  arch::free(A_d);
  delete[] A_h;
}

// getri_bufferSize
template<typename T>
inline static void getri_bufferSize(int n, device_pointer<T> a, int lda, int& lwork)
{
  // gpu uses getrs to invert matrix, which requires n*n workspace
  lwork = n * n;
}

// write separate query function to avoid hack!!!
template<typename T, typename R, typename I>
inline static void getri(int n,
                         device_pointer<T> a,
                         int lda,
                         device_pointer<I> piv,
                         device_pointer<R> work,
                         int n1,
                         int& status)
{
  if (n1 < n * n)
    throw std::runtime_error("Error: getri<GPU_MEMORY_POINTER_TYPE> required buffer space of n*n.");
  if (lda != n)
    throw std::runtime_error("Error: getri<GPU_MEMORY_POINTER_TYPE> required lda = n.");

  // info isn't returned from ?getrs.
  int* info;
  arch::malloc((void**)&info, sizeof(int), "lapack_hip_gpu_ptr::getri");

  kernels::set_identity(n, n, to_address(work), n);
  //std::cout << (*work).num_elements() << " " << (*a).num_elements() << std::endl;
  if (rocblas_status_success !=
      rocsolver::rocsolver_getrs(*a.handles.rocsolver_handle_, rocblas_operation_none, n, n, to_address(a), lda,
                                 to_address(piv), to_address(work), n, info))
    throw std::runtime_error("Error: rocsolver_getrs returned error code.");
  arch::memcopy(to_address(a), to_address(work), n * n * sizeof(T), arch::memcopyD2D);
  //arch::memcopy(&status,info,sizeof(int),arch::memcopyD2H);
  status = 0;
  arch::free(info);
}

// getriBatched
template<typename T, typename I>
inline static void getriBatched(int n,
                                device_pointer<T>* a,
                                int lda,
                                device_pointer<I> piv,
                                device_pointer<T>* ainv,
                                int ldc,
                                device_pointer<I> info,
                                int batchSize)
{
  T **A_d, **C_d;
  T **A_h, **C_h;
  A_h = new T*[batchSize];
  C_h = new T*[batchSize];
  for (int i = 0; i < batchSize; i++)
  {
    kernels::set_identity(n, n, to_address(ainv[i]), ldc);
    A_h[i] = to_address(a[i]);
    C_h[i] = to_address(ainv[i]);
  }
  arch::malloc((void**)&A_d, batchSize * sizeof(*A_h));
  arch::malloc((void**)&C_d, batchSize * sizeof(*C_h));
  arch::memcopy(A_d, A_h, batchSize * sizeof(*A_h), arch::memcopyH2D);
  arch::memcopy(C_d, C_h, batchSize * sizeof(*C_h), arch::memcopyH2D);
  hipblasStatus_t status = hipblas::hipblas_getriBatched(*(a[0]).handles.hipblas_handle, HIPBLAS_OP_N, n, n, A_d, lda,
                                                         to_address(piv), C_d, ldc, to_address(info), batchSize);
  if (HIPBLAS_STATUS_SUCCESS != status)
    throw std::runtime_error("Error: hipblas_getri returned error code.");
  arch::free(A_d);
  arch::free(C_d);
  delete[] A_h;
  delete[] C_h;
}

// matinveBatched
template<typename T1, typename T2, typename I>
inline static void matinvBatched(int n,
                                 device_pointer<T1>* a,
                                 int lda,
                                 device_pointer<T2>* ainv,
                                 int lda_inv,
                                 device_pointer<I> info,
                                 int batchSize)
{
  T1 **A_d, **A_h;
  T2 **C_d, **C_h;
  A_h = new T1*[batchSize];
  C_h = new T2*[batchSize];
  for (int i = 0; i < batchSize; i++)
  {
    A_h[i] = to_address(a[i]);
    C_h[i] = to_address(ainv[i]);
  }
  arch::malloc((void**)&A_d, batchSize * sizeof(*A_h));
  arch::malloc((void**)&C_d, batchSize * sizeof(*C_h));
  arch::memcopy(A_d, A_h, batchSize * sizeof(*A_h), arch::memcopyH2D);
  arch::memcopy(C_d, C_h, batchSize * sizeof(*C_h), arch::memcopyH2D);
  hipblasStatus_t status = hipblas::hipblas_matinvBatched(*(a[0]).handles.hipblas_handle, n, A_d, lda, C_d, lda_inv,
                                                          to_address(info), batchSize);
  if (HIPBLAS_STATUS_SUCCESS != status)
    throw std::runtime_error("Error: hipblas_matinv returned error code.");
  arch::free(A_d);
  arch::free(C_d);
  delete[] A_h;
  delete[] C_h;
}

// geqrf
template<typename T>
inline static void geqrf_bufferSize(int m, int n, device_pointer<T> a, int lda, int& lwork)
{
  if (rocblas_status_success !=
      rocsolver::rocsolver_geqrf_bufferSize(*a.handles.rocsolver_handle_, m, n, to_address(a), lda, &lwork))
    throw std::runtime_error("Error: rocsolver_geqrf_bufferSize returned error code.");
}

template<typename T>
inline static void geqrf(int M,
                         int N,
                         device_pointer<T> A,
                         const int LDA,
                         device_pointer<T> TAU,
                         device_pointer<T> WORK,
                         int LWORK,
                         int& INFO)
{
  // allocating here for now
  //size_t dim = std::min(M,N);
  //std::vector<T> piv(dim,0.0);
  //T* dpiv;
  //if(hipSuccess != arch::malloc(&dpiv,sizeof(T)*dim)) {
  //std::cerr << " Error gqr: Error allocating piv on GPU." << std::endl;
  //throw std::runtime_error("Error: arch::malloc returned error code.");
  //}
  rocsolverStatus_t status =
      rocsolver::rocsolver_geqrf(*A.handles.rocsolver_handle_, M, N, to_address(A), LDA, to_address(TAU));
  //arch::memcopy(piv.data(),dpiv,sizeof(T)*dim,arch::memcopyD2H);
  if (rocblas_status_success != status)
  {
    std::cerr << " hipblas_geqrf status, info: " << status << std::endl;
    std::cerr.flush();
    throw std::runtime_error("Error: hipblas_geqrf returned error code.");
  }
  //arch::free(dpiv);
  INFO = 0;
}

// gelqf
template<typename T>
inline static void gelqf_bufferSize(int m, int n, device_pointer<T> a, int lda, int& lwork)
{
  lwork = 0;
}

template<typename T>
inline static void gelqf(int M,
                         int N,
                         device_pointer<T> A,
                         const int LDA,
                         device_pointer<T> TAU,
                         device_pointer<T> WORK,
                         int LWORK,
                         int& INFO)
{
  throw std::runtime_error("Error: gelqf not implemented in HIP backend. \n");
}

// gqr
template<typename T>
static void gqr_bufferSize(int m, int n, int k, device_pointer<T> a, int lda, int& lwork)
{
  if (rocblas_status_success !=
      rocsolver::rocsolver_gqr_bufferSize(*a.handles.rocsolver_handle_, m, n, k, to_address(a), lda, &lwork))
    throw std::runtime_error("Error: rocsolver_gqr_bufferSize returned error code.");
}

template<typename T>
void static gqr(int M,
                int N,
                int K,
                device_pointer<T> A,
                const int LDA,
                device_pointer<T> TAU,
                device_pointer<T> WORK,
                int LWORK,
                int& INFO)
{
  // allocating here for now
  rocsolverStatus_t status =
      rocsolver::rocsolver_gqr(*A.handles.rocsolver_handle_, M, N, K, to_address(A), LDA, to_address(TAU));
  if (rocblas_status_success != status)
  {
    int st;
    std::cerr << " hipblas_gqr status, info: " << status << std::endl;
    std::cerr.flush();
    throw std::runtime_error("Error: hipblas_gqr returned error code.");
  }
  // Not returned from rocm
  INFO = 0;
  //arch::free(piv);
}

template<typename T, typename I>
void static gqrStrided(int M,
                       int N,
                       int K,
                       device_pointer<T> A,
                       const int LDA,
                       const int Astride,
                       device_pointer<T> TAU,
                       const int Tstride,
                       device_pointer<T> WORK,
                       int LWORK,
                       device_pointer<I> info,
                       int batchSize)
{
  rocsolverStatus_t status =
      rocsolver::rocsolver_gqr_strided(*A.handles.rocsolver_handle_, M, N, K, to_address(A), LDA, Astride,
                                       to_address(TAU), Tstride, to_address(WORK), LWORK, to_address(info), batchSize);
  if (rocblas_status_success != status)
  {
    std::cerr << " hipblas_gqr_strided status: " << status << std::endl;
    std::cerr.flush();
    throw std::runtime_error("Error: hipblas_gqr_strided returned error code.");
  }
}

// glq
template<typename T>
static void glq_bufferSize(int m, int n, int k, device_pointer<T> a, int lda, int& lwork)
{
  lwork = 0;
}

template<typename T>
void static glq(int M,
                int N,
                int K,
                device_pointer<T> A,
                const int LDA,
                device_pointer<T> TAU,
                device_pointer<T> WORK,
                int LWORK,
                int& INFO)
{
  throw std::runtime_error("Error: glq not implemented in HIP backend. \n");
}

// batched geqrf
//template<typename T, typename I>
//inline static void geqrfBatched(int M, int N, device_pointer<T> *A, const int LDA, device_pointer<T>* TAU, device_pointer<I> info, int batchSize)
//{
//T **B_h = new T*[2*batchSize];
//T **A_h(B_h);
//T **T_h(B_h+batchSize);
//for(int i=0; i<batchSize; i++)
//A_h[i] = to_address(A[i]);
//for(int i=0; i<batchSize; i++)
//T_h[i] = to_address(TAU[i]);
//T **B_d;
//std::vector<int> inf(batchSize);
//arch::malloc((void **)&B_d,  2*batchSize*sizeof(*B_h));
//arch::memcopy(B_d, B_h, 2*batchSize*sizeof(*B_h), arch::memcopyH2D);
//T **A_d(B_d);
//T **T_d(B_d+batchSize);
//int pstride = std::min(M,N);
//rocsolverStatus_t status = rocsolver::rocsolver_geqrf_batched(
//*(A[0]).handles.rocsolver_handle_, M, N, A_d, LDA, T_d,
//pstride, batchSize);
//if(rocblas_status_success != status)
//throw std::runtime_error("Error: hipblas_geqrfBatched returned error code.");
//arch::free(B_d);
//delete [] B_h;
//}

template<typename T, typename I>
inline static void geqrfStrided(int M,
                                int N,
                                device_pointer<T> A,
                                const int LDA,
                                const int Astride,
                                device_pointer<T> TAU,
                                const int Tstride,
                                device_pointer<I> info,
                                int batchSize)
{
  T** A_h = new T*[batchSize];
  for (int i = 0; i < batchSize; i++)
    A_h[i] = to_address(A) + i * Astride;
  T** A_d;
  arch::malloc((void**)&A_d, batchSize * sizeof(*A_h));
  arch::memcopy(A_d, A_h, batchSize * sizeof(*A_h), arch::memcopyH2D);
  rocsolverStatus_t status = rocsolver::rocsolver_geqrf_batched(*A.handles.rocsolver_handle_, M, N, A_d, LDA,
                                                                to_address(TAU), Tstride, batchSize);
  if (rocblas_status_success != status)
    throw std::runtime_error("Error: hipblas_geqrfStrided returned error code.");
  arch::free(A_d);
  delete[] A_h;
}

// gesvd_bufferSize
template<typename T>
inline static void gesvd_bufferSize(const int m, const int n, device_pointer<T> a, int& lwork)
{
  rocsolver::rocsolver_gesvd_bufferSize(*a.handles.rocsolver_handle_, m, n, to_address(a), &lwork);
}

template<typename T, typename R>
inline static void gesvd(char jobU,
                         char jobVT,
                         const int m,
                         const int n,
                         device_pointer<T>&& A,
                         int lda,
                         device_pointer<R>&& S,
                         device_pointer<T>&& U,
                         int ldu,
                         device_pointer<T>&& VT,
                         int ldvt,
                         device_pointer<T>&& W,
                         int lw,
                         int& st)
{
  int* devSt;
  arch::malloc((void**)&devSt, sizeof(int));
  rocsolverStatus_t status =
      rocsolver::rocsolver_gesvd(*A.handles.rocsolver_handle_, jobU, jobVT, m, n, to_address(A), lda, to_address(S),
                                 to_address(U), ldu, to_address(VT), ldvt, to_address(W), lw, devSt);
  arch::memcopy(&st, devSt, sizeof(int), arch::memcopyD2H);
  if (rocblas_status_success != status)
  {
    std::cerr << " hipblas_gesvd status, info: " << status << " " << st << std::endl;
    std::cerr.flush();
    throw std::runtime_error("Error: hipblas_gesvd returned error code.");
  }
  arch::free(devSt);
}

template<typename T, typename R>
inline static void gesvd(char jobU,
                         char jobVT,
                         const int m,
                         const int n,
                         device_pointer<T>&& A,
                         int lda,
                         device_pointer<R>&& S,
                         device_pointer<T>&& U,
                         int ldu,
                         device_pointer<T>&& VT,
                         int ldvt,
                         device_pointer<T>&& W,
                         int lw,
                         device_pointer<R>&& RW,
                         int& st)
{
  int* devSt;
  arch::malloc((void**)&devSt, sizeof(int));
  rocsolverStatus_t status =
      rocsolver::rocsolver_gesvd(*A.handles.rocsolver_handle_, jobU, jobVT, m, n, to_address(A), lda, to_address(S),
                                 to_address(U), ldu, to_address(VT), ldvt, to_address(W), lw, devSt);
  arch::memcopy(&st, devSt, sizeof(int), arch::memcopyD2H);
  if (rocblas_status_success != status)
  {
    std::cerr << " hipblas_gesvd status, info: " << status << " " << st << std::endl;
    std::cerr.flush();
    throw std::runtime_error("Error: hipblas_gesvd returned error code.");
  }
  arch::free(devSt);
}


} // namespace device

#endif
