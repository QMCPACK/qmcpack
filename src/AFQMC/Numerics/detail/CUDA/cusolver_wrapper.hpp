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

#ifndef CUSOLVER_FUNCTIONDEFS_H
#define CUSOLVER_FUNCTIONDEFS_H

#include <cassert>
#include <cuda_runtime.h>
#include "cusolverDn.h"
#include "AFQMC/Memory/CUDA/cuda_utilities.h"

namespace cusolver
{
using qmc_cuda::cublasOperation;

inline cusolverStatus_t cusolver_getrf_bufferSize(cusolverDnHandle_t handle,
                                                  int m,
                                                  int n,
                                                  float* A,
                                                  int lda,
                                                  int* Lwork)
{
  cusolverStatus_t success = cusolverDnSgetrf_bufferSize(handle, m, n, A, lda, Lwork);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_getrf_bufferSize(cusolverDnHandle_t handle,
                                                  int m,
                                                  int n,
                                                  std::complex<float>* A,
                                                  int lda,
                                                  int* Lwork)
{
  cusolverStatus_t success = cusolverDnCgetrf_bufferSize(handle, m, n, reinterpret_cast<cuComplex*>(A), lda, Lwork);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_getrf_bufferSize(cusolverDnHandle_t handle,
                                                  int m,
                                                  int n,
                                                  double* A,
                                                  int lda,
                                                  int* Lwork)
{
  cusolverStatus_t success = cusolverDnDgetrf_bufferSize(handle, m, n, A, lda, Lwork);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_getrf_bufferSize(cusolverDnHandle_t handle,
                                                  int m,
                                                  int n,
                                                  std::complex<double>* A,
                                                  int lda,
                                                  int* Lwork)
{
  cusolverStatus_t success =
      cusolverDnZgetrf_bufferSize(handle, m, n, reinterpret_cast<cuDoubleComplex*>(A), lda, Lwork);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_getrf(cusolverDnHandle_t handle,
                                       int m,
                                       int n,
                                       float* A,
                                       int lda,
                                       float* Work,
                                       int* devIpiv,
                                       int* devInfo)
{
  cusolverStatus_t success = cusolverDnSgetrf(handle, m, n, A, lda, Work, devIpiv, devInfo);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_getrf(cusolverDnHandle_t handle,
                                       int m,
                                       int n,
                                       double* A,
                                       int lda,
                                       double* Work,
                                       int* devIpiv,
                                       int* devInfo)
{
  cusolverStatus_t success = cusolverDnDgetrf(handle, m, n, A, lda, Work, devIpiv, devInfo);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_getrf(cusolverDnHandle_t handle,
                                       int m,
                                       int n,
                                       std::complex<float>* A,
                                       int lda,
                                       std::complex<float>* Work,
                                       int* devIpiv,
                                       int* devInfo)
{
  cusolverStatus_t success = cusolverDnCgetrf(handle, m, n, reinterpret_cast<cuComplex*>(A), lda,
                                             reinterpret_cast<cuComplex*>(Work), devIpiv, devInfo);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_getrf(cusolverDnHandle_t handle,
                                       int m,
                                       int n,
                                       std::complex<double>* A,
                                       int lda,
                                       std::complex<double>* Work,
                                       int* devIpiv,
                                       int* devInfo)
{
  cusolverStatus_t success = cusolverDnZgetrf(handle, m, n, reinterpret_cast<cuDoubleComplex*>(A), lda,
                                             reinterpret_cast<cuDoubleComplex*>(Work), devIpiv, devInfo);
  cudaDeviceSynchronize();
  return success;
}


// getrs
inline cusolverStatus_t cusolver_getrs(cusolverDnHandle_t handle,
                                       cublasOperation_t trans,
                                       int n,
                                       int nrhs,
                                       const float* A,
                                       int lda,
                                       const int* devIpiv,
                                       float* B,
                                       int ldb,
                                       int* devInfo)
{
  cusolverStatus_t success = cusolverDnSgetrs(handle, trans, n, nrhs, A, lda, devIpiv, B, ldb, devInfo);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_getrs(cusolverDnHandle_t handle,
                                       cublasOperation_t trans,
                                       int n,
                                       int nrhs,
                                       const double* A,
                                       int lda,
                                       const int* devIpiv,
                                       double* B,
                                       int ldb,
                                       int* devInfo)
{
  cusolverStatus_t success = cusolverDnDgetrs(handle, trans, n, nrhs, A, lda, devIpiv, B, ldb, devInfo);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_getrs(cusolverDnHandle_t handle,
                                       cublasOperation_t trans,
                                       int n,
                                       int nrhs,
                                       const std::complex<float>* A,
                                       int lda,
                                       const int* devIpiv,
                                       std::complex<float>* B,
                                       int ldb,
                                       int* devInfo)
{
  cusolverStatus_t success = cusolverDnCgetrs(handle, trans, n, nrhs, reinterpret_cast<cuComplex const*>(A), lda,
                                             devIpiv, reinterpret_cast<cuComplex*>(B), ldb, devInfo);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_getrs(cusolverDnHandle_t handle,
                                       cublasOperation_t trans,
                                       int n,
                                       int nrhs,
                                       const std::complex<double>* A,
                                       int lda,
                                       const int* devIpiv,
                                       std::complex<double>* B,
                                       int ldb,
                                       int* devInfo)
{
  cusolverStatus_t success = cusolverDnZgetrs(handle, trans, n, nrhs, reinterpret_cast<cuDoubleComplex const*>(A), lda,
                                             devIpiv, reinterpret_cast<cuDoubleComplex*>(B), ldb, devInfo);
  cudaDeviceSynchronize();
  return success;
}

//geqrf_bufferSize
inline cusolverStatus_t cusolver_geqrf_bufferSize(cusolverDnHandle_t handle,
                                                  int m,
                                                  int n,
                                                  float* A,
                                                  int lda,
                                                  int* Lwork)
{
  cusolverStatus_t success = cusolverDnSgeqrf_bufferSize(handle, m, n, A, lda, Lwork);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_geqrf_bufferSize(cusolverDnHandle_t handle,
                                                  int m,
                                                  int n,
                                                  double* A,
                                                  int lda,
                                                  int* Lwork)
{
  cusolverStatus_t success = cusolverDnDgeqrf_bufferSize(handle, m, n, A, lda, Lwork);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_geqrf_bufferSize(cusolverDnHandle_t handle,
                                                  int m,
                                                  int n,
                                                  std::complex<float>* A,
                                                  int lda,
                                                  int* Lwork)
{
  cusolverStatus_t success = cusolverDnCgeqrf_bufferSize(handle, m, n, reinterpret_cast<cuComplex*>(A), lda, Lwork);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_geqrf_bufferSize(cusolverDnHandle_t handle,
                                                  int m,
                                                  int n,
                                                  std::complex<double>* A,
                                                  int lda,
                                                  int* Lwork)
{
  cusolverStatus_t success =
      cusolverDnZgeqrf_bufferSize(handle, m, n, reinterpret_cast<cuDoubleComplex*>(A), lda, Lwork);
  cudaDeviceSynchronize();
  return success;
}

//geqrf
inline cusolverStatus_t cusolver_geqrf(cusolverDnHandle_t handle,
                                       int m,
                                       int n,
                                       float* A,
                                       int lda,
                                       float* TAU,
                                       float* Workspace,
                                       int Lwork,
                                       int* devInfo)
{
  cusolverStatus_t success = cusolverDnSgeqrf(handle, m, n, A, lda, TAU, Workspace, Lwork, devInfo);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_geqrf(cusolverDnHandle_t handle,
                                       int m,
                                       int n,
                                       double* A,
                                       int lda,
                                       double* TAU,
                                       double* Workspace,
                                       int Lwork,
                                       int* devInfo)
{
  cusolverStatus_t success = cusolverDnDgeqrf(handle, m, n, A, lda, TAU, Workspace, Lwork, devInfo);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_geqrf(cusolverDnHandle_t handle,
                                       int m,
                                       int n,
                                       std::complex<float>* A,
                                       int lda,
                                       std::complex<float>* TAU,
                                       std::complex<float>* Workspace,
                                       int Lwork,
                                       int* devInfo)
{
  cusolverStatus_t success =
      cusolverDnCgeqrf(handle, m, n, reinterpret_cast<cuComplex*>(A), lda, reinterpret_cast<cuComplex*>(TAU),
                       reinterpret_cast<cuComplex*>(Workspace), Lwork, devInfo);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_geqrf(cusolverDnHandle_t handle,
                                       int m,
                                       int n,
                                       std::complex<double>* A,
                                       int lda,
                                       std::complex<double>* TAU,
                                       std::complex<double>* Workspace,
                                       int Lwork,
                                       int* devInfo)
{
  cusolverStatus_t success = cusolverDnZgeqrf(handle, m, n, reinterpret_cast<cuDoubleComplex*>(A), lda,
                                             reinterpret_cast<cuDoubleComplex*>(TAU),
                                             reinterpret_cast<cuDoubleComplex*>(Workspace), Lwork, devInfo);
  cudaDeviceSynchronize();
  return success;
}


//gqr_bufferSize
// MAM: There is an inconsistency between the documentation online and the header files
// in cuda.9.2.148 in the declaration of: cusolverDnXungqr_bufferSize and cusolverDnXorgqr_bufferSize
//  I'm going to hack it right here, under the assumption that *tau is not used at all.
inline cusolverStatus_t cusolver_gqr_bufferSize(cusolverDnHandle_t handle,
                                                int m,
                                                int n,
                                                int k,
                                                const float* A,
                                                int lda,
                                                int* lwork)
{
  cusolverStatus_t success = cusolverDnSorgqr_bufferSize(handle, m, n, k, A, lda, A, lwork);
  // HACK
  //              cusolverDnSorgqr_bufferSize(handle,m,n,k,A,lda,lwork);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_gqr_bufferSize(cusolverDnHandle_t handle,
                                                int m,
                                                int n,
                                                int k,
                                                const double* A,
                                                int lda,
                                                int* lwork)
{
  cusolverStatus_t success = cusolverDnDorgqr_bufferSize(handle, m, n, k, A, lda, A, lwork);
  // HACK
  //              cusolverDnDorgqr_bufferSize(handle,m,n,k,A,lda,lwork);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_gqr_bufferSize(cusolverDnHandle_t handle,
                                                int m,
                                                int n,
                                                int k,
                                                const std::complex<float>* A,
                                                int lda,
                                                int* lwork)
{
  cusolverStatus_t success = cusolverDnCungqr_bufferSize(handle, m, n, k, reinterpret_cast<cuComplex const*>(A), lda,
                                                        reinterpret_cast<cuComplex const*>(A), lwork);
  // HACK
  //                                          lwork);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_gqr_bufferSize(cusolverDnHandle_t handle,
                                                int m,
                                                int n,
                                                int k,
                                                const std::complex<double>* A,
                                                int lda,
                                                int* lwork)
{
  cusolverStatus_t success = cusolverDnZungqr_bufferSize(handle, m, n, k, reinterpret_cast<cuDoubleComplex const*>(A),
                                                        lda, reinterpret_cast<cuDoubleComplex const*>(A), lwork);
  // HACK
  //                                          lwork);
  cudaDeviceSynchronize();
  return success;
}

//gqr
inline cusolverStatus_t cusolver_gqr(cusolverDnHandle_t handle,
                                     int m,
                                     int n,
                                     int k,
                                     float* A,
                                     int lda,
                                     const float* tau,
                                     float* work,
                                     int lwork,
                                     int* devInfo)
{
  cusolverStatus_t success = cusolverDnSorgqr(handle, m, n, k, A, lda, tau, work, lwork, devInfo);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_gqr(cusolverDnHandle_t handle,
                                     int m,
                                     int n,
                                     int k,
                                     double* A,
                                     int lda,
                                     const double* tau,
                                     double* work,
                                     int lwork,
                                     int* devInfo)
{
  cusolverStatus_t success = cusolverDnDorgqr(handle, m, n, k, A, lda, tau, work, lwork, devInfo);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_gqr(cusolverDnHandle_t handle,
                                     int m,
                                     int n,
                                     int k,
                                     std::complex<float>* A,
                                     int lda,
                                     const std::complex<float>* tau,
                                     std::complex<float>* work,
                                     int lwork,
                                     int* devInfo)
{
  cusolverStatus_t success =
      cusolverDnCungqr(handle, m, n, k, reinterpret_cast<cuComplex*>(A), lda, reinterpret_cast<cuComplex const*>(tau),
                       reinterpret_cast<cuComplex*>(work), lwork, devInfo);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_gqr(cusolverDnHandle_t handle,
                                     int m,
                                     int n,
                                     int k,
                                     std::complex<double>* A,
                                     int lda,
                                     const std::complex<double>* tau,
                                     std::complex<double>* work,
                                     int lwork,
                                     int* devInfo)
{
  cusolverStatus_t success = cusolverDnZungqr(handle, m, n, k, reinterpret_cast<cuDoubleComplex*>(A), lda,
                                             reinterpret_cast<cuDoubleComplex const*>(tau),
                                             reinterpret_cast<cuDoubleComplex*>(work), lwork, devInfo);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_gqr_strided(cusolverDnHandle_t handle,
                                             const int m,
                                             const int n,
                                             const int k,
                                             std::complex<double>* A,
                                             const int lda,
                                             const int Astride,
                                             const std::complex<double>* tau,
                                             const int tstride,
                                             std::complex<double>* work,
                                             int lwork,
                                             int* devInfo,
                                             int batchsize)
{
  using qmc_cuda::afqmc_cuda_streams;
  if (afqmc_cuda_streams.size() < batchsize)
  {
    int n0 = afqmc_cuda_streams.size();
    for (int n = n0; n < batchsize; n++)
    {
      afqmc_cuda_streams.emplace_back(cudaStream_t{});
      cudaStreamCreate(&(afqmc_cuda_streams.back()));
    }
  }

  cudaStream_t s0;
  qmc_cuda::cusolver_check(cusolverDnGetStream(handle, &s0), "cusolverDnGetStream");

  for (int i = 0; i < batchsize; i++)
  {
    qmc_cuda::cusolver_check(cusolverDnSetStream(handle, afqmc_cuda_streams[i]), "cusolverDnSetStream");
    cusolverStatus_t success =
        cusolverDnZungqr(handle, m, n, k, reinterpret_cast<cuDoubleComplex*>(A) + i * Astride, lda,
                         reinterpret_cast<cuDoubleComplex const*>(tau) + i * tstride,
                         reinterpret_cast<cuDoubleComplex*>(work) + i * lwork, lwork, devInfo + i);
    qmc_cuda::cusolver_check(success, "cusolver_gqr_strided_status");
  }
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "cusolver_gqr_strided_sync");
  qmc_cuda::cuda_check(cudaGetLastError(), "cusolver_gqr_strided_error");
  qmc_cuda::cusolver_check(cusolverDnSetStream(handle, s0), "cusolverDnSetStream");

  return CUSOLVER_STATUS_SUCCESS;
}

//gesvd_bufferSize
inline cusolverStatus_t cusolver_gesvd_bufferSize(cusolverDnHandle_t handle, int m, int n, float* A, int* Lwork)
{
  cusolverStatus_t success = cusolverDnSgesvd_bufferSize(handle, m, n, Lwork);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_gesvd_bufferSize(cusolverDnHandle_t handle, int m, int n, double* A, int* Lwork)
{
  cusolverStatus_t success = cusolverDnDgesvd_bufferSize(handle, m, n, Lwork);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_gesvd_bufferSize(cusolverDnHandle_t handle,
                                                  int m,
                                                  int n,
                                                  std::complex<float>* A,
                                                  int* Lwork)
{
  cusolverStatus_t success = cusolverDnCgesvd_bufferSize(handle, m, n, Lwork);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_gesvd_bufferSize(cusolverDnHandle_t handle,
                                                  int m,
                                                  int n,
                                                  std::complex<double>* A,
                                                  int* Lwork)
{
  cusolverStatus_t success = cusolverDnZgesvd_bufferSize(handle, m, n, Lwork);
  cudaDeviceSynchronize();
  return success;
}

//gesvd
inline cusolverStatus_t cusolver_gesvd(cusolverDnHandle_t handle,
                                       signed char jobu,
                                       signed char jobvt,
                                       int m,
                                       int n,
                                       float* A,
                                       int lda,
                                       float* S,
                                       float* U,
                                       int ldu,
                                       float* VT,
                                       int ldvt,
                                       float* work,
                                       int lwork,
                                       int* devInfo)
{
  cusolverStatus_t success =
      cusolverDnSgesvd(handle, jobu, jobvt, m, n, A, lda, S, U, ldu, VT, ldvt, work, lwork, nullptr, devInfo);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_gesvd(cusolverDnHandle_t handle,
                                       signed char jobu,
                                       signed char jobvt,
                                       int m,
                                       int n,
                                       double* A,
                                       int lda,
                                       double* S,
                                       double* U,
                                       int ldu,
                                       double* VT,
                                       int ldvt,
                                       double* work,
                                       int lwork,
                                       int* devInfo)
{
  cusolverStatus_t success =
      cusolverDnDgesvd(handle, jobu, jobvt, m, n, A, lda, S, U, ldu, VT, ldvt, work, lwork, nullptr, devInfo);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_gesvd(cusolverDnHandle_t handle,
                                       signed char jobu,
                                       signed char jobvt,
                                       int m,
                                       int n,
                                       std::complex<float>* A,
                                       int lda,
                                       float* S,
                                       std::complex<float>* U,
                                       int ldu,
                                       std::complex<float>* VT,
                                       int ldvt,
                                       std::complex<float>* work,
                                       int lwork,
                                       int* devInfo)
{
  cusolverStatus_t success = cusolverDnCgesvd(handle, jobu, jobvt, m, n, reinterpret_cast<cuComplex*>(A), lda, S,
                                             reinterpret_cast<cuComplex*>(U), ldu, reinterpret_cast<cuComplex*>(VT),
                                             ldvt, reinterpret_cast<cuComplex*>(work), lwork, nullptr, devInfo);
  cudaDeviceSynchronize();
  return success;
}

inline cusolverStatus_t cusolver_gesvd(cusolverDnHandle_t handle,
                                       signed char jobu,
                                       signed char jobvt,
                                       int m,
                                       int n,
                                       std::complex<double>* A,
                                       int lda,
                                       double* S,
                                       std::complex<double>* U,
                                       int ldu,
                                       std::complex<double>* VT,
                                       int ldvt,
                                       std::complex<double>* work,
                                       int lwork,
                                       int* devInfo)
{
  cusolverStatus_t success =
      cusolverDnZgesvd(handle, jobu, jobvt, m, n, reinterpret_cast<cuDoubleComplex*>(A), lda, S,
                       reinterpret_cast<cuDoubleComplex*>(U), ldu, reinterpret_cast<cuDoubleComplex*>(VT), ldvt,
                       reinterpret_cast<cuDoubleComplex*>(work), lwork, nullptr, devInfo);
  cudaDeviceSynchronize();
  return success;
}

} // namespace cusolver

#endif
