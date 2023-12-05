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

#ifndef HIPSOLVER_FUNCTIONDEFS_H
#define HIPSOLVER_FUNCTIONDEFS_H

#include <cassert>
#include <hip/hip_runtime.h>
#include <rocsolver/rocsolver.h>
#include "AFQMC/Memory/HIP/hip_utilities.h"
#include "AFQMC/Numerics/detail/CPU/lapack_cpu.hpp"

// Abusing hip/rocm names
namespace rocsolver
{
using qmc_hip::hipblasOperation;
using qmc_hip::rocsolverHandle_t;
using qmc_hip::rocsolverStatus_t;

// TODO: FDM These don't exist in rocsolver.
inline rocsolverStatus_t rocsolver_getrf_bufferSize(rocsolverHandle_t handle,
                                                    int m,
                                                    int n,
                                                    float* A,
                                                    int lda,
                                                    int* Lwork)
{
  rocsolverStatus_t success;
  success = rocblas_status_success;
  *Lwork  = 0;
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_getrf_bufferSize(rocsolverHandle_t handle,
                                                    int m,
                                                    int n,
                                                    std::complex<float>* A,
                                                    int lda,
                                                    int* Lwork)
{
  rocsolverStatus_t success;
  success = rocblas_status_success;
  *Lwork  = 0;
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_getrf_bufferSize(rocsolverHandle_t handle,
                                                    int m,
                                                    int n,
                                                    double* A,
                                                    int lda,
                                                    int* Lwork)
{
  rocsolverStatus_t success;
  success = rocblas_status_success;
  *Lwork  = 0;
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_getrf_bufferSize(rocsolverHandle_t handle,
                                                    int m,
                                                    int n,
                                                    std::complex<double>* A,
                                                    int lda,
                                                    int* Lwork)
{
  rocsolverStatus_t success;
  success = rocblas_status_success;
  *Lwork  = 0;
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_getrf(rocsolverHandle_t handle,
                                         int m,
                                         int n,
                                         float* A,
                                         int lda,
                                         float* Work,
                                         int* devIpiv,
                                         int* devInfo)
{
  rocsolverStatus_t success = rocsolver_sgetrf(handle, m, n, A, lda, devIpiv, devInfo);
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_getrf(rocsolverHandle_t handle,
                                         int m,
                                         int n,
                                         double* A,
                                         int lda,
                                         double* Work,
                                         int* devIpiv,
                                         int* devInfo)
{
  rocsolverStatus_t success = rocsolver_dgetrf(handle, m, n, A, lda, devIpiv, devInfo);
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_getrf(rocsolverHandle_t handle,
                                         int m,
                                         int n,
                                         std::complex<float>* A,
                                         int lda,
                                         std::complex<float>* Work,
                                         int* devIpiv,
                                         int* devInfo)
{
  rocsolverStatus_t success =
      rocsolver_cgetrf(handle, m, n, reinterpret_cast<rocblas_float_complex*>(A), lda, devIpiv, devInfo);
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_getrf(rocsolverHandle_t handle,
                                         int m,
                                         int n,
                                         std::complex<double>* A,
                                         int lda,
                                         std::complex<double>* Work,
                                         int* devIpiv,
                                         int* devInfo)
{
  rocsolverStatus_t success =
      rocsolver_zgetrf(handle, m, n, reinterpret_cast<rocblas_double_complex*>(A), lda, devIpiv, devInfo);
  hipDeviceSynchronize();
  return success;
}


// getrs
// TODO: FDM rocsolver changes type of A to float *A, rather than const float *A
inline rocsolverStatus_t rocsolver_getrs(rocsolverHandle_t handle,
                                         rocblas_operation_ trans,
                                         int n,
                                         int nrhs,
                                         float* A,
                                         int lda,
                                         const int* devIpiv,
                                         float* B,
                                         int ldb,
                                         int* devInfo)
{
  rocsolverStatus_t success = rocsolver_sgetrs(handle, trans, n, nrhs, A, lda, devIpiv, B, ldb);
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_getrs(rocsolverHandle_t handle,
                                         rocblas_operation_ trans,
                                         int n,
                                         int nrhs,
                                         double* A,
                                         int lda,
                                         const int* devIpiv,
                                         double* B,
                                         int ldb,
                                         int* devInfo)
{
  rocsolverStatus_t success = rocsolver_dgetrs(handle, trans, n, nrhs, A, lda, devIpiv, B, ldb);
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_getrs(rocsolverHandle_t handle,
                                         rocblas_operation_ trans,
                                         int n,
                                         int nrhs,
                                         std::complex<float>* A,
                                         int lda,
                                         const int* devIpiv,
                                         std::complex<float>* B,
                                         int ldb,
                                         int* devInfo)
{
  rocsolverStatus_t success = rocsolver_cgetrs(handle, trans, n, nrhs, reinterpret_cast<rocblas_float_complex*>(A), lda,
                                               devIpiv, reinterpret_cast<rocblas_float_complex*>(B), ldb);
  // Info isn't returned from ?getrs.
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_getrs(rocsolverHandle_t handle,
                                         rocblas_operation_ trans,
                                         int n,
                                         int nrhs,
                                         std::complex<double>* A,
                                         int lda,
                                         const int* devIpiv,
                                         std::complex<double>* B,
                                         int ldb,
                                         int* devInfo)
{
  rocsolverStatus_t success = rocsolver_zgetrs(handle, trans, n, nrhs, reinterpret_cast<rocblas_double_complex*>(A),
                                               lda, devIpiv, reinterpret_cast<rocblas_double_complex*>(B), ldb);
  hipDeviceSynchronize();
  return success;
}

//geqrf_bufferSize
inline rocsolverStatus_t rocsolver_geqrf_bufferSize(rocsolverHandle_t handle,
                                                    int m,
                                                    int n,
                                                    float* A,
                                                    int lda,
                                                    int* Lwork)
{
  rocsolverStatus_t success;
  //throw std::runtime_error("Error: geqrf_bufferSize does not exist.");
  *Lwork  = 0;
  success = rocblas_status_success;
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_geqrf_bufferSize(rocsolverHandle_t handle,
                                                    int m,
                                                    int n,
                                                    double* A,
                                                    int lda,
                                                    int* Lwork)
{
  //rocsolverStatus_t success =
  //rocsolver_dgeqrf_bufferSize(handle,m,n,A,lda,Lwork);
  rocsolverStatus_t success;
  success = rocblas_status_success;
  *Lwork  = 0;
  //throw std::runtime_error("Error: geqrf_bufferSize does not exist.");
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_geqrf_bufferSize(rocsolverHandle_t handle,
                                                    int m,
                                                    int n,
                                                    std::complex<float>* A,
                                                    int lda,
                                                    int* Lwork)
{
  //rocsolverStatus_t success =
  //rocsolver_cgeqrf_bufferSize(handle,m,n,reinterpret_cast<rocblas_float_complex *>(A),lda,Lwork);
  rocsolverStatus_t success;
  success = rocblas_status_success;
  *Lwork  = 0;
  //throw std::runtime_error("Error: geqrf_bufferSize does not exist.");
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_geqrf_bufferSize(rocsolverHandle_t handle,
                                                    int m,
                                                    int n,
                                                    std::complex<double>* A,
                                                    int lda,
                                                    int* Lwork)
{
  //rocsolverStatus_t success =
  //rocsolver_zgeqrf_bufferSize(handle,m,n,reinterpret_cast<rocblas_double_complex *>(A),lda,Lwork);
  rocsolverStatus_t success;
  success = rocblas_status_success;
  *Lwork  = 0;
  //throw std::runtime_error("Error: geqrf_bufferSize does not exist.");
  hipDeviceSynchronize();
  return success;
}

//geqrf
inline rocsolverStatus_t rocsolver_geqrf(rocsolverHandle_t handle, int m, int n, float* A, int lda, float* ipiv)
{
  rocsolverStatus_t success = rocsolver_sgeqrf(handle, m, n, A, lda, ipiv);
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_geqrf(rocsolverHandle_t handle, int m, int n, double* A, int lda, double* ipiv)
{
  rocsolverStatus_t success = rocsolver_dgeqrf(handle, m, n, A, lda, ipiv);
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_geqrf(rocsolverHandle_t handle,
                                         int m,
                                         int n,
                                         std::complex<float>* A,
                                         int lda,
                                         std::complex<float>* ipiv)
{
  rocsolverStatus_t success = rocsolver_cgeqrf(handle, m, n, reinterpret_cast<rocblas_float_complex*>(A), lda,
                                               reinterpret_cast<rocblas_float_complex*>(ipiv));
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_geqrf(rocsolverHandle_t handle,
                                         int m,
                                         int n,
                                         std::complex<double>* A,
                                         int lda,
                                         std::complex<double>* ipiv)
{
  rocsolverStatus_t success = rocsolver_zgeqrf(handle, m, n, reinterpret_cast<rocblas_double_complex*>(A), lda,
                                               reinterpret_cast<rocblas_double_complex*>(ipiv));
  hipDeviceSynchronize();
  return success;
}


//gqr_bufferSize
// TODO: FDM these don't exist in rocsolver
inline rocsolverStatus_t rocsolver_gqr_bufferSize(rocsolverHandle_t handle,
                                                  int m,
                                                  int n,
                                                  int k,
                                                  const float* A,
                                                  int lda,
                                                  int* lwork)
{
  rocsolverStatus_t success;
  success = rocblas_status_success;
  *lwork  = m;
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_gqr_bufferSize(rocsolverHandle_t handle,
                                                  int m,
                                                  int n,
                                                  int k,
                                                  const double* A,
                                                  int lda,
                                                  int* lwork)
{
  //rocsolverStatus_t success =
  //rocsolver_dorgqr_bufferSize(handle,m,n,k,A,lda,A,lwork);
  rocsolverStatus_t success;
  *lwork  = m;
  success = rocblas_status_success;
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_gqr_bufferSize(rocsolverHandle_t handle,
                                                  int m,
                                                  int n,
                                                  int k,
                                                  const std::complex<float>* A,
                                                  int lda,
                                                  int* lwork)
{
  rocsolverStatus_t success;
  success = rocblas_status_success;
  *lwork  = m;
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_gqr_bufferSize(rocsolverHandle_t handle,
                                                  int m,
                                                  int n,
                                                  int k,
                                                  const std::complex<double>* A,
                                                  int lda,
                                                  int* lwork)
{
  rocsolverStatus_t success;
  success = rocblas_status_success;
  *lwork  = m;
  hipDeviceSynchronize();
  return success;
}

//gqr
inline rocsolverStatus_t rocsolver_gqr(rocsolverHandle_t handle, int m, int n, int k, float* A, int lda, float* ipiv)
{
  rocsolverStatus_t success = rocsolver_sorgqr(handle, m, n, k, A, lda, ipiv);
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_gqr(rocsolverHandle_t handle, int m, int n, int k, double* A, int lda, double* ipiv)
{
  rocsolverStatus_t success = rocsolver_dorgqr(handle, m, n, k, A, lda, ipiv);
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_gqr(rocsolverHandle_t handle,
                                       int m,
                                       int n,
                                       int k,
                                       std::complex<float>* A,
                                       int lda,
                                       std::complex<float>* ipiv)
{
  rocsolverStatus_t success = rocsolver_cungqr(handle, m, n, k, reinterpret_cast<rocblas_float_complex*>(A), lda,
                                               reinterpret_cast<rocblas_float_complex*>(ipiv));
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_gqr(rocsolverHandle_t handle,
                                       int m,
                                       int n,
                                       int k,
                                       std::complex<double>* A,
                                       int lda,
                                       std::complex<double>* ipiv)
{
  rocsolverStatus_t success = rocsolver_zungqr(handle, m, n, k, reinterpret_cast<rocblas_double_complex*>(A), lda,
                                               reinterpret_cast<rocblas_double_complex*>(ipiv));
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_gqr_strided(rocsolverHandle_t handle,
                                               const int m,
                                               const int n,
                                               const int k,
                                               std::complex<double>* A,
                                               const int lda,
                                               const int Astride,
                                               std::complex<double>* tau,
                                               const int tstride,
                                               std::complex<double>* work,
                                               int lwork,
                                               int* devInfo,
                                               int batchsize)
{
  // TODO FDM: double free with emplace_back(hipStream_t)?
  //using qmc_hip::afqmc_hip_streams;
  //auto t = hipStream_t{};
  //afqmc_hip_streams.push_back(hipStream_t{});
  //hipStreamCreate(&(afqmc_hip_streams.back()));
  //hipStreamDestroy(afqmc_hip_streams.back());
  //if(afqmc_hip_streams.size() < batchsize) {
  //int n0=afqmc_hip_streams.size();
  //for(int n=n0; n<batchsize; n++) {
  //auto t = hipStream_t{};
  //afqmc_hip_streams.emplace_back(hipStream_t{});
  //hipStreamCreate(&(afqmc_hip_streams.back()));
  //}
  //}

  //hipStream_t s0;
  //qmc_hip::hipsolver_check(rocsolver_get_stream(handle,&s0), "rocsolver_get_stream");

  rocsolverStatus_t success = rocblas_status_success;
  for (int i = 0; i < batchsize; i++)
  {
    //qmc_hip::hipsolver_check(rocsolver_set_stream(handle,afqmc_hip_streams[i]), "rocsolver_get_stream");
    success = rocsolver_zungqr(handle, m, n, k, reinterpret_cast<rocblas_double_complex*>(A) + i * Astride, lda,
                               reinterpret_cast<rocblas_double_complex*>(tau) + i * tstride);
    qmc_hip::hip_check(hipDeviceSynchronize(), "rocsolver_gqr_strided");
    qmc_hip::hipsolver_check(success, "rocsolver_gqr_strided");
  }
  qmc_hip::hip_check(hipDeviceSynchronize(), "rocsolver_gqr_strided");
  qmc_hip::hip_check(hipGetLastError(), "rocsolver_gqr_strided");
  //qmc_hip::hipsolver_check(rocsolver_set_stream(handle,s0), "rocsolver_setStream");

  return success;
}

// geqrf_batched
inline rocsolverStatus_t rocsolver_geqrf_batched(rocsolverHandle_t handle,
                                                 int m,
                                                 int n,
                                                 std::complex<double>** A,
                                                 int lda,
                                                 std::complex<double>* ipiv,
                                                 int strideP,
                                                 int batch_count)
{
  rocsolverStatus_t success =
      rocsolver_zgeqrf_batched(handle, m, n, reinterpret_cast<rocblas_double_complex* const*>(A), lda,
                               reinterpret_cast<rocblas_double_complex*>(ipiv), strideP, batch_count);
  hipDeviceSynchronize();
  return success;
}
inline rocsolverStatus_t rocsolver_geqrf_batched(rocsolverHandle_t handle,
                                                 int m,
                                                 int n,
                                                 std::complex<float>** A,
                                                 int lda,
                                                 std::complex<float>* ipiv,
                                                 int strideP,
                                                 int batch_count)
{
  rocsolverStatus_t success =
      rocsolver_cgeqrf_batched(handle, m, n, reinterpret_cast<rocblas_float_complex* const*>(A), lda,
                               reinterpret_cast<rocblas_float_complex*>(ipiv), strideP, batch_count);
  hipDeviceSynchronize();
  return success;
}
inline rocsolverStatus_t rocsolver_geqrf_batched(rocsolverHandle_t handle,
                                                 int m,
                                                 int n,
                                                 float** A,
                                                 int lda,
                                                 float* ipiv,
                                                 int strideP,
                                                 int batch_count)
{
  rocsolverStatus_t success = rocsolver_sgeqrf_batched(handle, m, n, A, lda, ipiv, strideP, batch_count);
  hipDeviceSynchronize();
  return success;
}
inline rocsolverStatus_t rocsolver_geqrf_batched(rocsolverHandle_t handle,
                                                 int m,
                                                 int n,
                                                 double** A,
                                                 int lda,
                                                 double* ipiv,
                                                 int strideP,
                                                 int batch_count)
{
  rocsolverStatus_t success = rocsolver_dgeqrf_batched(handle, m, n, A, lda, ipiv, strideP, batch_count);
  hipDeviceSynchronize();
  return success;
}

//gesvd_bufferSize
inline rocsolverStatus_t rocsolver_gesvd_bufferSize(rocsolverHandle_t handle, int m, int n, float* A, int* Lwork)
{
  //rocsolverStatus_t success =
  //rocsolver_sgesvd_bufferSize(handle,m,n,Lwork);
  rocsolverStatus_t success;
  success = rocblas_status_success;
  *Lwork  = 0;
  //throw std::runtime_error("Error: rocsolver_gesvd_bufferSize not implemented in rocsolver.");
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_gesvd_bufferSize(rocsolverHandle_t handle, int m, int n, double* A, int* Lwork)
{
  //rocsolverStatus_t success =
  //rocsolver_dgesvd_bufferSize(handle,m,n,Lwork);
  rocsolverStatus_t success;
  //throw std::runtime_error("Error: rocsolver_gesvd_bufferSize not implemented in rocsolver.");
  success = rocblas_status_success;
  *Lwork  = 0;
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_gesvd_bufferSize(rocsolverHandle_t handle,
                                                    int m,
                                                    int n,
                                                    std::complex<float>* A,
                                                    int* Lwork)
{
  //rocsolverStatus_t success =
  //rocsolver_cgesvd_bufferSize(handle,m,n,Lwork);
  rocsolverStatus_t success;
  //throw std::runtime_error("Error: rocsolver_gesvd_bufferSize not implemented in rocsolver.");
  success = rocblas_status_success;
  *Lwork  = 0;
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_gesvd_bufferSize(rocsolverHandle_t handle,
                                                    int m,
                                                    int n,
                                                    std::complex<double>* A,
                                                    int* Lwork)
{
  //rocsolverStatus_t success =
  //rocsolver_zgesvd_bufferSize(handle,m,n,Lwork);
  rocsolverStatus_t success;
  //throw std::runtime_error("Error: rocsolver_gesvd_bufferSize not implemented in rocsolver.");
  success = rocblas_status_success;
  *Lwork  = 0;
  hipDeviceSynchronize();
  return success;
}

//gesvd
inline rocsolverStatus_t rocsolver_gesvd(rocsolverHandle_t handle,
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
  //rocsolverStatus_t success =
  //rocsolver_sgesvd(handle,jobu,jobvt,m,n,A,lda,S,U,ldu,VT,ldvt,work,lwork,nullptr,devInfo);
  rocsolverStatus_t success;
  throw std::runtime_error("Error: rocsolver_gesvd not implemented in rocsolver.");
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_gesvd(rocsolverHandle_t handle,
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
  //rocsolverStatus_t success =
  //rocsolver_dgesvd(handle,jobu,jobvt,m,n,A,lda,S,U,ldu,VT,ldvt,work,lwork,nullptr,devInfo);
  rocsolverStatus_t success;
  throw std::runtime_error("Error: rocsolver_gesvd not implemented in rocsolver.");
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_gesvd(rocsolverHandle_t handle,
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
  //rocsolverStatus_t success =
  //rocsolver_cgesvd(handle,jobu,jobvt,m,n,
  //reinterpret_cast<rocblas_float_complex *>(A),lda,S,
  //reinterpret_cast<rocblas_float_complex *>(U),ldu,
  //reinterpret_cast<rocblas_float_complex *>(VT),ldvt,
  //reinterpret_cast<rocblas_float_complex *>(work),lwork,
  //nullptr,devInfo);
  rocsolverStatus_t success;
  throw std::runtime_error("Error: rocsolver_gesvd not implemented in rocsolver.");
  hipDeviceSynchronize();
  return success;
}

inline rocsolverStatus_t rocsolver_gesvd(rocsolverHandle_t handle,
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
  //rocsolverStatus_t success =
  //rocsolver_zgesvd(handle,jobu,jobvt,m,n,
  //reinterpret_cast<rocblas_double_complex *>(A),lda,S,
  //reinterpret_cast<rocblas_double_complex *>(U),ldu,
  //reinterpret_cast<rocblas_double_complex *>(VT),ldvt,
  //reinterpret_cast<rocblas_double_complex *>(work),lwork,
  //nullptr,devInfo);
  rocsolverStatus_t success;
  throw std::runtime_error("Error: rocsolver_gesvd not implemented in rocsolver.");
  hipDeviceSynchronize();
  return success;
}

} // namespace rocsolver

#endif
