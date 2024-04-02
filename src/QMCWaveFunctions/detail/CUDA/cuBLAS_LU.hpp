//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_CUBLAS_LU_HPP
#define QMCPLUSPLUS_CUBLAS_LU_HPP

#include <complex>
#include <type_traits>
#include "config.h"
#include <CUDA/CUDAruntime.hpp>
#ifndef QMC_CUDA2HIP
#include <cublas_v2.h>
#include <cuComplex.h>
#else
#include <hipblas/hipblas.h>
#include <hip/hip_complex.h>
#include <ROCm/hipBLAS.hpp>
#endif

/** \file
 *  At the qmcplusplus cuBLAS_LU level all *, **, *[] are assumed to be to device
 *  addresses.
 */
namespace qmcplusplus
{
namespace cuBLAS_LU
{
/** Takes PsiM in column major layout and uses LU factorization to compute the log determinant and invPsiM.
 *  This is the call the QMCPACK should use.
 *
 *  \param[inout] Ms -       device pointers to pointers to Ms on input and to LU matrices on output
 *  \param[out]   Cs -     device pointers to memory space same size as M which over written with invM
 *  \param[in]    pivots -   pointer to n * nw ints allocated in device memory for pivots array.
 *  \param[in]    host_infos - pointer to nw ints allocated in pinned host memory for factorization infos
 *  \param[in]    infos -    pointer to nw ints allocated in device memory factorization infos
 *  \param[out]   log_dets - pointer device memory for nw log determinant values to be returned will be zeroed. 
 *  \param[in]    batch_size - if this changes over run a huge performance hit will be taken as memory allocation syncs device.
 *
 *  The host infos is an exception to this that may be changed in the future. The logic for this should probably be in
 *  the next class up. This would obviously split the computeInverseAndDetLog_batched call.
 */
template<typename T>
void computeInverseAndDetLog_batched(cublasHandle_t& h_cublas,
                                     cudaStream_t& hstream,
                                     const int n,
                                     const int lda,
                                     T* Ms[],
                                     T* Cs[],
                                     T* LU_diags,
                                     int* pivots,
                                     int* host_infos,
                                     int* infos,
                                     std::complex<double>* log_dets,
                                     const int batch_size);

template<typename T>
void computeGetrf_batched(cublasHandle_t& h_cublas,
                          cudaStream_t& hstream,
                          const int n,
                          const int lda,
                          T* Ms[],
                          int* pivots,
                          int* host_infos,
                          int* infos,
                          const int batch_size);

template<typename T>
void computeLogDet_batched(cudaStream_t& hstream,
                           const int n,
                           const int lda,
                           T** Ms,
                           const int* pivots,
                           std::complex<double>* logdets,
                           const int batch_size);

template<typename T>
void computeGetri_batched(cublasHandle_t& h_cublas,
                          cudaStream_t& hstream,
                          const int n,
                          const int lda,
                          T* Ms[],
                          T* Cs[],
                          int* pivots,
                          int* host_infos,
                          int* infos,
                          const int batch_size);

extern template void computeInverseAndDetLog_batched<double>(cublasHandle_t& h_cublas,
                                                             cudaStream_t& hstream,
                                                             const int n,
                                                             const int lda,
                                                             double* Ms[],
                                                             double* Cs[],
                                                             double* LU_diags,
                                                             int* pivots,
                                                             int* host_infos,
                                                             int* infos,
                                                             std::complex<double>* log_dets,
                                                             const int batch_size);

extern template void computeInverseAndDetLog_batched<std::complex<double>>(cublasHandle_t& h_cublas,
                                                                           cudaStream_t& hstream,
                                                                           const int n,
                                                                           const int lda,
                                                                           std::complex<double>* Ms[],
                                                                           std::complex<double>* Cs[],
                                                                           std::complex<double>* LU_diags,
                                                                           int* pivots,
                                                                           int* host_infos,
                                                                           int* infos,
                                                                           std::complex<double>* log_dets,
                                                                           const int batch_size);

} // namespace cuBLAS_LU
} // namespace qmcplusplus
#endif
