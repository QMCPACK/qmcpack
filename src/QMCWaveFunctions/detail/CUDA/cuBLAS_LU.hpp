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
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <cuComplex.h>

namespace qmcplusplus
{
namespace cuBLAS_LU
{
void computeInverseAndDetLog_batched(cublasHandle_t& h_cublas,
                                     cudaStream_t& hstream,
                                     const int n,
                                     const int lda,
                                     double** Ms,
                                     double** Cs,
                                     double* LU_diags,
                                     int* pivots,
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

void computeGetri_batched(cublasHandle_t& h_cublas,
                          const int n,
                          const int lda,
                          double* Ms[],
                          double* Cs[],
                          int* pivots,
                          int* infos,
                          const int batch_size);

} // namespace cuBLAS_LU
} // namespace qmcplusplus
#endif
