//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Jakub Kurzak, jakurzak@amd.com, Advanced Micro Devices, Inc
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jakub Kurzak, jakurzak@amd.com, Advanced Micro Devices, Inc
//////////////////////////////////////////////////////////////////////////////////////


#ifndef HIPBLAS_HPP
#define HIPBLAS_HPP

#include <hipblas.h>
#include <hip/hip_complex.h>
#include "cuda2hip.h"

//------------------------------------------------------------------------------
hipblasStatus_t
hipblasCgemmBatched(hipblasHandle_t handle,
                    hipblasOperation_t transa,
                    hipblasOperation_t transb,
                    int m,
                    int n,
                    int k,
                    const hipComplex *alpha,
                    const hipComplex *const Aarray[],
                    int lda,
                    const hipComplex *const Barray[],
                    int ldb,
                    const hipComplex *beta,
                    hipComplex *const Carray[],
                    int ldc,
                    int batchCount);

hipblasStatus_t
hipblasZgemmBatched(hipblasHandle_t handle,
                    hipblasOperation_t transa,
                    hipblasOperation_t transb,
                    int m,
                    int n,
                    int k,
                    const hipDoubleComplex *alpha,
                    const hipDoubleComplex *const Aarray[],
                    int lda,
                    const hipDoubleComplex *const Barray[],
                    int ldb,
                    const hipDoubleComplex *beta,
                    hipDoubleComplex *const Carray[],
                    int ldc,
                    int batchCount);

//------------------------------------------------------------------------------
hipblasStatus_t
hipblasSgetrfBatched_(cublasHandle_t handle,
                      int n,
                      float *const A[],
                      int lda,
                      int *P,
                      int *info,
                      int batchSize);

hipblasStatus_t
hipblasDgetrfBatched_(cublasHandle_t handle,
                      int n,
                      double *const A[],
                      int lda,
                      int *P,
                      int *info,
                      int batchSize);

hipblasStatus_t
hipblasCgetrfBatched_(cublasHandle_t handle,
                      int n,
                      cuComplex *const A[],
                      int lda,
                      int *P,
                      int *info,
                      int batchSize);

hipblasStatus_t
hipblasZgetrfBatched_(cublasHandle_t handle,
                      int n,
                      cuDoubleComplex *const A[],
                      int lda,
                      int *P,
                      int *info,
                      int batchSize);

//------------------------------------------------------------------------------
hipblasStatus_t
hipblasSgetriBatched_(cublasHandle_t handle,
                      int n,
                      const float *const A[],
                      int lda,
                      const int *P,
                      float *const C[],
                      int ldc,
                      int *info,
                      int batchSize);

hipblasStatus_t
hipblasDgetriBatched_(cublasHandle_t handle,
                      int n,
                      const double *const A[],
                      int lda,
                      const int *P,
                      double *const C[],
                      int ldc,
                      int *info,
                      int batchSize);

hipblasStatus_t
hipblasCgetriBatched_(cublasHandle_t handle,
                      int n,
                      const cuComplex *const A[],
                      int lda,
                      const int *P,
                      cuComplex *const C[],
                      int ldc,
                      int *info,
                      int batchSize);

hipblasStatus_t
hipblasZgetriBatched_(cublasHandle_t handle,
                      int n,
                      const cuDoubleComplex *const A[],
                      int lda,
                      const int *P,
                      cuDoubleComplex *const C[],
                      int ldc,
                      int *info,
                      int batchSize);

#endif /* HIPBLAS_HPP */
