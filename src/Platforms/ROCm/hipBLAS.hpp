//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
// Copyright(C) 2021 Advanced Micro Devices, Inc. All rights reserved.
//
// File developed by: Jakub Kurzak, jakurzak@amd.com, Advanced Micro Devices, Inc.
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jakub Kurzak, jakurzak@amd.com, Advanced Micro Devices, Inc.
//////////////////////////////////////////////////////////////////////////////////////


#ifndef HIPBLAS_HPP
#define HIPBLAS_HPP

#include <hipblas/hipblas.h>
#include <hip/hip_complex.h>

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
hipblasSgetrfBatched_(hipblasHandle_t handle,
                      int n,
                      float *const A[],
                      int lda,
                      int *P,
                      int *info,
                      int batchSize);

hipblasStatus_t
hipblasDgetrfBatched_(hipblasHandle_t handle,
                      int n,
                      double *const A[],
                      int lda,
                      int *P,
                      int *info,
                      int batchSize);

hipblasStatus_t
hipblasCgetrfBatched_(hipblasHandle_t handle,
                      int n,
                      hipComplex *const A[],
                      int lda,
                      int *P,
                      int *info,
                      int batchSize);

hipblasStatus_t
hipblasZgetrfBatched_(hipblasHandle_t handle,
                      int n,
                      hipDoubleComplex *const A[],
                      int lda,
                      int *P,
                      int *info,
                      int batchSize);

//------------------------------------------------------------------------------
hipblasStatus_t
hipblasSgetriBatched_(hipblasHandle_t handle,
                      int n,
                      const float *const A[],
                      int lda,
                      const int *P,
                      float *const C[],
                      int ldc,
                      int *info,
                      int batchSize);

hipblasStatus_t
hipblasDgetriBatched_(hipblasHandle_t handle,
                      int n,
                      const double *const A[],
                      int lda,
                      const int *P,
                      double *const C[],
                      int ldc,
                      int *info,
                      int batchSize);

hipblasStatus_t
hipblasCgetriBatched_(hipblasHandle_t handle,
                      int n,
                      const hipComplex *const A[],
                      int lda,
                      const int *P,
                      hipComplex *const C[],
                      int ldc,
                      int *info,
                      int batchSize);

hipblasStatus_t
hipblasZgetriBatched_(hipblasHandle_t handle,
                      int n,
                      const hipDoubleComplex *const A[],
                      int lda,
                      const int *P,
                      hipDoubleComplex *const C[],
                      int ldc,
                      int *info,
                      int batchSize);

#endif /* HIPBLAS_HPP */
