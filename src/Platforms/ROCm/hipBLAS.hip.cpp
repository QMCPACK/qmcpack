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


#include "hipBLAS.hpp"
#include <rocsolver/rocsolver.h>

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
                    int batchCount)
{
    return hipblasCgemmBatched(handle,
                               transa,
                               transb,
                               m,
                               n,
                               k,
                               (const hipblasComplex *)alpha,
                               (const hipblasComplex *const *)Aarray,
                               lda,
                               (const hipblasComplex *const *)Barray,
                               ldb,
                               (const hipblasComplex *)beta,
                               (hipblasComplex *const *)Carray,
                               ldc,
                               batchCount);
}

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
                    int batchCount)
{
    return hipblasZgemmBatched(handle,
                               transa,
                               transb,
                               m,
                               n,
                               k,
                               (const hipblasDoubleComplex *)alpha,
                               (const hipblasDoubleComplex *const *)Aarray,
                               lda,
                               (const hipblasDoubleComplex *const *)Barray,
                               ldb,
                               (const hipblasDoubleComplex *)beta,
                               (hipblasDoubleComplex *const *)Carray,
                               ldc,
                               batchCount);
}

//------------------------------------------------------------------------------
hipblasStatus_t
hipblasSgetrfBatched_(hipblasHandle_t handle,
                      int n,
                      float *const A[],
                      int lda,
                      int *P,
                      int *info,
                      int batchSize)
{
    if (P != nullptr) {
        return (hipblasStatus_t)rocsolver_sgetrf_batched(
            (rocblas_handle)handle,
            (const rocblas_int)n,
            (const rocblas_int)n,
            (float* const *)A,
            (const rocblas_int)lda,
            (rocblas_int*)P,
            (const rocblas_stride)n,
            (rocblas_int*)info,
            (const rocblas_int)batchSize);
    }
    else {
        return (hipblasStatus_t)rocsolver_sgetrf_npvt_batched(
            (rocblas_handle)handle,
            (const rocblas_int)n,
            (const rocblas_int)n,
            (float* const *)A,
            (const rocblas_int)lda,
            (rocblas_int*)info,
            (const rocblas_int)batchSize);
    }
}

hipblasStatus_t
hipblasDgetrfBatched_(hipblasHandle_t handle,
                      int n,
                      double *const A[],
                      int lda,
                      int *P,
                      int *info,
                      int batchSize)
{
    if (P != nullptr) {
        return (hipblasStatus_t)rocsolver_dgetrf_batched(
            (rocblas_handle)handle,
            (const rocblas_int)n,
            (const rocblas_int)n,
            (double* const *)A,
            (const rocblas_int)lda,
            (rocblas_int*)P,
            (const rocblas_stride)n,
            (rocblas_int*)info,
            (const rocblas_int)batchSize);
    }
    else {
        return (hipblasStatus_t)rocsolver_dgetrf_npvt_batched(
            (rocblas_handle)handle,
            (const rocblas_int)n,
            (const rocblas_int)n,
            (double* const *)A,
            (const rocblas_int)lda,
            (rocblas_int*)info,
            (const rocblas_int)batchSize);
    }
}

hipblasStatus_t
hipblasCgetrfBatched_(hipblasHandle_t handle,
                      int n,
                      hipComplex *const A[],
                      int lda,
                      int *P,
                      int *info,
                      int batchSize)
{
    if (P != nullptr) {
        return (hipblasStatus_t)rocsolver_cgetrf_batched(
            (rocblas_handle)handle,
            (const rocblas_int)n,
            (const rocblas_int)n,
            (rocblas_float_complex* const *)A,
            (const rocblas_int)lda,
            (rocblas_int*)P,
            (const rocblas_stride)n,
            (rocblas_int*)info,
            (const rocblas_int)batchSize);
    }
    else {
        return (hipblasStatus_t)rocsolver_cgetrf_npvt_batched(
            (rocblas_handle)handle,
            (const rocblas_int)n,
            (const rocblas_int)n,
            (rocblas_float_complex* const *)A,
            (const rocblas_int)lda,
            (rocblas_int*)info,
            (const rocblas_int)batchSize);
    }
}

hipblasStatus_t
hipblasZgetrfBatched_(hipblasHandle_t handle,
                      int n,
                      hipDoubleComplex *const A[],
                      int lda,
                      int *P,
                      int *info,
                      int batchSize)
{
    if (P != nullptr) {
        return (hipblasStatus_t)rocsolver_zgetrf_batched(
            (rocblas_handle)handle,
            (const rocblas_int)n,
            (const rocblas_int)n,
            (rocblas_double_complex* const *)A,
            (const rocblas_int)lda,
            (rocblas_int*)P,
            (const rocblas_stride)n,
            (rocblas_int*)info,
            (const rocblas_int)batchSize);
    }
    else {
        return (hipblasStatus_t)rocsolver_zgetrf_npvt_batched(
            (rocblas_handle)handle,
            (const rocblas_int)n,
            (const rocblas_int)n,
            (rocblas_double_complex* const *)A,
            (const rocblas_int)lda,
            (rocblas_int*)info,
            (const rocblas_int)batchSize);
    }
}

//------------------------------------------------------------------------------
__global__
static void identity(int n, float* const a_array[], int lda)
{
    float* const a = a_array[blockIdx.x];
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; i += blockDim.x)
            if (i+threadIdx.x < n)
                a[i+threadIdx.x + j*lda] =
                    (i+threadIdx.x == j) ? 1.0f : 0.0f;
}

__global__
static void identity(int n, double* const a_array[], int lda)
{
    double* const a = a_array[blockIdx.x];
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; i += blockDim.x)
            if (i+threadIdx.x < n)
                a[i+threadIdx.x + j*lda] =
                    (i+threadIdx.x == j) ? 1.0 : 0.0;
}

__global__
static void identity(int n, rocblas_float_complex* const a_array[], int lda)
{
    rocblas_float_complex one = {1.0f, 0.0f};
    rocblas_float_complex zero = {0.0f, 0.0f};
    rocblas_float_complex* const a = a_array[blockIdx.x];
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; i += blockDim.x)
            if (i+threadIdx.x < n)
                a[i+threadIdx.x + j*lda] =
                    (i+threadIdx.x == j) ? one : zero;
}

__global__
static void identity(int n, rocblas_double_complex* const a_array[], int lda)
{
    rocblas_double_complex one = {1.0, 0.0};
    rocblas_double_complex zero = {0.0, 0.0};
    rocblas_double_complex* const a = a_array[blockIdx.x];
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; i += blockDim.x)
            if (i+threadIdx.x < n)
                a[i+threadIdx.x + j*lda] =
                    (i+threadIdx.x == j) ? one : zero;
}

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
                      int batchSize)
{

    if (P != nullptr) {
        return hipblasSgetriBatched(handle,
                                    n,
                                    (float *const *)A,
                                    lda,
                                    (int *)P,
                                    (float *const *)C,
                                    ldc,
                                    info,
                                    batchSize);
    }
    else {
        identity<<<dim3(batchSize, 1, 1), dim3(64, 1, 1), 0, 0>>>(n, C, ldc);

        float one = 1.0f;
        rocblas_strsm_batched((rocblas_handle)handle,
                              rocblas_side_left,
                              rocblas_fill_lower,
                              rocblas_operation_none,
                              rocblas_diagonal_unit,
                              (rocblas_int)n,
                              (rocblas_int)n,
                              &one,
                              (const float* const *)A,
                              (rocblas_int)lda,
                              (float* const *)C,
                              (rocblas_int)ldc,
                              (rocblas_int)batchSize);

        rocblas_strsm_batched((rocblas_handle)handle,
                              rocblas_side_left,
                              rocblas_fill_upper,
                              rocblas_operation_none,
                              rocblas_diagonal_non_unit,
                              (rocblas_int)n,
                              (rocblas_int)n,
                              &one,
                              (const float* const *)A,
                              (rocblas_int)lda,
                              (float* const *)C,
                              (rocblas_int)ldc,
                              (rocblas_int)batchSize);

        return HIPBLAS_STATUS_SUCCESS;
    }
}

hipblasStatus_t
hipblasDgetriBatched_(hipblasHandle_t handle,
                      int n,
                      const double *const A[],
                      int lda,
                      const int *P,
                      double *const C[],
                      int ldc,
                      int *info,
                      int batchSize)
{
    if (P != nullptr) {
        return hipblasDgetriBatched(handle,
                                    n,
                                    (double *const *)A,
                                    lda,
                                    (int *)P,
                                    (double *const *)C,
                                    ldc,
                                    info,
                                    batchSize);
    }
    else {
        identity<<<dim3(batchSize, 1, 1), dim3(64, 1, 1), 0, 0>>>(n, C, ldc);

        double one = 1.0;
        rocblas_dtrsm_batched((rocblas_handle)handle,
                              rocblas_side_left,
                              rocblas_fill_lower,
                              rocblas_operation_none,
                              rocblas_diagonal_unit,
                              (rocblas_int)n,
                              (rocblas_int)n,
                              &one,
                              (const double* const *)A,
                              (rocblas_int)lda,
                              (double* const *)C,
                              (rocblas_int)ldc,
                              (rocblas_int)batchSize);

        rocblas_dtrsm_batched((rocblas_handle)handle,
                              rocblas_side_left,
                              rocblas_fill_upper,
                              rocblas_operation_none,
                              rocblas_diagonal_non_unit,
                              (rocblas_int)n,
                              (rocblas_int)n,
                              &one,
                              (const double* const *)A,
                              (rocblas_int)lda,
                              (double* const *)C,
                              (rocblas_int)ldc,
                              (rocblas_int)batchSize);

        return HIPBLAS_STATUS_SUCCESS;
    }
}

hipblasStatus_t
hipblasCgetriBatched_(hipblasHandle_t handle,
                      int n,
                      const hipComplex *const A[],
                      int lda,
                      const int *P,
                      hipComplex *const C[],
                      int ldc,
                      int *info,
                      int batchSize)
{
    if (P != nullptr) {
        return hipblasCgetriBatched(handle,
                                    n,
                                    (hipblasComplex *const *)A,
                                    lda,
                                    (int *)P,
                                    (hipblasComplex *const *)C,
                                    ldc,
                                    info,
                                    batchSize);
    }
    else {
        identity<<<dim3(batchSize, 1, 1), dim3(64, 1, 1), 0, 0>>>(
            n, (rocblas_float_complex *const *)C, ldc);

        rocblas_float_complex one = {1.0f, 0.0f};
        rocblas_ctrsm_batched((rocblas_handle)handle,
                              rocblas_side_left,
                              rocblas_fill_lower,
                              rocblas_operation_none,
                              rocblas_diagonal_unit,
                              (rocblas_int)n,
                              (rocblas_int)n,
                              &one,
                              (const rocblas_float_complex* const *)A,
                              (rocblas_int)lda,
                              (rocblas_float_complex* const *)C,
                              (rocblas_int)ldc,
                              (rocblas_int)batchSize);

        rocblas_ctrsm_batched((rocblas_handle)handle,
                              rocblas_side_left,
                              rocblas_fill_upper,
                              rocblas_operation_none,
                              rocblas_diagonal_non_unit,
                              (rocblas_int)n,
                              (rocblas_int)n,
                              &one,
                              (const rocblas_float_complex* const *)A,
                              (rocblas_int)lda,
                              (rocblas_float_complex* const *)C,
                              (rocblas_int)ldc,
                              (rocblas_int)batchSize);

        return HIPBLAS_STATUS_SUCCESS;
    }
}

hipblasStatus_t
hipblasZgetriBatched_(hipblasHandle_t handle,
                      int n,
                      const hipDoubleComplex *const A[],
                      int lda,
                      const int *P,
                      hipDoubleComplex *const C[],
                      int ldc,
                      int *info,
                      int batchSize)
{
    if (P != nullptr) {
        return hipblasZgetriBatched(handle,
                                    n,
                                    (hipblasDoubleComplex *const *)A,
                                    lda,
                                    (int *)P,
                                    (hipblasDoubleComplex *const *)C,
                                    ldc,
                                    info,
                                    batchSize);
    }
    else {
        identity<<<dim3(batchSize, 1, 1), dim3(64, 1, 1), 0, 0>>>(
            n, (rocblas_double_complex *const *)C, ldc);

        rocblas_double_complex one = {1.0, 0.0};
        rocblas_ztrsm_batched((rocblas_handle)handle,
                              rocblas_side_left,
                              rocblas_fill_lower,
                              rocblas_operation_none,
                              rocblas_diagonal_unit,
                              (rocblas_int)n,
                              (rocblas_int)n,
                              &one,
                              (const rocblas_double_complex* const *)A,
                              (rocblas_int)lda,
                              (rocblas_double_complex* const *)C,
                              (rocblas_int)ldc,
                              (rocblas_int)batchSize);

        rocblas_ztrsm_batched((rocblas_handle)handle,
                              rocblas_side_left,
                              rocblas_fill_upper,
                              rocblas_operation_none,
                              rocblas_diagonal_non_unit,
                              (rocblas_int)n,
                              (rocblas_int)n,
                              &one,
                              (const rocblas_double_complex* const *)A,
                              (rocblas_int)lda,
                              (rocblas_double_complex* const *)C,
                              (rocblas_int)ldc,
                              (rocblas_int)batchSize);

        return HIPBLAS_STATUS_SUCCESS;
    }
}
