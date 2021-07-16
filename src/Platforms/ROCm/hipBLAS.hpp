
#ifndef HIPBLAS_HPP
#define HIPBLAS_HPP

#include <hipblas.h>
#include <hip/hip_complex.h>

static inline hipblasStatus_t
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

static inline hipblasStatus_t
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

static inline hipblasStatus_t
hipblasCgetrfBatched(cublasHandle_t handle,
                     int n,
                     cuComplex *const A[],
                     int lda,
                     int *P,
                     int *info,
                     int batchSize)
{
    return hipblasCgetrfBatched(handle,
                                n,
                                (hipblasComplex *const *)A,
                                lda,
                                P,
                                info,
                                batchSize);
}

static inline hipblasStatus_t
hipblasZgetrfBatched(cublasHandle_t handle,
                     int n,
                     cuDoubleComplex *const A[],
                     int lda,
                     int *P,
                     int *info,
                     int batchSize)
{
    return hipblasZgetrfBatched(handle,
                                n,
                                (hipblasDoubleComplex *const *)A,
                                lda,
                                P,
                                info,
                                batchSize);
}

static inline hipblasStatus_t
hipblasCgetriBatched(cublasHandle_t handle,
                    int n,
                    const cuComplex *const A[],
                    int lda,
                    const int *P,
                    cuComplex *const C[],
                    int ldc,
                    int *info,
                    int batchSize)
{
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

static inline hipblasStatus_t
hipblasZgetriBatched(cublasHandle_t handle,
                    int n,
                    const cuDoubleComplex *const A[],
                    int lda,
                    const int *P,
                    cuDoubleComplex *const C[],
                    int ldc,
                    int *info,
                    int batchSize)
{
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

#endif /* HIPBLAS_HPP */
