
#ifndef CUDA2HIP_H
#define CUDA2HIP_H

#define CUBLAS_OP_N                     HIPBLAS_OP_N
#define CUBLAS_STATUS_ALLOC_FAILED      HIPBLAS_STATUS_ALLOC_FAILED
#define CUBLAS_STATUS_ARCH_MISMATCH     HIPBLAS_STATUS_ARCH_MISMATCH
#define CUBLAS_STATUS_EXECUTION_FAILED  HIPBLAS_STATUS_EXECUTION_FAILED
#define CUBLAS_STATUS_INTERNAL_ERROR    HIPBLAS_STATUS_INTERNAL_ERROR
#define CUBLAS_STATUS_INVALID_VALUE     HIPBLAS_STATUS_INVALID_VALUE
#define CUBLAS_STATUS_LICENSE_ERROR     HIPBLAS_STATUS_LICENSE_ERROR
#define CUBLAS_STATUS_MAPPING_ERROR     HIPBLAS_STATUS_MAPPING_ERROR
#define CUBLAS_STATUS_NOT_INITIALIZED   HIPBLAS_STATUS_NOT_INITIALIZED
#define CUBLAS_STATUS_NOT_SUPPORTED     HIPBLAS_STATUS_NOT_SUPPORTED
#define CUBLAS_STATUS_SUCCESS           HIPBLAS_STATUS_SUCCESS

#define cublasCgemmBatched      hipblasCgemmBatched
#define cublasCgetrfBatched     hipblasCgetrfBatched
#define cublasCgetrfBatched     hipblasCgetrfBatched
#define cublasCgetriBatched     hipblasCgetriBatched
#define cublasCgetriBatched     hipblasCgetriBatched
#define cublasComplex           hipblasComplex
#define cublasCreate            hipblasCreate
#define cublasDestroy           hipblasDestroy
#define cublasDgemmBatched      hipblasDgemmBatched
#define cublasDgetrfBatched     hipblasDgetrfBatched
#define cublasDgetriBatched     hipblasDgetriBatched
#define cublasDgetriBatched     hipblasDgetriBatched
#define cublasDoubleComplex     hipblasDoubleComplex
#define cublasHandle_t          hipblasHandle_t
#define cublasSgemmBatched      hipblasSgemmBatched
#define cublasSgetrfBatched     hipblasSgetrfBatched
#define cublasSgetriBatched     hipblasSgetriBatched
#define cublasSgetriBatched     hipblasSgetriBatched
#define cublasStatus_t          hipblasStatus_t
#define cublasZgemmBatched      hipblasZgemmBatched
#define cublasZgetrfBatched     hipblasZgetrfBatched
#define cublasZgetrfBatched     hipblasZgetrfBatched
#define cublasZgetriBatched     hipblasZgetriBatched
#define cublasZgetriBatched     hipblasZgetriBatched

#define cuComplex                   hipComplex
#define cudaAddressModeClamp        hipAddressModeClamp
#define cudaArray                   hipArray
#define cudaBindTextureToArray      hipBindTextureToArray
#define cudaChannelFormatDesc       hipChannelFormatDesc
#define cudaChannelFormatKindFloat  hipChannelFormatKindFloat
#define cudaCreateChannelDesc       hipCreateChannelDesc
#define cudaDeviceProp              hipDeviceProp_t
#define cudaDeviceReset             hipDeviceReset
#define cudaDeviceSynchronize       hipDeviceSynchronize
#define cudaError_t                 hipError_t
#define cudaEvent_t                 hipEvent_t
#define cudaEventCreate             hipEventCreate
#define cudaEventCreateWithFlags    hipEventCreateWithFlags
#define cudaEventDestroy            hipEventDestroy
#define cudaEventDisableTiming      hipEventDisableTiming
#define cudaEventElapsedTime        hipEventElapsedTime
#define cudaEventRecord             hipEventRecord
#define cudaEventSynchronize        hipEventSynchronize
#define cudaFilterModeLinear        hipFilterModeLinear
#define cudaFree                    hipFree
#define cudaFreeHost                hipHostFree
#define cudaGetDevice               hipGetDevice
#define cudaGetDeviceCount          hipGetDeviceCount
#define cudaGetDeviceProperties     hipGetDeviceProperties
#define cudaGetErrorString          hipGetErrorString
#define cudaGetLastError            hipGetLastError
#define cudaHostAlloc               hipHostMalloc
#define cudaHostAllocMapped         hipHostMallocMapped
#define cudaIpcMemHandle_t          hipIpcMemHandle_t
#define cudaMalloc                  hipMalloc
#define cudaMallocArray             hipMallocArray
#define cudaMallocManaged           hipMallocManaged
#define cudaMemAdvise               hipMemAdvise
#define cudaMemAdviseSetAccessedBy  hipMemAdviseSetAccessedBy
#define cudaMemAdviseSetReadMostly  hipMemAdviseSetReadMostly
#define cudaMemAttachGlobal         hipMemAttachGlobal
#define cudaMemcpy                  hipMemcpy
#define cudaMemcpyAsync             hipMemcpyAsync
#define cudaMemcpyDeviceToDevice    hipMemcpyDeviceToDevice
#define cudaMemcpyDeviceToHost      hipMemcpyDeviceToHost
#define cudaMemcpyHostToDevice      hipMemcpyHostToDevice
#define cudaMemcpyHostToHost        hipMemcpyHostToHost
#define cudaMemcpyToArrayAsync      hipMemcpyToArray
#define cudaMemcpyToSymbol          hipMemcpyToSymbol
#define cudaMemcpyToSymbolAsync     hipMemcpyToSymbolAsync
#define cudaMemPrefetchAsync        hipMemPrefetchAsync
#define cudaReadModeElementType     hipReadModeElementType
#define cudaSetDevice               hipSetDevice
#define cudaStream_t                hipStream_t
#define cudaStreamCreate            hipStreamCreate
#define cudaStreamDestroy           hipStreamDestroy
#define cudaStreamSynchronize       hipStreamSynchronize
#define cudaStreamWaitEvent         hipStreamWaitEvent
#define cudaSuccess                 hipSuccess
#define cuDoubleComplex             hipDoubleComplex
#define make_cuComplex              make_hipComplex
#define make_cuDoubleComplex        make_hipDoubleComplex

#define cudaDeviceSetLimit(limit, falue) ;

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

#endif /* CUDA2HIP_H */
