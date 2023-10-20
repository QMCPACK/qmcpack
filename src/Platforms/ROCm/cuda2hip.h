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


#ifndef CUDA2HIP_H
#define CUDA2HIP_H

// cuBLAS to hipBLAS
#define CUBLAS_OP_N                     HIPBLAS_OP_N
#define CUBLAS_OP_T                     HIPBLAS_OP_T
#define CUBLAS_STATUS_ALLOC_FAILED      HIPBLAS_STATUS_ALLOC_FAILED
#define CUBLAS_STATUS_ARCH_MISMATCH     HIPBLAS_STATUS_ARCH_MISMATCH
#define CUBLAS_STATUS_EXECUTION_FAILED  HIPBLAS_STATUS_EXECUTION_FAILED
#define CUBLAS_STATUS_INTERNAL_ERROR    HIPBLAS_STATUS_INTERNAL_ERROR
#define CUBLAS_STATUS_INVALID_VALUE     HIPBLAS_STATUS_INVALID_VALUE
//#define CUBLAS_STATUS_LICENSE_ERROR     HIPBLAS_STATUS_LICENSE_ERROR
#define CUBLAS_STATUS_MAPPING_ERROR     HIPBLAS_STATUS_MAPPING_ERROR
#define CUBLAS_STATUS_NOT_INITIALIZED   HIPBLAS_STATUS_NOT_INITIALIZED
#define CUBLAS_STATUS_NOT_SUPPORTED     HIPBLAS_STATUS_NOT_SUPPORTED
#define CUBLAS_STATUS_SUCCESS           HIPBLAS_STATUS_SUCCESS

#define cublasComplex           hipblasComplex
#define cublasDoubleComplex     hipblasDoubleComplex
#define cublasHandle_t          hipblasHandle_t
#define cublasStatus_t          hipblasStatus_t
#define cublasCreate            hipblasCreate
#define cublasDestroy           hipblasDestroy
#define cublasSetStream         hipblasSetStream
#define cublasGetStream         hipblasGetStream
#define cublasOperation_t       hipblasOperation_t
#define cublasCgeam             hipblasCgeam
#define cublasCgemm             hipblasCgemm
#define cublasCgemmBatched      hipblasCgemmBatched
#define cublasCgetrfBatched     hipblasCgetrfBatched_
#define cublasCgetriBatched     hipblasCgetriBatched_
#define cublasDgeam             hipblasDgeam
#define cublasDgemm             hipblasDgemm
#define cublasDgemmBatched      hipblasDgemmBatched
#define cublasDgetrfBatched     hipblasDgetrfBatched_
#define cublasDgetriBatched     hipblasDgetriBatched_
#define cublasSgeam             hipblasSgeam
#define cublasSgemm             hipblasSgemm
#define cublasSgemmBatched      hipblasSgemmBatched
#define cublasSgetrfBatched     hipblasSgetrfBatched_
#define cublasSgetriBatched     hipblasSgetriBatched_
#define cublasZgeam             hipblasZgeam
#define cublasZgemm             hipblasZgemm
#define cublasZgemmBatched      hipblasZgemmBatched
#define cublasZgetrfBatched     hipblasZgetrfBatched_
#define cublasZgetriBatched     hipblasZgetriBatched_

// CUDA runtime to HIP runtime
#define cuComplex                       hipComplex
#define cuDoubleComplex                 hipDoubleComplex
#define cuCrealf                        hipCrealf
#define cuCimagf                        hipCimagf
#define cuCreal                         hipCreal
#define cuCimag                         hipCimag
#define make_cuComplex                  make_hipComplex
#define make_cuDoubleComplex            make_hipDoubleComplex
#define cudaAddressModeClamp            hipAddressModeClamp
#define cudaArray                       hipArray
#define cudaBindTextureToArray          hipBindTextureToArray
#define cudaChannelFormatDesc           hipChannelFormatDesc
#define cudaChannelFormatKindFloat      hipChannelFormatKindFloat
#define cudaCreateChannelDesc           hipCreateChannelDesc
#define cudaDeviceProp                  hipDeviceProp_t
#define cudaDeviceReset                 hipDeviceReset
#define cudaDeviceSynchronize           hipDeviceSynchronize
#define cudaError_t                     hipError_t
#define cudaEvent_t                     hipEvent_t
#define cudaEventCreate                 hipEventCreate
#define cudaEventCreateWithFlags        hipEventCreateWithFlags
#define cudaEventDestroy                hipEventDestroy
#define cudaEventDisableTiming          hipEventDisableTiming
#define cudaEventElapsedTime            hipEventElapsedTime
#define cudaEventRecord                 hipEventRecord
#define cudaEventSynchronize            hipEventSynchronize
#define cudaFilterModeLinear            hipFilterModeLinear
#define cudaFree                        hipFree
#define cudaFreeHost                    hipHostFree
#define cudaGetDevice                   hipGetDevice
#define cudaGetDeviceCount              hipGetDeviceCount
#define cudaGetDeviceProperties         hipGetDeviceProperties
#define cudaGetErrorName                hipGetErrorName
#define cudaGetErrorString              hipGetErrorString
#define cudaGetLastError                hipGetLastError
#define cudaPeekAtLastError             hipPeekAtLastError
#define cudaHostAlloc                   hipHostMalloc
#define cudaHostAllocMapped             hipHostMallocMapped
#define cudaPointerGetAttributes        hipPointerGetAttributes
#define cudaPointerAttributes           hipPointerAttribute_t
#define cudaMemoryTypeHost              hipMemoryTypeHost
#define cudaMemoryTypeDevice            hipMemoryTypeDevice
#define cudaIpcGetMemHandle             hipIpcGetMemHandle
#define cudaIpcMemHandle_t              hipIpcMemHandle_t
#define cudaIpcMemLazyEnablePeerAccess  hipIpcMemLazyEnablePeerAccess
#define cudaIpcOpenMemHandle            hipIpcOpenMemHandle
#define cudaMalloc                      hipMalloc
#define cudaMallocArray                 hipMallocArray
#define cudaMallocHost                  hipHostMalloc
#if defined(QMC_DISABLE_HIP_HOST_REGISTER)
#define cudaHostRegister(ptr, size, flags) hipSuccess
#define cudaHostUnregister(ptr) hipSuccess
#else
#define cudaHostRegister                hipHostRegister
#define cudaHostUnregister              hipHostUnregister
#endif
#define cudaHostRegisterDefault         hipHostRegisterDefault
#define cudaMallocManaged               hipMallocManaged
#define cudaMemAdvise                   hipMemAdvise
#define cudaMemAdviseSetAccessedBy      hipMemAdviseSetAccessedBy
#define cudaMemAdviseSetReadMostly      hipMemAdviseSetReadMostly
#define cudaMemAttachGlobal             hipMemAttachGlobal
#define cudaMemcpy                      hipMemcpy
#define cudaMemcpyAsync                 hipMemcpyAsync
#define cudaMemcpyDeviceToDevice        hipMemcpyDeviceToDevice
#define cudaMemcpyDeviceToHost          hipMemcpyDeviceToHost
#define cudaMemcpyHostToDevice          hipMemcpyHostToDevice
#define cudaMemcpyHostToHost            hipMemcpyHostToHost
#define cudaMemcpyToArrayAsync          hipMemcpyToArray
#define cudaMemcpyToSymbol              hipMemcpyToSymbol
#define cudaMemcpyToSymbolAsync         hipMemcpyToSymbolAsync
#define cudaMemset                      hipMemset
#define cudaMemGetInfo                  hipMemGetInfo
#define cudaMemPrefetchAsync            hipMemPrefetchAsync
#define cudaReadModeElementType         hipReadModeElementType
#define cudaSetDevice                   hipSetDevice
#define cudaStream_t                    hipStream_t
#define cudaStreamCreate                hipStreamCreate
#define cudaStreamDestroy               hipStreamDestroy
#define cudaStreamSynchronize           hipStreamSynchronize
#define cudaStreamWaitEvent             hipStreamWaitEvent
#define cudaSuccess                     hipSuccess

#define cudaDeviceSetLimit(limit, value) ;

#endif /* CUDA2HIP_H */
