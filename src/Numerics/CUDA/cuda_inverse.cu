//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//		      Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign 
//    		      Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ying Wai Li, yingwaili@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by:  Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////




// ============ Matrix inversion using cuBLAS library ============ //
//
// To compile as a standalone test program:
//
// 1. Make sure libcublas.so is in the search path
// 2. cd to build/ directory
// 3. For real numbers, compile with
//    nvcc -o cuda_inverse -arch=sm_35 -lcublas -DCUDA_TEST_MAIN
//         ../src/Numerics/CUDA/cuda_inverse.cu
// 
//    For complex numbers, compile with
//    nvcc -o cuda_inverse -arch=sm_35 -lcublas -DCUDA_TEST_MAIN
//         -DQMC_COMPLEX=1 ../src/Numerics/CUDA/cuda_inverse.cu
//
// =============================================================== //

#include <cstdio>
#include <unistd.h>
#include <sstream>
#include <vector>
#include <iostream>
#include <complex>
#include <cuda.h>
#include <cublas_v2.h>
#include <cuComplex.h>

#define CONVERT_BS 256
#define INVERSE_BS 16


void
callAndCheckError (cudaError_t cudaFunc, const int line)
{
  if (cudaFunc != cudaSuccess)
  {
    fprintf(stderr, "CUDA error in %s, line %d \n", __FILE__, line); 
    fprintf(stderr, "CUDA error message : %s \n", cudaGetErrorString(cudaFunc));
    fflush(stderr);
    abort();
  }
}

void
callAndCheckError(cublasStatus_t cublasFunc, const int line)
{
  if (cublasFunc != CUBLAS_STATUS_SUCCESS)
  {
    fprintf(stderr, "CUBLAS error in %s, line %d \n", __FILE__, line);
    fprintf(stderr, "CUBLAS error message: ");
    switch (cublasFunc)
    {
       case CUBLAS_STATUS_NOT_INITIALIZED:
         fprintf(stderr, "CUBLAS_STATUS_NOT_INITIALIZED\n");
         break;
       case CUBLAS_STATUS_ALLOC_FAILED:
         fprintf(stderr, "CUBLAS_STATUS_ALLOC_FAILED\n");
         break;
       case CUBLAS_STATUS_INVALID_VALUE:
         fprintf(stderr, "CUBLAS_STATUS_INVALID_VALUE\n");
         break;
       case CUBLAS_STATUS_ARCH_MISMATCH:
         fprintf(stderr, "CUBLAS_STATUS_ARCH_MISMATCH\n");
         break;
       case CUBLAS_STATUS_MAPPING_ERROR:
         fprintf(stderr, "CUBLAS_STATUS_MAPPING_ERROR\n");
         break;
       case CUBLAS_STATUS_EXECUTION_FAILED:
         fprintf(stderr, "CUBLAS_STATUS_EXECUTION_FAILED\n");
         break;
       case CUBLAS_STATUS_INTERNAL_ERROR:
         fprintf(stderr, "CUBLAS_STATUS_INTERNAL_ERROR\n");
         break;
#if (CUDA_VERSION >= 6050)
       case CUBLAS_STATUS_NOT_SUPPORTED:
         fprintf(stderr, "CUBLAS_STATUS_NOT_SUPPORTED\n");
         break;
       case CUBLAS_STATUS_LICENSE_ERROR:
         fprintf(stderr, "CUBLAS_STATUS_LICENSE_ERROR\n");
         break;
#endif
       default:
         fprintf(stderr, "unknown\n");
    }
    fflush(stderr);
    abort();
  }
}

// Convert matrix elements from one type (Tsrc) in the source matrix to
// another type (Tdest) and put them in the destination matrix
// (assumed src and dest have the same dimensions)
template <typename Tdest, typename Tsrc>
__global__ void
convert (Tdest **dest_list, Tsrc **src_list, int len)
{
  __shared__ Tsrc *mysrc;
  __shared__ Tdest *mydest;
  if (threadIdx.x == 0)
  {
    mysrc = src_list[blockIdx.y];
    mydest = dest_list[blockIdx.y];
  }
  __syncthreads();
  int i = blockIdx.x * CONVERT_BS + threadIdx.x;
  if (i < len)
    mydest[i] = (Tdest) mysrc[i];
}

// Convert for complex numbers
template <typename Tdest, typename Tdest2, typename Tsrc>
__global__ void
convert_complex (Tdest **dest_list, Tsrc **src_list, int len)
{
  __shared__ Tsrc *mysrc;
  __shared__ Tdest *mydest;
  if (threadIdx.x == 0)
  {
    mysrc = src_list[blockIdx.y];
    mydest = dest_list[blockIdx.y];
  }
  __syncthreads();
  int i = blockIdx.x * CONVERT_BS + threadIdx.x;
  if (i < len) {
    mydest[i].x = (Tdest2) mysrc[i].x;
    mydest[i].y = (Tdest2) mysrc[i].y;
  }
}

// C = A - B
template <typename T>
__global__ void
subtract (T **C, T **A, T **B, int len)
{
  __shared__ T *myA, *myB, *myC;
  if (threadIdx.x == 0)
  {
    myA = A[blockIdx.y];
    myB = B[blockIdx.y];
    myC = C[blockIdx.y];
  }
  __syncthreads();
  int i = blockIdx.x * CONVERT_BS + threadIdx.x;
  if (i < len)
    myC[i] = myA[i] - myB[i];
}

/** Calculate Lemma Matrix: I_k + V' * ( A^(-1) * U )
  * for each walker
  * -> returns L-U decomposed lemma matrix for easy determinant calculations and inverse calculation later
  */
void
cublas_lemma_mats (cublasHandle_t handle,
                   float *AList_d[], float *AWorkList_d[],
                   float *AinvList_d[], float *AinvkList_d[], float *U_d[],
                   float *lemma_d[], float *AinvUList_d[],
                   int k, int N, int nw, int RowStride)
{
  float one=1.0;
  float zero=0.0;
  // Calculate Lemma Matrix
  // V^-1 * A^(-1) * U
  // per walker: [k x N] * [N x k] = [k x k]
  callAndCheckError( cublasSgemmBatched( handle, CUBLAS_OP_N, CUBLAS_OP_N, k, k, N,
                                         &one,
                                         (const float**)AinvkList_d, RowStride,
                                         (const float**)U_d, RowStride, &zero,
                                         lemma_d, k,
                                         nw), __LINE__ );
  // Calculate - A^-1*dU
  dim3 dimBlockConvert (CONVERT_BS);
  dim3 dimGridConvert ((k*RowStride + (CONVERT_BS-1)) / CONVERT_BS, nw);
  // Calculate -dU=U(old)-U(new)
  subtract <<< dimGridConvert, dimBlockConvert >>> (AWorkList_d, AList_d, U_d, k*RowStride);
  // -A^(-1) * dU
  // per walker: [N x N] * [N x k] = [N x k]
#ifndef AINVU_TRANSPOSE
  callAndCheckError( cublasSgemmBatched( handle, CUBLAS_OP_N, CUBLAS_OP_N, N, k, N,
                                         &one,
                                         (const float**)AinvList_d, RowStride,
                                         (const float**)AWorkList_d, RowStride, &zero,
                                         AinvUList_d, RowStride,
                                         nw), __LINE__ );
#else
  // calculate AinvU as row major
  // per walker: [N x k]^T * [N x N]^T = [k x N] * [N x N] = [k x N]
  callAndCheckError( cublasSgemmBatched( handle, CUBLAS_OP_T, CUBLAS_OP_T, k, N, N,
                                         &one,
                                         (const float**)AWorkList_d, RowStride,
                                         (const float**)AinvList_d, RowStride, &zero,
                                         AinvUList_d, k,
                                         nw), __LINE__ );
#endif
//  cudaDeviceSynchronize();
}

void
cublas_ainv_row (cublasHandle_t handle,
                 float *AinvkList_d[], float *AWorkList_d[], float *AinvList_d[],
                 int k, int N, int nw, int RowStride)
{
  float one=1.0;
  // A^-1 - { A^-1 * dU  * Lemma^-1 } * { V' * A^(-1) }
  // per walker: [1 x N] - [1 x k] * [k x N] = [1 x N]
  callAndCheckError( cublasSgemmBatched( handle, CUBLAS_OP_N, CUBLAS_OP_N, 1, N, k,
                                         &one,
                                         (const float**)AWorkList_d, 1,
                                         (const float**)AinvkList_d, RowStride, &one,
                                         AinvList_d, 1,
                                         nw), __LINE__ );
}

void
cublas_ainv_row (cublasHandle_t handle,
                 double *AinvkList_d[], double *AWorkList_d[], double *AinvList_d[],
                 int k, int N, int nw, int RowStride)
{
  double one=1.0;
  // A^-1 - { A^-1 * dU  * Lemma^-1 } * { V' * A^(-1) }
  // per walker: [1 x N] - [1 x k] * [k x N] = [1 x N]
  callAndCheckError( cublasDgemmBatched( handle, CUBLAS_OP_N, CUBLAS_OP_N, 1, N, k,
                                         &one,
                                         (const double**)AWorkList_d, 1,
                                         (const double**)AinvkList_d, RowStride, &one,
                                         AinvList_d, 1,
                                         nw), __LINE__ );
}


// #define DEBUG_DELAYED
// #define USE_TRSM

void
cublas_smw_update (cublasHandle_t handle,
                   float *AinvkList_d[], float *AinvList_d[],
                   float *AinvUList_d[], float *AWorkList_d[],
                   float *lemma_inv[], float *lemma_lu[],
                   int *infoArray,
                   int k, int kd, int M, int N, int nw, int RowStride)
{
#ifdef DEBUG_DELAYED
  fprintf(stderr,"*** Sherman-Morrison-Woodbury Update (k = %i, %i walkers) ***\n",k,nw);
#endif
  int pitch=RowStride;
  if(M==1) pitch=1;
  float one=1.0;

  // LU decomposition needs to be updated
  callAndCheckError( cublasSgetrfBatched( handle, k, lemma_lu, kd, NULL,
                                          infoArray, nw), __LINE__ );

#ifdef USE_TRSM
  if(M==1)
  {
    // {-A^-1 * dU } * Lemma^(-1) => solve for y: Lemma y * (L * U) = (y * L) * U = -A^-1 * dU
    // z * U = -A^-1 *dU
    callAndCheckError( cublasStrsmBatched( handle, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT,
                                           M, k, &one, (const float**) lemma_lu, kd,
                                           AWorkList_d, pitch, nw), __LINE__ );
    // y * L = z => y = {-A^-1 * dU } * Lemma^(-1)
    callAndCheckError( cublasStrsmBatched( handle, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_UNIT,
                                           M, k, &one, (const float**) lemma_lu, kd,
                                           AWorkList_d, pitch, nw), __LINE__ );
    // A^-1 + { -A^-1 * dU *  Lemma^-1 } * { V' * A^(-1) }
    // per walker: [1 x N] - [1 x k] * [k x N] = [1 x N]
    callAndCheckError( cublasSgemmBatched( handle, CUBLAS_OP_N, CUBLAS_OP_N, M, N, k,
                                           &one,
                                           (const float**)AWorkList_d, M,
                                           (const float**)AinvkList_d, RowStride, &one,
                                           AinvList_d, pitch,
                                           nw), __LINE__ );
  } else
  {
    // Lemma^(-1) * V' * A^(-1) => solve for y: Lemma (L * U) * y = L * (U * y) = V' * A^(-1)
    // L * z = V' * A^(-1)
    callAndCheckError( cublasStrsmBatched( handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_UNIT,
                                           k, N, &one, (const float**) lemma_lu, kd,
                                           AWorkList_d, k, nw), __LINE__ );
    // U * y = z => y = Lemma^(-1) * V' * A^(-1)
    callAndCheckError( cublasStrsmBatched( handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT,
                                           k, N, &one, (const float**) lemma_lu, kd,
                                           AWorkList_d, k, nw), __LINE__ );
    // A^-1 + { -A^-1 * dU } * { Lemma^-1 * V' * A^(-1) }
    // per walker: [M x N] - [M x k] * [k x N] = [M x N]
    callAndCheckError( cublasSgemmBatched( handle, CUBLAS_OP_N, CUBLAS_OP_N, M, N, k,
                                           &one,
                                           (const float**)AinvUList_d, RowStride,
                                           (const float**)AWorkList_d, k, &one,
                                           AinvList_d, pitch,
                                           nw), __LINE__ );
  }
#else
  float zero=0.0;
  // Calculate Lemma Inverse and store it in lemma_d
  // per walker: [k x k]^-1 = [k x k]
  callAndCheckError( cublasSgetriBatched( handle, k, (const float**) lemma_lu, kd, NULL,
                                          lemma_inv, k, infoArray, nw), __LINE__ );
  // Calculate new A inverse using Sherman-Morrison-Woodbury formula
  if(M==1) // row update can use different order to save flops
  {
    // { -A^-1 * dU } * Lemma^-1
    // per walker: [M x k] * [k x k] = [M x k]
    callAndCheckError( cublasSgemmBatched( handle, CUBLAS_OP_N, CUBLAS_OP_N, M, k, k,
                                           &one,
                                           (const float**)AinvUList_d, RowStride,
                                           (const float**)lemma_inv, k, &zero,
                                           AWorkList_d, M,
                                           nw), __LINE__ );
    // A^-1 - { A^-1 * dU  * Lemma^-1 } * { V' * A^(-1) }
    // per walker: [M x N] - [M x k] * [k x N] = [M x N]
    callAndCheckError( cublasSgemmBatched( handle, CUBLAS_OP_N, CUBLAS_OP_N, M, N, k,
                                           &one,
                                           (const float**)AWorkList_d, M,
                                           (const float**)AinvkList_d, RowStride, &one,
                                           AinvList_d, pitch,
                                           nw), __LINE__ );
  } else
  {
    // Need to use this matrix order for the overall update as AinvList and AinvkList have overlapping memory
    // Lemma^-1 * V' * A^(-1)
    // per walker: [k x k] * [k x N] = [k x N]
    callAndCheckError( cublasSgemmBatched( handle, CUBLAS_OP_N, CUBLAS_OP_N, k, N, k,
                                           &one,
                                           (const float**)lemma_inv, k,
                                           (const float**)AinvkList_d, RowStride, &zero,
                                           AWorkList_d, k,
                                           nw), __LINE__ );
    // A^-1 + { -A^-1 * dU } * { Lemma^-1 * V' * A^(-1) }
    // per walker: [M x N] - [M x k] * [k x N] = [M x N]
    callAndCheckError( cublasSgemmBatched( handle, CUBLAS_OP_N, CUBLAS_OP_N, M, N, k,
                                           &one,
                                           (const float**)AinvUList_d, RowStride,
                                           (const float**)AWorkList_d, k, &one,
                                           AinvList_d, pitch,
                                           nw), __LINE__ );
  }
#endif
}

/** Calculate Lemma Matrix: I_k + V' * ( A^(-1) * U )
  * for each walker
  * -> returns L-U decomposed lemma matrix for easy determinant calculations and inverse calculation later
  */
void
cublas_lemma_mats (cublasHandle_t handle,
                   double *AList_d[], double *AWorkList_d[],
                   double *AinvList_d[], double *AinvkList_d[], double *U_d[],
                   double *lemma_d[], double *AinvUList_d[],
                   int k, int N, int nw, int RowStride)
{
  double one=1.0;
  double zero=0.0;
  // Calculate Lemma Matrix
  // V^-1 * A^(-1) * U
  // per walker: [k x N] * [N x k] = [k x k]
  callAndCheckError( cublasDgemmBatched( handle, CUBLAS_OP_N, CUBLAS_OP_N, k, k, N,
                                         &one,
                                         (const double**)AinvkList_d, RowStride,
                                         (const double**)U_d, RowStride, &zero,
                                         lemma_d, k,
                                         nw), __LINE__ );
  // Calculate - A^-1*dU
  dim3 dimBlockConvert (CONVERT_BS);
  dim3 dimGridConvert ((k*RowStride + (CONVERT_BS-1)) / CONVERT_BS, nw);
  // Calculate -dU=U(old)-U(new)
  subtract <<< dimGridConvert, dimBlockConvert >>> (AWorkList_d, AList_d, U_d, k*RowStride);
  // -A^(-1) * dU
  // per walker: [N x N] * [N x k] = [N x k]
  callAndCheckError( cublasDgemmBatched( handle, CUBLAS_OP_N, CUBLAS_OP_N, N, k, N,
                                         &one,
                                         (const double**)AinvList_d, RowStride,
                                         (const double**)AWorkList_d, RowStride, &zero,
                                         AinvUList_d, RowStride,
                                         nw), __LINE__ );
  // no synchronization needed here as this function call will always be followed by lemma lu calculation
}

void
cublas_smw_update (cublasHandle_t handle,
                   double *AinvkList_d[], double *AinvList_d[], 
                   double *AinvUList_d[], double *AWorkList_d[],
                   double *lemma_inv[], double *lemma_lu[],
                   int *infoArray,
                   int k, int kd, int M, int N, int nw, int RowStride)
{
#ifdef DEBUG_DELAYED
  fprintf(stderr,"*** Sherman-Morrison-Woodbury Update (k = %i, %i walkers) ***\n",k,nw);
#endif
  int pitch=RowStride;
  if(M==1) pitch=1;
  double one=1.0;

  // LU decomposition needs to be updated
  callAndCheckError( cublasDgetrfBatched( handle, k, lemma_lu, kd, NULL,
                                          infoArray, nw), __LINE__ );

#ifdef USE_TRSM
  if(M==1)
  {
    // {-A^-1 * dU } * Lemma^(-1) => solve for y: Lemma y * (L * U) = (y * L) * U = -A^-1 * dU
    // z * U = -A^-1 *dU
    callAndCheckError( cublasDtrsmBatched( handle, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT,
                                           M, k, &one, (const double**) lemma_lu, kd,
                                           AWorkList_d, pitch, nw), __LINE__ );
    // y * L = z => y = {-A^-1 * dU } * Lemma^(-1)
    callAndCheckError( cublasDtrsmBatched( handle, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_UNIT,
                                           M, k, &one, (const double**) lemma_lu, kd,
                                           AWorkList_d, pitch, nw), __LINE__ );
    // A^-1 + { -A^-1 * dU *  Lemma^-1 } * { V' * A^(-1) }
    // per walker: [1 x N] - [1 x k] * [k x N] = [1 x N]
    callAndCheckError( cublasDgemmBatched( handle, CUBLAS_OP_N, CUBLAS_OP_N, M, N, k,
                                           &one,
                                           (const double**)AWorkList_d, M,
                                           (const double**)AinvkList_d, RowStride, &one,
                                           AinvList_d, pitch,
                                           nw), __LINE__ );
  } else
  {
    // Lemma^(-1) * V' * A^(-1) => solve for y: Lemma (L * U) * y = L * (U * y) = V' * A^(-1)
    // L * z = V' * A^(-1)
    callAndCheckError( cublasDtrsmBatched( handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_UNIT,
                                           k, N, &one, (const double**) lemma_lu, kd,
                                           AWorkList_d, k, nw), __LINE__ );
    // U * y = z => y = Lemma^(-1) * V' * A^(-1)
    callAndCheckError( cublasDtrsmBatched( handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT,
                                           k, N, &one, (const double**) lemma_lu, kd,
                                           AWorkList_d, k, nw), __LINE__ );
    // A^-1 + { -A^-1 * dU } * { Lemma^-1 * V' * A^(-1) }
    // per walker: [M x N] - [M x k] * [k x N] = [M x N]
    callAndCheckError( cublasDgemmBatched( handle, CUBLAS_OP_N, CUBLAS_OP_N, M, N, k,
                                           &one,
                                           (const double**)AinvUList_d, RowStride,
                                           (const double**)AWorkList_d, k, &one,
                                           AinvList_d, pitch,
                                           nw), __LINE__ );
  }
#else
  double zero=0.0;
  // Calculate Lemma Inverse and store it in lemma_d
  // per walker: [k x k]^-1 = [k x k]
  callAndCheckError( cublasDgetriBatched( handle, k, (const double**) lemma_lu, kd, NULL,
                                          lemma_inv, k, infoArray, nw), __LINE__ );
  // Calculate new A inverse using Sherman-Morrison-Woodbury formula
  if(M==1) // row update can use different order to save flops
  {
    // { -A^-1 * dU } * Lemma^-1
    // per walker: [M x k] * [k x k] = [M x k]
    callAndCheckError( cublasDgemmBatched( handle, CUBLAS_OP_N, CUBLAS_OP_N, M, k, k,
                                           &one,
                                           (const double**)AinvUList_d, RowStride,
                                           (const double**)lemma_inv, k, &zero,
                                           AWorkList_d, M,
                                           nw), __LINE__ );
    // A^-1 - { A^-1 * dU  * Lemma^-1 } * { V' * A^(-1) }
    // per walker: [M x N] - [M x k] * [k x N] = [M x N]
    callAndCheckError( cublasDgemmBatched( handle, CUBLAS_OP_N, CUBLAS_OP_N, M, N, k,
                                           &one,
                                           (const double**)AWorkList_d, M,
                                           (const double**)AinvkList_d, RowStride, &one,
                                           AinvList_d, pitch,
                                           nw), __LINE__ );
  } else
  {
    // Need to use this matrix order for the overall update as AinvList and AinvkList have overlapping memory
    // Lemma^-1 * V' * A^(-1)
    // per walker: [k x k] * [k x N] = [k x N]
    callAndCheckError( cublasDgemmBatched( handle, CUBLAS_OP_N, CUBLAS_OP_N, k, N, k,
                                           &one,
                                           (const double**)lemma_inv, k,
                                           (const double**)AinvkList_d, RowStride, &zero,
                                           AWorkList_d, k,
                                           nw), __LINE__ );
    // A^-1 + { -A^-1 * dU } * { Lemma^-1 * V' * A^(-1) }
    // per walker: [M x N] - [M x k] * [k x N] = [M x N]
    callAndCheckError( cublasDgemmBatched( handle, CUBLAS_OP_N, CUBLAS_OP_N, M, N, k,
                                           &one,
                                           (const double**)AinvUList_d, RowStride,
                                           (const double**)AWorkList_d, k, &one,
                                           AinvList_d, pitch,
                                           nw), __LINE__ );
  }
#endif
}

// Four matrix inversion functions
// 1. for float matrices
//    useHigherPrecision = false --> single precision operations
//    useHigherPrecision = true  --> double precision operations  (default)
void
cublas_inverse (cublasHandle_t handle, 
                float *Alist_d[], float *Ainvlist_d[],
                float *AWorklist_d[], float *AinvWorklist_d[],
                int *PivotArray, int *infoArray,
                int N, int rowStride, int numMats,
                bool useHigherPrecision)
{

  // Info array tells if a matrix inversion is successful
  // = 0 : successful
  // = k : U(k,k) = 0; inversion failed 

  // If double precision operations are desired...
  if (useHigherPrecision)
  {

    // (i)   convert elements in Alist from float to double, put them in AWorklist
    dim3 dimBlockConvert (CONVERT_BS);
    dim3 dimGridConvert ((N*rowStride + (CONVERT_BS-1)) / CONVERT_BS, numMats);
    convert <<< dimGridConvert, dimBlockConvert >>> ((double**)AWorklist_d, Alist_d, N*rowStride);

    // (ii)  call cublas to do matrix inversion
    //       LU decomposition
    callAndCheckError( cublasDgetrfBatched( handle, N, (double**)AWorklist_d, rowStride, PivotArray,
                                            infoArray, numMats), __LINE__ );

    //       Inversion
#if (CUDA_VERSION >= 6050)
    callAndCheckError( cublasDgetriBatched( handle, N, (const double**)AWorklist_d, rowStride, PivotArray,
                                            (double**)AinvWorklist_d, rowStride, infoArray+numMats, numMats), __LINE__ );
#else
    callAndCheckError( cublasDgetriBatched( handle, N, (double**)AWorklist_d, rowStride, PivotArray,
                                            (double**)AinvWorklist_d, rowStride, infoArray+numMats, numMats), __LINE__ );
#endif

    // (iii) convert results back to single precision
    convert <<< dimGridConvert, dimBlockConvert >>> (Ainvlist_d, (double**)AinvWorklist_d, N*rowStride);

  }
  // else, carry out single precision operations
  else
  {
    // Call cublas to do matrix inversion
    // LU decomposition
    callAndCheckError( cublasSgetrfBatched( handle, N, Alist_d, rowStride, PivotArray,
                                            infoArray, numMats), __LINE__ );
  
    // Inversion
#if (CUDA_VERSION >= 6050)
    callAndCheckError( cublasSgetriBatched( handle, N, (const float**) Alist_d, rowStride, PivotArray,
                                            Ainvlist_d, rowStride, infoArray+numMats, numMats), __LINE__ );
#else
    callAndCheckError( cublasSgetriBatched( handle, N, Alist_d, rowStride, PivotArray,
                                            Ainvlist_d, rowStride, infoArray+numMats, numMats), __LINE__ );
#endif
  }

  cudaDeviceSynchronize();
}

// 2. for double matrices
void
cublas_inverse (cublasHandle_t handle, 
                double *Alist_d[], double *Ainvlist_d[],
                double *AWorklist_d[], double *AinvWorklist_d[],
                int *PivotArray, int *infoArray,
                int N, int rowStride, int numMats,
                bool useHigherPrecision)
{

  // Info array tells if a matrix inversion is successful
  // = 0 : successful
  // = k : U(k,k) = 0; inversion failed 

  // (i)   copy all the elements of Alist to AWorklist
  dim3 dimBlockConvert (CONVERT_BS);
  dim3 dimGridConvert ((N*rowStride + (CONVERT_BS-1)) / CONVERT_BS, numMats);
  convert <<< dimGridConvert, dimBlockConvert >>> (AWorklist_d, Alist_d, N*rowStride);
  
  // (ii)  call cublas functions to do inversion
  //       LU decomposition
  callAndCheckError( cublasDgetrfBatched( handle, N, AWorklist_d, rowStride, PivotArray,
                                          infoArray, numMats), __LINE__ );

  //       Inversion
#if (CUDA_VERSION >= 6050)
  callAndCheckError( cublasDgetriBatched( handle, N, (const double**) AWorklist_d, rowStride, PivotArray,
                                          Ainvlist_d, rowStride, infoArray+numMats, numMats), __LINE__ );
#else
  callAndCheckError( cublasDgetriBatched( handle, N, AWorklist_d, rowStride, PivotArray,
                                          Ainvlist_d, rowStride, infoArray+numMats, numMats), __LINE__ );
#endif

  cudaDeviceSynchronize();
}

// 3. for complex float matrices
//    useHigherPrecision = false --> single precision operations
//    useHigherPrecision = true  --> double precision operations  (default)
void
cublas_inverse (cublasHandle_t handle, 
                std::complex<float> *Alist_d[], std::complex<float> *Ainvlist_d[],
                std::complex<float> *AWorklist_d[], std::complex<float> *AinvWorklist_d[],
                int *PivotArray, int *infoArray,
                int N, int rowStride, int numMats,
                bool useHigherPrecision)
{

  // Info array tells if a matrix inversion is successful
  // = 0 : successful
  // = k : U(k,k) = 0; inversion failed 

  // If double precision operations are desired...
  if (useHigherPrecision)
  {

    // (i)   convert elements in Alist from float to double, put them in AWorklist
    dim3 dimBlockConvert (CONVERT_BS);
    dim3 dimGridConvert ((N*rowStride + (CONVERT_BS-1)) / CONVERT_BS, numMats);
    convert_complex<cuDoubleComplex, double, cuComplex> <<< dimGridConvert, dimBlockConvert >>> ((cuDoubleComplex**)AWorklist_d, (cuComplex**)Alist_d, N*rowStride);

    // (ii)  call cublas to do matrix inversion
    //       LU decomposition
    callAndCheckError( cublasZgetrfBatched( handle, N, (cuDoubleComplex**)AWorklist_d, rowStride, PivotArray, infoArray, numMats), __LINE__ );
    //       Inversion
#if (CUDA_VERSION >= 6050)
    callAndCheckError( cublasZgetriBatched( handle, N, (const cuDoubleComplex**)AWorklist_d, rowStride, PivotArray, (cuDoubleComplex**)AinvWorklist_d, rowStride, infoArray+numMats, numMats), __LINE__ );
#else
    callAndCheckError( cublasZgetriBatched( handle, N, (cuDoubleComplex**)AWorklist_d, rowStride, PivotArray, (cuDoubleComplex**)AinvWorklist_d, rowStride, infoArray+numMats, numMats), __LINE__ );
#endif

    // (iii) convert results back to single precision
    convert_complex<cuComplex, float, cuDoubleComplex> <<< dimGridConvert, dimBlockConvert >>> ((cuComplex**)Ainvlist_d, (cuDoubleComplex**)AinvWorklist_d, N*rowStride);

  }
  // else, carry out single precision operations
  else
  {
    // Call cublas to do matrix inversion
    // LU decomposition
    callAndCheckError( cublasCgetrfBatched( handle, N, (cuComplex**)Alist_d, rowStride, PivotArray,
                                            infoArray, numMats), __LINE__ );
  
    // Inversion
#if (CUDA_VERSION >= 6050)
    callAndCheckError( cublasCgetriBatched( handle, N, (const cuComplex**)Alist_d, rowStride, PivotArray, (cuComplex**)Ainvlist_d, rowStride, infoArray+numMats, numMats), __LINE__ );
#else
    callAndCheckError( cublasCgetriBatched( handle, N, (cuComplex**)Alist_d, rowStride, PivotArray, (cuComplex**)Ainvlist_d, rowStride, infoArray+numMats, numMats), __LINE__ );
#endif
  }

  cudaDeviceSynchronize();
}

// 4. for complex double matrices
void
cublas_inverse (cublasHandle_t handle, 
                std::complex<double> *Alist_d[], std::complex<double> *Ainvlist_d[],
                std::complex<double> *AWorklist_d[], std::complex<double> *AinvWorklist_d[],
                int *PivotArray, int *infoArray,
                int N, int rowStride, int numMats,
                bool useHigherPrecision)
{

  // Info array tells if a matrix inversion is successful
  // = 0 : successful
  // = k : U(k,k) = 0; inversion failed 

  // (i)   copy all the elements of Alist to AWorklist
  dim3 dimBlockConvert (CONVERT_BS);
  dim3 dimGridConvert ((N*rowStride + (CONVERT_BS-1)) / CONVERT_BS, numMats);
  convert_complex<cuDoubleComplex, double, cuDoubleComplex> <<< dimGridConvert, dimBlockConvert >>> ((cuDoubleComplex**)AWorklist_d, (cuDoubleComplex**)Alist_d, N*rowStride);

  // (ii)  call cublas to do matrix inversion
  //       LU decomposition
  callAndCheckError( cublasZgetrfBatched( handle, N, (cuDoubleComplex**)AWorklist_d, rowStride, PivotArray, infoArray, numMats), __LINE__ );
  //       Inversion
#if (CUDA_VERSION >= 6050)
  callAndCheckError( cublasZgetriBatched( handle, N, (const cuDoubleComplex**)AWorklist_d, rowStride, PivotArray, (cuDoubleComplex**)Ainvlist_d, rowStride, infoArray+numMats, numMats), __LINE__ );
#else
  callAndCheckError( cublasZgetriBatched( handle, N, (cuDoubleComplex**)AWorklist_d, rowStride, PivotArray, (cuDoubleComplex**)Ainvlist_d, rowStride, infoArray+numMats, numMats), __LINE__ );
#endif

  cudaDeviceSynchronize();
}



//////////////////////////////////////////////////////
//                  Test routines                   //
//////////////////////////////////////////////////////

#ifdef CUDA_TEST_MAIN

template<typename T>
void
test_cublas_inverse(int matSize, int numMats)
{

  // Initialize cublas
  cublasHandle_t handle;
  callAndCheckError( cublasCreate(&handle), __LINE__);

  srand48((long) 12394);
  int N = matSize;
  int row_stride = (matSize+15) / 16 * 16;
  T **Alist, **AWorklist;
  T **Alist_d, **AWorklist_d;
  T **Clist, **CWorklist;
  T **Clist_d, **CWorklist_d;

  // Allocate arrays of pointers (one set on host, one set on device)
  // pointing to the starting address (on device) of each matrix and its buffer
  // (similar to DiracDeterminantCUDA)
  Alist = (T**) malloc(numMats * sizeof(T*));
  callAndCheckError( cudaMalloc((void**) &Alist_d, numMats * sizeof(T*)), __LINE__ );

  AWorklist = (T**) malloc(numMats * sizeof(T*));
  callAndCheckError( cudaMalloc((void**) &AWorklist_d, numMats * sizeof(T*)), __LINE__ );

  Clist = (T**) malloc(numMats * sizeof(T*));
  callAndCheckError( cudaMalloc((void**) &Clist_d, numMats * sizeof(T*)), __LINE__ );

  CWorklist = (T**) malloc(numMats * sizeof(T*));
  callAndCheckError( cudaMalloc((void**) &CWorklist_d, numMats * sizeof(T*)), __LINE__ );

  // Generate matrices filled with random numbers
  T* A = (T*) malloc(sizeof(T) * numMats * N * row_stride);

  for (int j=0; j<numMats; j++)
    for (int i=0; i<N*row_stride; i++) {
#ifndef QMC_COMPLEX
        A[j*N*row_stride+i] = 1.0 * (drand48() - 0.5);
#else
        A[j*N*row_stride+i] = T(1.0 * (drand48() - 0.5), 1.0 * (drand48() - 0.5));
#endif
    }

  // Allocate memory on device for each matrix
  for (int mat=0; mat<numMats; mat++)
  {
    callAndCheckError( cudaMalloc((void**) &(Alist[mat]), N * row_stride * sizeof(T)), __LINE__ );

    callAndCheckError( cudaMemcpyAsync(Alist[mat], &A[mat*N*row_stride],
                                       N * row_stride * sizeof(T),
                                       cudaMemcpyHostToDevice), __LINE__ );

    callAndCheckError( cudaMalloc((void**) &(AWorklist[mat]), 2 * N * row_stride * sizeof(T)), __LINE__ );

    callAndCheckError( cudaMalloc((void**) &(Clist[mat]), N * row_stride * sizeof(T)), __LINE__ );

    callAndCheckError( cudaMalloc((void**) &(CWorklist[mat]), 2 * N * row_stride * sizeof(T)), __LINE__ );
  }

  // Copy the starting address of each matrix
  callAndCheckError( cudaMemcpyAsync (Alist_d, Alist, numMats * sizeof(T*),
                                      cudaMemcpyHostToDevice), __LINE__ ); 

  callAndCheckError( cudaMemcpyAsync (AWorklist_d, AWorklist, numMats * sizeof(T*),
                                      cudaMemcpyHostToDevice), __LINE__ );

  callAndCheckError( cudaMemcpyAsync (Clist_d, Clist, numMats * sizeof(T*),
                                      cudaMemcpyHostToDevice), __LINE__ );

  callAndCheckError( cudaMemcpyAsync (CWorklist_d, CWorklist, numMats * sizeof(T*),
                                      cudaMemcpyHostToDevice), __LINE__ );

  cudaDeviceSynchronize();
 
  clock_t start = clock();
 
  // Call cublas functions to do inversion
  cublas_inverse (handle, Alist_d, Clist_d, AWorklist_d, CWorklist_d, N, row_stride, numMats, true);
  cudaDeviceSynchronize();
  
  clock_t end = clock();
  double t = double(end-start) / double(CLOCKS_PER_SEC) / double(numMats);
  double rate = 1.0 / t;
  fprintf (stderr, "Rate is %1.3f matrix inversions per second.\n",
           rate);

  // Copy A^(-1) back to host memory Ainv; one matrix at a time
  // Calculate error of A^(-1)A from unit matrix I
  for (int mat=0; mat<numMats; mat++)
  {
    T Ainv[N*row_stride];
    callAndCheckError( cudaMemcpy(Ainv, Clist[mat], N * row_stride * sizeof(T),
                                  cudaMemcpyDeviceToHost), __LINE__ );

    double error = 0.0;
    for (int i=0; i<N; i++)
      for (int j=0; j<N; j++)
      {
        T val = 0.0;
        for (int k=0; k<N; k++)
          val += Ainv[i*row_stride+k] * A[mat*N*row_stride+k*row_stride+j];
        double diff = (i==j) ? (1.0f - std::real(val)) : std::real(val);
        error += diff * diff;
      }
    fprintf (stderr, "error = %1.8e\n", sqrt(error/(double)(N*N)));
  }

  // Finalize cublas
  callAndCheckError( cublasDestroy(handle), __LINE__ );

  // Free resources on both host and device
  for (int mat=0; mat<numMats; mat++)
  {
    cudaFree(Alist[mat]);
    cudaFree(Clist[mat]);
    cudaFree(AWorklist[mat]);
    cudaFree(CWorklist[mat]);
  }
  cudaFree(Alist_d);
  cudaFree(Clist_d);
  cudaFree(AWorklist_d);
  cudaFree(CWorklist_d);
  free(Alist);
  free(Clist);
  free(AWorklist);
  free(CWorklist);
  free(A);

  // Reset device. Required for memory leak debugging
  cudaDeviceReset();

}


int main(int argc, char** argv)
{
  int matSize = 0;
  int numMats = 0;

  if (argc == 3) {
    matSize = atoi(argv[1]);
    numMats = atoi(argv[2]);
  }
  else {
    printf("Usage: ./cuda_inverse [matrix size] [number of matrices]\n");
    exit(1);
  }

#ifndef QMC_COMPLEX
  test_cublas_inverse<double>(matSize, numMats);
  test_cublas_inverse<float>(matSize, numMats);
#else
  test_cublas_inverse<std::complex<double> >(matSize, numMats);
  test_cublas_inverse<std::complex<float> >(matSize, numMats);
#endif

  return 0;
}

#endif
