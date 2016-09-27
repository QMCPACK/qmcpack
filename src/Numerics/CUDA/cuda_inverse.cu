//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, StoneRidge Inc.
//		      Jeremy McMinnis, jmcminis@gmail.com, Navar Inc. 
//    		      Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//                    Ying Wai Li, yingwaili@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by:  Ken Esler, kpesler@gmail.com, StoneRidge Inc.
//////////////////////////////////////////////////////////////////////////////////////




// ============ Matrix inversion using cuBLAS library ============ //
//
// To compile as a standalone test program:
//
// 1. Make sure libcublas.so is in the search path
// 2. cd to build/ directory
// 3. nvcc -o cuda_inverse -arch=sm_35 -lcublas -DCUDA_TEST_MAIN
//         ../src/Numerics/CUDA/cuda_inverse.cu
//
// =============================================================== //

#include <cstdio>
#include <unistd.h>
#include <sstream>
#include <vector>
#include <iostream>
#include <cuda.h>
#include <cublas_v2.h>

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

// Two matrix inversion functions
// 1. for float matrices
//    useHigherPrecision = false --> single precision operations
//    useHigherPrecision = true  --> double precision operations
void
cublas_inverse (cublasHandle_t handle, 
                float *Alist_d[], float *Ainvlist_d[],
                float *AWorklist_d[], float *AinvWorklist_d[],
                int N, int rowStride, int numMats,
                bool useHigherPrecision)
{

  // Info array tells if a matrix inversion is successful
  // = 0 : successful
  // = k : U(k,k) = 0; inversion failed 
  int *infoArray;
  callAndCheckError( cudaMalloc((void**) &infoArray, numMats * sizeof(int)), __LINE__ ); 

  // If double precision operations are desired...
  if (useHigherPrecision)
  {

    // (i)   convert elements in Alist from float to double, put them in AWorklist
    dim3 dimBlockConvert (CONVERT_BS);
    dim3 dimGridConvert ((N*rowStride + (CONVERT_BS-1)) / CONVERT_BS, numMats);
    convert <<< dimGridConvert, dimBlockConvert >>> ((double**)AWorklist_d, Alist_d, N*rowStride);

    // (ii)  call cublas to do matrix inversion
    //       LU decomposition
    callAndCheckError( cublasDgetrfBatched( handle, N, (double**)AWorklist_d, rowStride, NULL, 
                                            infoArray, numMats), __LINE__ );
  
    //       Inversion
#if (CUDA_VERSION >= 6050)
    callAndCheckError( cublasDgetriBatched( handle, N, (const double**)AWorklist_d, rowStride, NULL, 
                                            (double**)AinvWorklist_d, rowStride, infoArray, numMats), __LINE__ );
#else
    callAndCheckError( cublasDgetriBatched( handle, N, (double**)AWorklist_d, rowStride, NULL, 
                                            (double**)AinvWorklist_d, rowStride, infoArray, numMats), __LINE__ );
#endif

    // (iii) convert results back to single precision
    convert <<< dimGridConvert, dimBlockConvert >>> (Ainvlist_d, (double**)AinvWorklist_d, N*rowStride);

  }
  // else, carry out single precision operations
  else
  {
    // Call cublas to do matrix inversion
    // LU decomposition
    callAndCheckError( cublasSgetrfBatched( handle, N, Alist_d, rowStride, NULL, 
                                            infoArray, numMats), __LINE__ );
  
    // Inversion
#if (CUDA_VERSION >= 6050)
    callAndCheckError( cublasSgetriBatched( handle, N, (const float**) Alist_d, rowStride, NULL, 
                                            Ainvlist_d, rowStride, infoArray, numMats), __LINE__ );
#else
    callAndCheckError( cublasSgetriBatched( handle, N, Alist_d, rowStride, NULL, 
                                            Ainvlist_d, rowStride, infoArray, numMats), __LINE__ );
#endif
  }

  cudaDeviceSynchronize();

  // Free resources
  cudaFree(infoArray);

}

// 2. for double matrices
void
cublas_inverse (cublasHandle_t handle, 
                double *Alist_d[], double *Ainvlist_d[],
                double *AWorklist_d[], double *AinvWorklist_d[],
                int N, int rowStride, int numMats,
                bool useHigherPrecision)
{

  // Info array tells if a matrix inversion is successful
  // = 0 : successful
  // = k : U(k,k) = 0; inversion failed 
  int *infoArray;
  callAndCheckError( cudaMalloc((void**) &infoArray, numMats * sizeof(int)), __LINE__ );

  // (i)   copy all the elements of Alist to AWorklist
  dim3 dimBlockConvert (CONVERT_BS);
  dim3 dimGridConvert ((N*rowStride + (CONVERT_BS-1)) / CONVERT_BS, numMats);
  convert <<< dimGridConvert, dimBlockConvert >>> (AWorklist_d, Alist_d, N*rowStride);
  
  // (ii)  call cublas functions to do inversion
  //       LU decomposition
  callAndCheckError( cublasDgetrfBatched( handle, N, AWorklist_d, rowStride, NULL, 
                                          infoArray, numMats), __LINE__ );

  //       Inversion
#if (CUDA_VERSION >= 6050)
  callAndCheckError( cublasDgetriBatched( handle, N, (const double**) AWorklist_d, rowStride, NULL, 
                                          Ainvlist_d, rowStride, infoArray, numMats), __LINE__ );
#else
  callAndCheckError( cublasDgetriBatched( handle, N, AWorklist_d, rowStride, NULL, 
                                          Ainvlist_d, rowStride, infoArray, numMats), __LINE__ );
#endif

  cudaDeviceSynchronize();

  cudaFree(infoArray);

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
    for (int i=0; i<N*row_stride; i++)
      A[j*N*row_stride+i] = 1.0 * (drand48() - 0.5);

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
        double val = 0.0;
        for (int k=0; k<N; k++)
          val += Ainv[i*row_stride+k] * A[mat*N*row_stride+k*row_stride+j];
        double diff = (i==j) ? (1.0f - val) : val;
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

  test_cublas_inverse<double>(matSize, numMats);
  test_cublas_inverse<float>(matSize, numMats);

  return 0;
}

#endif
