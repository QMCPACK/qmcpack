#include <cstdio>
#include <unistd.h>
#include <sstream>
#include <vector>
#include <iostream>
#include <complex>
#include <cuda.h>
#include <cublas_v2.h>
#include <cuComplex.h>
#include <thrust/complex.h>

#define BLOCKSIZE 256

// #define DEBUG_DELAYED
// #define USE_TRSM

template<typename T>
__global__ void finish_lemma_calc(T** ainvu, T** lemma, int k, int kstart, int N, int rowstride)
{
  __shared__ T *my_ainvu, *my_lemma;
  if (threadIdx.x == 0)
  {
    my_lemma = lemma[blockIdx.y];
    my_ainvu = ainvu[blockIdx.y] + kstart;
  }
  __syncthreads();
  int i = blockIdx.x * BLOCKSIZE + threadIdx.x;
  int r = i / k;
  int s = i % k;
  if (i < k * k)
  {
    my_lemma[r * k + s] = -my_ainvu[r * rowstride + s];
    // A^-1*A_k is only 1.0 for index i*rowstride+i+kstart per column and can be added now that the value has been read
    if (r == s)
      my_ainvu[r * rowstride + r] += 1.0;
  }
}

/** Calculate Lemma Matrix: I_k + V' * ( A^(-1) * U )
  * for each walker
  * -> returns L-U decomposed lemma matrix for easy determinant calculations and inverse calculation later
  */
void cublas_lemma_mats(cublasHandle_t handle,
                       float* AinvList_d[],
                       float* U_d[],
                       float* lemma_d[],
                       float* AinvUList_d[],
                       int k,
                       int kstart,
                       int N,
                       int nw,
                       int RowStride)
{
  float mone = -1.0;
  float zero = 0.0;
  // -A^(-1) * U
  // per walker: [N x N] * [N x k] = [N x k]
  cublasSgemmBatched(handle,
                     CUBLAS_OP_N,
                     CUBLAS_OP_N,
                     N,
                     k,
                     N,
                     &mone,
                     (const float**)AinvList_d,
                     RowStride,
                     (const float**)U_d,
                     RowStride,
                     &zero,
                     AinvUList_d,
                     RowStride,
                     nw);
  // Calculate Lemma Matrix using above calculations and calculate -A^-1*dU
  dim3 dimBlockConvert(BLOCKSIZE);
  dim3 dimGridConvert((k * k + (BLOCKSIZE - 1)) / BLOCKSIZE, nw);
  finish_lemma_calc<<<dimGridConvert, dimBlockConvert>>>(AinvUList_d, lemma_d, k, kstart, N, RowStride);
}

void cublas_ainv_row(cublasHandle_t handle,
                     float* AinvkList_d[],
                     float* AWorkList_d[],
                     float* AinvList_d[],
                     int k,
                     int N,
                     int nw,
                     int RowStride)
{
  float one = 1.0;
  // A^-1 - { A^-1 * dU  * Lemma^-1 } * { V' * A^(-1) }
  // per walker: [1 x N] - [1 x k] * [k x N] = [1 x N]
  cublasSgemmBatched(handle,
                     CUBLAS_OP_N,
                     CUBLAS_OP_N,
                     1,
                     N,
                     k,
                     &one,
                     (const float**)AWorkList_d,
                     1,
                     (const float**)AinvkList_d,
                     RowStride,
                     &one,
                     AinvList_d,
                     1,
                     nw);
}

void cublas_ainv_row(cublasHandle_t handle,
                     double* AinvkList_d[],
                     double* AWorkList_d[],
                     double* AinvList_d[],
                     int k,
                     int N,
                     int nw,
                     int RowStride)
{
  double one = 1.0;
  // A^-1 - { A^-1 * dU  * Lemma^-1 } * { V' * A^(-1) }
  // per walker: [1 x N] - [1 x k] * [k x N] = [1 x N]
  cublasDgemmBatched(handle,
                     CUBLAS_OP_N,
                     CUBLAS_OP_N,
                     1,
                     N,
                     k,
                     &one,
                     (const double**)AWorkList_d,
                     1,
                     (const double**)AinvkList_d,
                     RowStride,
                     &one,
                     AinvList_d,
                     1,
                     nw);
}


void cublas_smw_update(cublasHandle_t handle,
                       float* AinvkList_d[],
                       float* AinvList_d[],
                       float* AinvUList_d[],
                       float* AWorkList_d[],
                       float* lemma_inv[],
                       float* lemma_lu[],
                       int* PivotArray,
                       int* infoArray,
                       int k,
                       int kd,
                       int M,
                       int N,
                       int nw,
                       int RowStride)
{
#ifdef DEBUG_DELAYED
  fprintf(stderr, "*** Sherman-Morrison-Woodbury Update (k = %i, %i walkers) ***\n", k, nw);
#endif
  int pitch = RowStride;
  if (M == 1)
    pitch = 1;
  float one = 1.0;

  // LU decomposition needs to be updated
  cublasSgetrfBatched(handle, k, lemma_lu, kd, PivotArray, infoArray, nw);

#ifdef USE_TRSM
  if (M == 1)
  {
    // {-A^-1 * dU } * Lemma^(-1) => solve for y: Lemma y * (L * U) = (y * L) * U = -A^-1 * dU
    // z * U = -A^-1 *dU
    cublasStrsmBatched(handle,
                       CUBLAS_SIDE_RIGHT,
                       CUBLAS_FILL_MODE_UPPER,
                       CUBLAS_OP_N,
                       CUBLAS_DIAG_NON_UNIT,
                       M,
                       k,
                       &one,
                       (const float**)lemma_lu,
                       kd,
                       AWorkList_d,
                       pitch,
                       nw);
    // y * L = z => y = {-A^-1 * dU } * Lemma^(-1)
    cublasStrsmBatched(handle,
                       CUBLAS_SIDE_RIGHT,
                       CUBLAS_FILL_MODE_LOWER,
                       CUBLAS_OP_N,
                       CUBLAS_DIAG_UNIT,
                       M,
                       k,
                       &one,
                       (const float**)lemma_lu,
                       kd,
                       AWorkList_d,
                       pitch,
                       nw);
    // A^-1 + { -A^-1 * dU *  Lemma^-1 } * { V' * A^(-1) }
    // per walker: [1 x N] - [1 x k] * [k x N] = [1 x N]
    cublasSgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       M,
                       N,
                       k,
                       &one,
                       (const float**)AWorkList_d,
                       M,
                       (const float**)AinvkList_d,
                       RowStride,
                       &one,
                       AinvList_d,
                       pitch,
                       nw);
  }
  else
  {
    // Lemma^(-1) * V' * A^(-1) => solve for y: Lemma (L * U) * y = L * (U * y) = V' * A^(-1)
    // L * z = V' * A^(-1)
    cublasStrsmBatched(handle,
                       CUBLAS_SIDE_LEFT,
                       CUBLAS_FILL_MODE_LOWER,
                       CUBLAS_OP_N,
                       CUBLAS_DIAG_UNIT,
                       k,
                       N,
                       &one,
                       (const float**)lemma_lu,
                       kd,
                       AWorkList_d,
                       k,
                       nw);
    // U * y = z => y = Lemma^(-1) * V' * A^(-1)
    cublasStrsmBatched(handle,
                       CUBLAS_SIDE_LEFT,
                       CUBLAS_FILL_MODE_UPPER,
                       CUBLAS_OP_N,
                       CUBLAS_DIAG_NON_UNIT,
                       k,
                       N,
                       &one,
                       (const float**)lemma_lu,
                       kd,
                       AWorkList_d,
                       k,
                       nw);
    // A^-1 + { -A^-1 * dU } * { Lemma^-1 * V' * A^(-1) }
    // per walker: [M x N] - [M x k] * [k x N] = [M x N]
    cublasSgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       M,
                       N,
                       k,
                       &one,
                       (const float**)AinvUList_d,
                       RowStride,
                       (const float**)AWorkList_d,
                       k,
                       &one,
                       AinvList_d,
                       pitch,
                       nw);
  }
#else
  float zero = 0.0;
  // Calculate Lemma Inverse and store it in lemma_d
  // per walker: [k x k]^-1 = [k x k]
  cublasSgetriBatched(handle, k, (const float**)lemma_lu, kd, PivotArray, lemma_inv, k, infoArray, nw);
  // Calculate new A inverse using Sherman-Morrison-Woodbury formula
  if (M == 1) // row update can use different order to save flops
  {
    // { -A^-1 * dU } * Lemma^-1
    // per walker: [M x k] * [k x k] = [M x k]
    cublasSgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       M,
                       k,
                       k,
                       &one,
                       (const float**)AinvUList_d,
                       RowStride,
                       (const float**)lemma_inv,
                       k,
                       &zero,
                       AWorkList_d,
                       M,
                       nw);
    // A^-1 - { A^-1 * dU  * Lemma^-1 } * { V' * A^(-1) }
    // per walker: [M x N] - [M x k] * [k x N] = [M x N]
    cublasSgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       M,
                       N,
                       k,
                       &one,
                       (const float**)AWorkList_d,
                       M,
                       (const float**)AinvkList_d,
                       RowStride,
                       &one,
                       AinvList_d,
                       pitch,
                       nw);
  }
  else
  {
    // Need to use this matrix order for the overall update as AinvList and AinvkList have overlapping memory
    // Lemma^-1 * V' * A^(-1)
    // per walker: [k x k] * [k x N] = [k x N]
    cublasSgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       k,
                       N,
                       k,
                       &one,
                       (const float**)lemma_inv,
                       k,
                       (const float**)AinvkList_d,
                       RowStride,
                       &zero,
                       AWorkList_d,
                       k,
                       nw);
    // A^-1 + { -A^-1 * dU } * { Lemma^-1 * V' * A^(-1) }
    // per walker: [M x N] - [M x k] * [k x N] = [M x N]
    cublasSgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       M,
                       N,
                       k,
                       &one,
                       (const float**)AinvUList_d,
                       RowStride,
                       (const float**)AWorkList_d,
                       k,
                       &one,
                       AinvList_d,
                       pitch,
                       nw);
  }
#endif
}

/** Calculate Lemma Matrix: I_k + V' * ( A^(-1) * U )
  * for each walker
  * -> returns L-U decomposed lemma matrix for easy determinant calculations and inverse calculation later
  */
void cublas_lemma_mats(cublasHandle_t handle,
                       double* AinvList_d[],
                       double* U_d[],
                       double* lemma_d[],
                       double* AinvUList_d[],
                       int k,
                       int kstart,
                       int N,
                       int nw,
                       int RowStride)
{
  double mone = -1.0;
  double zero = 0.0;
  // -A^(-1) * U
  // per walker: [N x N] * [N x k] = [N x k]
  cublasDgemmBatched(handle,
                     CUBLAS_OP_N,
                     CUBLAS_OP_N,
                     N,
                     k,
                     N,
                     &mone,
                     (const double**)AinvList_d,
                     RowStride,
                     (const double**)U_d,
                     RowStride,
                     &zero,
                     AinvUList_d,
                     RowStride,
                     nw);
  // Calculate Lemma Matrix using above calculations and calculate -A^-1*dU
  dim3 dimBlockConvert(BLOCKSIZE);
  dim3 dimGridConvert((k * k + (BLOCKSIZE - 1)) / BLOCKSIZE, nw);
  finish_lemma_calc<<<dimGridConvert, dimBlockConvert>>>(AinvUList_d, lemma_d, k, kstart, N, RowStride);
}

void cublas_smw_update(cublasHandle_t handle,
                       double* AinvkList_d[],
                       double* AinvList_d[],
                       double* AinvUList_d[],
                       double* AWorkList_d[],
                       double* lemma_inv[],
                       double* lemma_lu[],
                       int* PivotArray,
                       int* infoArray,
                       int k,
                       int kd,
                       int M,
                       int N,
                       int nw,
                       int RowStride)
{
#ifdef DEBUG_DELAYED
  fprintf(stderr, "*** Sherman-Morrison-Woodbury Update (k = %i, %i walkers) ***\n", k, nw);
#endif
  int pitch = RowStride;
  if (M == 1)
    pitch = 1;
  double one = 1.0;

  // LU decomposition needs to be updated
  cublasDgetrfBatched(handle, k, lemma_lu, kd, PivotArray, infoArray, nw);

#ifdef USE_TRSM
  if (M == 1)
  {
    // {-A^-1 * dU } * Lemma^(-1) => solve for y: Lemma y * (L * U) = (y * L) * U = -A^-1 * dU
    // z * U = -A^-1 *dU
    cublasDtrsmBatched(handle,
                       CUBLAS_SIDE_RIGHT,
                       CUBLAS_FILL_MODE_UPPER,
                       CUBLAS_OP_N,
                       CUBLAS_DIAG_NON_UNIT,
                       M,
                       k,
                       &one,
                       (const double**)lemma_lu,
                       kd,
                       AWorkList_d,
                       pitch,
                       nw);
    // y * L = z => y = {-A^-1 * dU } * Lemma^(-1)
    cublasDtrsmBatched(handle,
                       CUBLAS_SIDE_RIGHT,
                       CUBLAS_FILL_MODE_LOWER,
                       CUBLAS_OP_N,
                       CUBLAS_DIAG_UNIT,
                       M,
                       k,
                       &one,
                       (const double**)lemma_lu,
                       kd,
                       AWorkList_d,
                       pitch,
                       nw);
    // A^-1 + { -A^-1 * dU *  Lemma^-1 } * { V' * A^(-1) }
    // per walker: [1 x N] - [1 x k] * [k x N] = [1 x N]
    cublasDgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       M,
                       N,
                       k,
                       &one,
                       (const double**)AWorkList_d,
                       M,
                       (const double**)AinvkList_d,
                       RowStride,
                       &one,
                       AinvList_d,
                       pitch,
                       nw);
  }
  else
  {
    // Lemma^(-1) * V' * A^(-1) => solve for y: Lemma (L * U) * y = L * (U * y) = V' * A^(-1)
    // L * z = V' * A^(-1)
    cublasDtrsmBatched(handle,
                       CUBLAS_SIDE_LEFT,
                       CUBLAS_FILL_MODE_LOWER,
                       CUBLAS_OP_N,
                       CUBLAS_DIAG_UNIT,
                       k,
                       N,
                       &one,
                       (const double**)lemma_lu,
                       kd,
                       AWorkList_d,
                       k,
                       nw);
    // U * y = z => y = Lemma^(-1) * V' * A^(-1)
    cublasDtrsmBatched(handle,
                       CUBLAS_SIDE_LEFT,
                       CUBLAS_FILL_MODE_UPPER,
                       CUBLAS_OP_N,
                       CUBLAS_DIAG_NON_UNIT,
                       k,
                       N,
                       &one,
                       (const double**)lemma_lu,
                       kd,
                       AWorkList_d,
                       k,
                       nw);
    // A^-1 + { -A^-1 * dU } * { Lemma^-1 * V' * A^(-1) }
    // per walker: [M x N] - [M x k] * [k x N] = [M x N]
    cublasDgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       M,
                       N,
                       k,
                       &one,
                       (const double**)AinvUList_d,
                       RowStride,
                       (const double**)AWorkList_d,
                       k,
                       &one,
                       AinvList_d,
                       pitch,
                       nw);
  }
#else
  double zero = 0.0;
  // Calculate Lemma Inverse and store it in lemma_d
  // per walker: [k x k]^-1 = [k x k]
  cublasDgetriBatched(handle, k, (const double**)lemma_lu, kd, PivotArray, lemma_inv, k, infoArray, nw);
  // Calculate new A inverse using Sherman-Morrison-Woodbury formula
  if (M == 1) // row update can use different order to save flops
  {
    // { -A^-1 * dU } * Lemma^-1
    // per walker: [M x k] * [k x k] = [M x k]
    cublasDgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       M,
                       k,
                       k,
                       &one,
                       (const double**)AinvUList_d,
                       RowStride,
                       (const double**)lemma_inv,
                       k,
                       &zero,
                       AWorkList_d,
                       M,
                       nw);
    // A^-1 - { A^-1 * dU  * Lemma^-1 } * { V' * A^(-1) }
    // per walker: [M x N] - [M x k] * [k x N] = [M x N]
    cublasDgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       M,
                       N,
                       k,
                       &one,
                       (const double**)AWorkList_d,
                       M,
                       (const double**)AinvkList_d,
                       RowStride,
                       &one,
                       AinvList_d,
                       pitch,
                       nw);
  }
  else
  {
    // Need to use this matrix order for the overall update as AinvList and AinvkList have overlapping memory
    // Lemma^-1 * V' * A^(-1)
    // per walker: [k x k] * [k x N] = [k x N]
    cublasDgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       k,
                       N,
                       k,
                       &one,
                       (const double**)lemma_inv,
                       k,
                       (const double**)AinvkList_d,
                       RowStride,
                       &zero,
                       AWorkList_d,
                       k,
                       nw);
    // A^-1 + { -A^-1 * dU } * { Lemma^-1 * V' * A^(-1) }
    // per walker: [M x N] - [M x k] * [k x N] = [M x N]
    cublasDgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       M,
                       N,
                       k,
                       &one,
                       (const double**)AinvUList_d,
                       RowStride,
                       (const double**)AWorkList_d,
                       k,
                       &one,
                       AinvList_d,
                       pitch,
                       nw);
  }
#endif
}

/** Calculate Lemma Matrix: I_k + V' * ( A^(-1) * U )
  * for each walker
  * -> returns L-U decomposed lemma matrix for easy determinant calculations and inverse calculation later
  */
void cublas_lemma_mats(cublasHandle_t handle,
                       std::complex<float>* AinvList_d[],
                       std::complex<float>* U_d[],
                       std::complex<float>* lemma_d[],
                       std::complex<float>* AinvUList_d[],
                       int k,
                       int kstart,
                       int N,
                       int nw,
                       int RowStride)
{
  cuComplex mone = make_cuComplex(-1.0f, 0.0f);
  cuComplex zero = make_cuComplex(0.0f, 0.0f);
  // -A^(-1) * U
  // per walker: [N x N] * [N x k] = [N x k]
  cublasCgemmBatched(handle,
                     CUBLAS_OP_N,
                     CUBLAS_OP_N,
                     N,
                     k,
                     N,
                     &mone,
                     (const cuComplex**)AinvList_d,
                     RowStride,
                     (const cuComplex**)U_d,
                     RowStride,
                     &zero,
                     (cuComplex**)AinvUList_d,
                     RowStride,
                     nw);
  // Calculate Lemma Matrix using above calculations and calculate -A^-1*dU
  dim3 dimBlockConvert(BLOCKSIZE);
  dim3 dimGridConvert((k * k + (BLOCKSIZE - 1)) / BLOCKSIZE, nw);
  finish_lemma_calc<<<dimGridConvert, dimBlockConvert>>>((thrust::complex<float>**)AinvUList_d,
                                                         (thrust::complex<float>**)lemma_d,
                                                         k,
                                                         kstart,
                                                         N,
                                                         RowStride);
}

void cublas_ainv_row(cublasHandle_t handle,
                     std::complex<float>* AinvkList_d[],
                     std::complex<float>* AWorkList_d[],
                     std::complex<float>* AinvList_d[],
                     int k,
                     int N,
                     int nw,
                     int RowStride)
{
  cuComplex one = make_cuComplex(1.0f, 0.0f);
  // A^-1 - { A^-1 * dU  * Lemma^-1 } * { V' * A^(-1) }
  // per walker: [1 x N] - [1 x k] * [k x N] = [1 x N]
  cublasCgemmBatched(handle,
                     CUBLAS_OP_N,
                     CUBLAS_OP_N,
                     1,
                     N,
                     k,
                     &one,
                     (const cuComplex**)AWorkList_d,
                     1,
                     (const cuComplex**)AinvkList_d,
                     RowStride,
                     &one,
                     (cuComplex**)AinvList_d,
                     1,
                     nw);
}

void cublas_smw_update(cublasHandle_t handle,
                       std::complex<float>* AinvkList_d[],
                       std::complex<float>* AinvList_d[],
                       std::complex<float>* AinvUList_d[],
                       std::complex<float>* AWorkList_d[],
                       std::complex<float>* lemma_inv[],
                       std::complex<float>* lemma_lu[],
                       int* PivotArray,
                       int* infoArray,
                       int k,
                       int kd,
                       int M,
                       int N,
                       int nw,
                       int RowStride)
{
#ifdef DEBUG_DELAYED
  fprintf(stderr, "*** Sherman-Morrison-Woodbury Update (k = %i, %i walkers) ***\n", k, nw);
#endif
  int pitch = RowStride;
  if (M == 1)
    pitch = 1;
  cuComplex one = make_cuComplex(1.0f, 0.0f);

  // LU decomposition needs to be updated
  cublasCgetrfBatched(handle, k, (cuComplex**)lemma_lu, kd, PivotArray, infoArray, nw);

#ifdef USE_TRSM
  if (M == 1)
  {
    // {-A^-1 * dU } * Lemma^(-1) => solve for y: Lemma y * (L * U) = (y * L) * U = -A^-1 * dU
    // z * U = -A^-1 *dU
    cublasCtrsmBatched(handle,
                       CUBLAS_SIDE_RIGHT,
                       CUBLAS_FILL_MODE_UPPER,
                       CUBLAS_OP_N,
                       CUBLAS_DIAG_NON_UNIT,
                       M,
                       k,
                       &one,
                       (const cuComplex**)lemma_lu,
                       kd,
                       (cuComplex**)AWorkList_d,
                       pitch,
                       nw);
    // y * L = z => y = {-A^-1 * dU } * Lemma^(-1)
    cublasCtrsmBatched(handle,
                       CUBLAS_SIDE_RIGHT,
                       CUBLAS_FILL_MODE_LOWER,
                       CUBLAS_OP_N,
                       CUBLAS_DIAG_UNIT,
                       M,
                       k,
                       &one,
                       (const cuComplex**)lemma_lu,
                       kd,
                       (cuComplex**)AWorkList_d,
                       pitch,
                       nw);
    // A^-1 + { -A^-1 * dU *  Lemma^-1 } * { V' * A^(-1) }
    // per walker: [1 x N] - [1 x k] * [k x N] = [1 x N]
    cublasCgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       M,
                       N,
                       k,
                       &one,
                       (const cuComplex**)AWorkList_d,
                       M,
                       (const cuComplex**)AinvkList_d,
                       RowStride,
                       &one,
                       (cuComplex**)AinvList_d,
                       pitch,
                       nw);
  }
  else
  {
    // Lemma^(-1) * V' * A^(-1) => solve for y: Lemma (L * U) * y = L * (U * y) = V' * A^(-1)
    // L * z = V' * A^(-1)
    cublasCtrsmBatched(handle,
                       CUBLAS_SIDE_LEFT,
                       CUBLAS_FILL_MODE_LOWER,
                       CUBLAS_OP_N,
                       CUBLAS_DIAG_UNIT,
                       k,
                       N,
                       &one,
                       (const cuComplex**)lemma_lu,
                       kd,
                       (cuComplex**)AWorkList_d,
                       k,
                       nw);
    // U * y = z => y = Lemma^(-1) * V' * A^(-1)
    cublasCtrsmBatched(handle,
                       CUBLAS_SIDE_LEFT,
                       CUBLAS_FILL_MODE_UPPER,
                       CUBLAS_OP_N,
                       CUBLAS_DIAG_NON_UNIT,
                       k,
                       N,
                       &one,
                       (const cuComplex**)lemma_lu,
                       kd,
                       (cuComplex**)AWorkList_d,
                       k,
                       nw);
    // A^-1 + { -A^-1 * dU } * { Lemma^-1 * V' * A^(-1) }
    // per walker: [M x N] - [M x k] * [k x N] = [M x N]
    cublasCgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       M,
                       N,
                       k,
                       &one,
                       (const cuComplex**)AinvUList_d,
                       RowStride,
                       (const cuComplex**)AWorkList_d,
                       k,
                       &one,
                       (cuComplex**)AinvList_d,
                       pitch,
                       nw);
  }
#else
  cuComplex zero = make_cuComplex(0.0f, 0.0f);
  // Calculate Lemma Inverse and store it in lemma_d
  // per walker: [k x k]^-1 = [k x k]
  cublasCgetriBatched(handle, k, (const cuComplex**)lemma_lu, kd, PivotArray, (cuComplex**)lemma_inv, k, infoArray, nw);
  // Calculate new A inverse using Sherman-Morrison-Woodbury formula
  if (M == 1) // row update can use different order to save flops
  {
    // { -A^-1 * dU } * Lemma^-1
    // per walker: [M x k] * [k x k] = [M x k]
    cublasCgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       M,
                       k,
                       k,
                       &one,
                       (const cuComplex**)AinvUList_d,
                       RowStride,
                       (const cuComplex**)lemma_inv,
                       k,
                       &zero,
                       (cuComplex**)AWorkList_d,
                       M,
                       nw);
    // A^-1 - { A^-1 * dU  * Lemma^-1 } * { V' * A^(-1) }
    // per walker: [M x N] - [M x k] * [k x N] = [M x N]
    cublasCgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       M,
                       N,
                       k,
                       &one,
                       (const cuComplex**)AWorkList_d,
                       M,
                       (const cuComplex**)AinvkList_d,
                       RowStride,
                       &one,
                       (cuComplex**)AinvList_d,
                       pitch,
                       nw);
  }
  else
  {
    // Need to use this matrix order for the overall update as AinvList and AinvkList have overlapping memory
    // Lemma^-1 * V' * A^(-1)
    // per walker: [k x k] * [k x N] = [k x N]
    cublasCgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       k,
                       N,
                       k,
                       &one,
                       (const cuComplex**)lemma_inv,
                       k,
                       (const cuComplex**)AinvkList_d,
                       RowStride,
                       &zero,
                       (cuComplex**)AWorkList_d,
                       k,
                       nw);
    // A^-1 + { -A^-1 * dU } * { Lemma^-1 * V' * A^(-1) }
    // per walker: [M x N] - [M x k] * [k x N] = [M x N]
    cublasCgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       M,
                       N,
                       k,
                       &one,
                       (const cuComplex**)AinvUList_d,
                       RowStride,
                       (const cuComplex**)AWorkList_d,
                       k,
                       &one,
                       (cuComplex**)AinvList_d,
                       pitch,
                       nw);
  }
#endif
}


/** Calculate Lemma Matrix: I_k + V' * ( A^(-1) * U )
  * for each walker
  * -> returns L-U decomposed lemma matrix for easy determinant calculations and inverse calculation later
  */
void cublas_lemma_mats(cublasHandle_t handle,
                       std::complex<double>* AinvList_d[],
                       std::complex<double>* U_d[],
                       std::complex<double>* lemma_d[],
                       std::complex<double>* AinvUList_d[],
                       int k,
                       int kstart,
                       int N,
                       int nw,
                       int RowStride)
{
  cuDoubleComplex mone = make_cuDoubleComplex(-1.0f, 0.0f);
  cuDoubleComplex zero = make_cuDoubleComplex(0.0f, 0.0f);
  // -A^(-1) * U
  // per walker: [N x N] * [N x k] = [N x k]
  cublasZgemmBatched(handle,
                     CUBLAS_OP_N,
                     CUBLAS_OP_N,
                     N,
                     k,
                     N,
                     &mone,
                     (const cuDoubleComplex**)AinvList_d,
                     RowStride,
                     (const cuDoubleComplex**)U_d,
                     RowStride,
                     &zero,
                     (cuDoubleComplex**)AinvUList_d,
                     RowStride,
                     nw);
  // Calculate Lemma Matrix using above calculations and calculate -A^-1*dU
  dim3 dimBlockConvert(BLOCKSIZE);
  dim3 dimGridConvert((k * k + (BLOCKSIZE - 1)) / BLOCKSIZE, nw);
  finish_lemma_calc<<<dimGridConvert, dimBlockConvert>>>((thrust::complex<double>**)AinvUList_d,
                                                         (thrust::complex<double>**)lemma_d,
                                                         k,
                                                         kstart,
                                                         N,
                                                         RowStride);
}

void cublas_ainv_row(cublasHandle_t handle,
                     std::complex<double>* AinvkList_d[],
                     std::complex<double>* AWorkList_d[],
                     std::complex<double>* AinvList_d[],
                     int k,
                     int N,
                     int nw,
                     int RowStride)
{
  cuDoubleComplex one = make_cuDoubleComplex(1.0f, 0.0f);
  // A^-1 - { A^-1 * dU  * Lemma^-1 } * { V' * A^(-1) }
  // per walker: [1 x N] - [1 x k] * [k x N] = [1 x N]
  cublasZgemmBatched(handle,
                     CUBLAS_OP_N,
                     CUBLAS_OP_N,
                     1,
                     N,
                     k,
                     &one,
                     (const cuDoubleComplex**)AWorkList_d,
                     1,
                     (const cuDoubleComplex**)AinvkList_d,
                     RowStride,
                     &one,
                     (cuDoubleComplex**)AinvList_d,
                     1,
                     nw);
}

void cublas_smw_update(cublasHandle_t handle,
                       std::complex<double>* AinvkList_d[],
                       std::complex<double>* AinvList_d[],
                       std::complex<double>* AinvUList_d[],
                       std::complex<double>* AWorkList_d[],
                       std::complex<double>* lemma_inv[],
                       std::complex<double>* lemma_lu[],
                       int* PivotArray,
                       int* infoArray,
                       int k,
                       int kd,
                       int M,
                       int N,
                       int nw,
                       int RowStride)
{
#ifdef DEBUG_DELAYED
  fprintf(stderr, "*** Sherman-Morrison-Woodbury Update (k = %i, %i walkers) ***\n", k, nw);
#endif
  int pitch = RowStride;
  if (M == 1)
    pitch = 1;
  cuDoubleComplex one = make_cuDoubleComplex(1.0f, 0.0f);

  // LU decomposition needs to be updated
  cublasZgetrfBatched(handle, k, (cuDoubleComplex**)lemma_lu, kd, PivotArray, infoArray, nw);

#ifdef USE_TRSM
  if (M == 1)
  {
    // {-A^-1 * dU } * Lemma^(-1) => solve for y: Lemma y * (L * U) = (y * L) * U = -A^-1 * dU
    // z * U = -A^-1 *dU
    cublasZtrsmBatched(handle,
                       CUBLAS_SIDE_RIGHT,
                       CUBLAS_FILL_MODE_UPPER,
                       CUBLAS_OP_N,
                       CUBLAS_DIAG_NON_UNIT,
                       M,
                       k,
                       &one,
                       (const cuDoubleComplex**)lemma_lu,
                       kd,
                       (cuDoubleComplex**)AWorkList_d,
                       pitch,
                       nw);
    // y * L = z => y = {-A^-1 * dU } * Lemma^(-1)
    cublasZtrsmBatched(handle,
                       CUBLAS_SIDE_RIGHT,
                       CUBLAS_FILL_MODE_LOWER,
                       CUBLAS_OP_N,
                       CUBLAS_DIAG_UNIT,
                       M,
                       k,
                       &one,
                       (const cuDoubleComplex**)lemma_lu,
                       kd,
                       (cuDoubleComplex**)AWorkList_d,
                       pitch,
                       nw);
    // A^-1 + { -A^-1 * dU *  Lemma^-1 } * { V' * A^(-1) }
    // per walker: [1 x N] - [1 x k] * [k x N] = [1 x N]
    cublasZgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       M,
                       N,
                       k,
                       &one,
                       (const cuDoubleComplex**)AWorkList_d,
                       M,
                       (const cuDoubleComplex**)AinvkList_d,
                       RowStride,
                       &one,
                       (cuDoubleComplex**)AinvList_d,
                       pitch,
                       nw);
  }
  else
  {
    // Lemma^(-1) * V' * A^(-1) => solve for y: Lemma (L * U) * y = L * (U * y) = V' * A^(-1)
    // L * z = V' * A^(-1)
    cublasZtrsmBatched(handle,
                       CUBLAS_SIDE_LEFT,
                       CUBLAS_FILL_MODE_LOWER,
                       CUBLAS_OP_N,
                       CUBLAS_DIAG_UNIT,
                       k,
                       N,
                       &one,
                       (const cuDoubleComplex**)lemma_lu,
                       kd,
                       (cuDoubleComplex**)AWorkList_d,
                       k,
                       nw);
    // U * y = z => y = Lemma^(-1) * V' * A^(-1)
    cublasZtrsmBatched(handle,
                       CUBLAS_SIDE_LEFT,
                       CUBLAS_FILL_MODE_UPPER,
                       CUBLAS_OP_N,
                       CUBLAS_DIAG_NON_UNIT,
                       k,
                       N,
                       &one,
                       (const cuDoubleComplex**)lemma_lu,
                       kd,
                       (cuDoubleComplex**)AWorkList_d,
                       k,
                       nw);
    // A^-1 + { -A^-1 * dU } * { Lemma^-1 * V' * A^(-1) }
    // per walker: [M x N] - [M x k] * [k x N] = [M x N]
    cublasZgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       M,
                       N,
                       k,
                       &one,
                       (const cuDoubleComplex**)AinvUList_d,
                       RowStride,
                       (const cuDoubleComplex**)AWorkList_d,
                       k,
                       &one,
                       (cuDoubleComplex**)AinvList_d,
                       pitch,
                       nw);
  }
#else
  cuDoubleComplex zero = make_cuDoubleComplex(0.0f, 0.0f);
  // Calculate Lemma Inverse and store it in lemma_d
  // per walker: [k x k]^-1 = [k x k]
  cublasZgetriBatched(handle,
                      k,
                      (const cuDoubleComplex**)lemma_lu,
                      kd,
                      PivotArray,
                      (cuDoubleComplex**)lemma_inv,
                      k,
                      infoArray,
                      nw);
  // Calculate new A inverse using Sherman-Morrison-Woodbury formula
  if (M == 1) // row update can use different order to save flops
  {
    // { -A^-1 * dU } * Lemma^-1
    // per walker: [M x k] * [k x k] = [M x k]
    cublasZgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       M,
                       k,
                       k,
                       &one,
                       (const cuDoubleComplex**)AinvUList_d,
                       RowStride,
                       (const cuDoubleComplex**)lemma_inv,
                       k,
                       &zero,
                       (cuDoubleComplex**)AWorkList_d,
                       M,
                       nw);
    // A^-1 - { A^-1 * dU  * Lemma^-1 } * { V' * A^(-1) }
    // per walker: [M x N] - [M x k] * [k x N] = [M x N]
    cublasZgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       M,
                       N,
                       k,
                       &one,
                       (const cuDoubleComplex**)AWorkList_d,
                       M,
                       (const cuDoubleComplex**)AinvkList_d,
                       RowStride,
                       &one,
                       (cuDoubleComplex**)AinvList_d,
                       pitch,
                       nw);
  }
  else
  {
    // Need to use this matrix order for the overall update as AinvList and AinvkList have overlapping memory
    // Lemma^-1 * V' * A^(-1)
    // per walker: [k x k] * [k x N] = [k x N]
    cublasZgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       k,
                       N,
                       k,
                       &one,
                       (const cuDoubleComplex**)lemma_inv,
                       k,
                       (const cuDoubleComplex**)AinvkList_d,
                       RowStride,
                       &zero,
                       (cuDoubleComplex**)AWorkList_d,
                       k,
                       nw);
    // A^-1 + { -A^-1 * dU } * { Lemma^-1 * V' * A^(-1) }
    // per walker: [M x N] - [M x k] * [k x N] = [M x N]
    cublasZgemmBatched(handle,
                       CUBLAS_OP_N,
                       CUBLAS_OP_N,
                       M,
                       N,
                       k,
                       &one,
                       (const cuDoubleComplex**)AinvUList_d,
                       RowStride,
                       (const cuDoubleComplex**)AWorkList_d,
                       k,
                       &one,
                       (cuDoubleComplex**)AinvList_d,
                       pitch,
                       nw);
  }
#endif
}


template<typename T>
__global__ void update_onemove(T** buff,
                               int newrow_off,
                               int row_off,
                               int newgl_off,
                               int gl_off,
                               int ainvu_off,
                               int lemma_off,
                               int lemmainv_off,
                               int awork_off,
                               int accepted,
                               int k,
                               int kstart,
                               int kdelay,
                               int rowstride)
{
  __shared__ T *my_row, *my_gl, *my_newrow, *my_newgl, *my_ainvu, *my_lemma, *my_lemmainv, *my_awork;
  int i  = blockIdx.x * BLOCKSIZE + threadIdx.x;
  int i4 = 4 * i;
  if (threadIdx.x == 0)
  {
    T* ptr    = buff[blockIdx.y];
    my_row    = ptr + row_off;
    my_gl     = ptr + gl_off;
    my_newrow = ptr + newrow_off;
    my_newgl  = ptr + newgl_off;
    // needed for updating A^-1*dU*Lemma^-1 for drift cases
    my_awork    = ptr + awork_off;
    my_lemmainv = ptr + lemmainv_off;
    // these are only needed for rejected moves
    my_ainvu = ptr + ainvu_off;
    my_lemma = ptr + lemma_off;
  }
  __syncthreads();
  // update A^-1*dU*Lemma^-1 for drift calculations (if needed)
  if (awork_off && (i <= k))
  {
    int k1       = k + 1;
    int ik1      = i * k1;
    int kk1      = k * k1;
    int ainv_off = kstart + k1;
    T value      = 0.0;
    T rejval     = 0.0;
    // need to setup lemmainv when it's not
    if (k == 0) // in other words, k=0 => i=0, k1=1, ik1=0, kk1=0
      my_lemmainv[0] = ((T)(1.0)) / my_lemma[0];
#pragma unroll
    for (int j = 0; j <= k; j++)
    {
      // ainvu matrix is Nxk
      T ainvu_val = my_ainvu[ainv_off + j * rowstride];
      T linv_ji   = my_lemmainv[ik1 + j];
      T linv_jk   = my_lemmainv[kk1 + j];
      value += linv_ji * ainvu_val;
      rejval += linv_jk * ainvu_val;
    }
    my_awork[i] = value;
    if (blockIdx.y >= accepted)
      my_awork[i] -= rejval * my_lemmainv[ik1 + k] /
          my_lemmainv[kk1 + k]; // scale factor lemmainv_ki/lemmainv_kk (remember that lemmainv pitch is k+1, not kdelay
  }
  if (i < rowstride)
  {
    if (blockIdx.y < accepted) // accepted moves
    {
      my_row[i]     = my_newrow[i];
      my_gl[i4]     = my_newgl[i4];
      my_gl[i4 + 1] = my_newgl[i4 + 1];
      my_gl[i4 + 2] = my_newgl[i4 + 2];
      my_gl[i4 + 3] = my_newgl[i4 + 3];
    }
    else // rejected moves
    {
      my_newrow[i] = my_row[i];
      // kth column needs to be set to zero
      my_ainvu[k * rowstride + i] = 0.0;
      my_newgl[i4]                = my_gl[i4];
      my_newgl[i4 + 1]            = my_gl[i4 + 1];
      my_newgl[i4 + 2]            = my_gl[i4 + 2];
      my_newgl[i4 + 3]            = my_gl[i4 + 3];
      if (i < kdelay)
        my_lemma[i] = 0.0;
      if (i == k)
        my_lemma[k] = 1.0;
    }
  }
}

void update_onemove(float* buff[],
                    int newrow_off,
                    int row_off,
                    int newgl_off,
                    int gl_off,
                    int ainvu_off,
                    int lemma_off,
                    int lemmainv_off,
                    int awork_off,
                    int accepted,
                    int k,
                    int kstart,
                    int kdelay,
                    int rowstride,
                    int num)
{
  dim3 dimBlock(BLOCKSIZE);
  dim3 dimGrid((rowstride + BLOCKSIZE - 1) / BLOCKSIZE, num);
  update_onemove<float><<<dimGrid, dimBlock>>>(buff,
                                               newrow_off,
                                               row_off,
                                               newgl_off,
                                               gl_off,
                                               ainvu_off,
                                               lemma_off,
                                               lemmainv_off,
                                               awork_off,
                                               accepted,
                                               k,
                                               kstart,
                                               kdelay,
                                               rowstride);
}

void update_onemove(double* buff[],
                    int newrow_off,
                    int row_off,
                    int newgl_off,
                    int gl_off,
                    int ainvu_off,
                    int lemma_off,
                    int lemmainv_off,
                    int awork_off,
                    int accepted,
                    int k,
                    int kstart,
                    int kdelay,
                    int rowstride,
                    int num)
{
  dim3 dimBlock(BLOCKSIZE);
  dim3 dimGrid((rowstride + BLOCKSIZE - 1) / BLOCKSIZE, num);
  update_onemove<double><<<dimGrid, dimBlock>>>(buff,
                                                newrow_off,
                                                row_off,
                                                newgl_off,
                                                gl_off,
                                                ainvu_off,
                                                lemma_off,
                                                lemmainv_off,
                                                awork_off,
                                                accepted,
                                                k,
                                                kstart,
                                                kdelay,
                                                rowstride);
}


template<typename T>
__global__ void multi_row_copy(T** dest, T** src, int len, int offset, int rows, int stride)
{
  __shared__ T *mysrc, *mydest;
  if (threadIdx.x == 0)
  {
    mysrc  = src[blockIdx.y];
    mydest = dest[blockIdx.y];
  }
  __syncthreads();
  int i = blockIdx.x * BLOCKSIZE + threadIdx.x;
  int j = i % rows;
  int k = i / rows;
  if (i < len)
    mydest[i] = mysrc[j + k * stride + offset];
}


void multi_row_copy(float* dest[], float* src[], int len, int offset, int rows, int stride, int num)
{
  dim3 dimBlock(BLOCKSIZE);
  dim3 dimGrid((len + BLOCKSIZE - 1) / BLOCKSIZE, num);
  multi_row_copy<float><<<dimGrid, dimBlock>>>(dest, src, len, offset, rows, stride);
}


void multi_row_copy(double* dest[], double* src[], int len, int offset, int rows, int stride, int num)
{
  dim3 dimBlock(BLOCKSIZE);
  dim3 dimGrid(len / BLOCKSIZE, num);
  if (len % BLOCKSIZE)
    dimGrid.x++;
  multi_row_copy<double><<<dimGrid, dimBlock>>>(dest, src, len, offset, rows, stride);
}

template<typename T>
__global__ void calc_lemma_column(T** ainv,
                                  T** newrow,
                                  T** lemma,
                                  T** ainvu,
                                  int k,
                                  int kd,
                                  int kstart,
                                  int N,
                                  int stride)
{
  __shared__ T *myainv, *mynewrow, *mylemma_col, *myainvu_col;
  int tid = threadIdx.x;
  if (tid == 0)
  { // inputs
    myainv   = ainv[blockIdx.y];
    mynewrow = newrow[blockIdx.y] + k * stride;
    // outputs
    mylemma_col = lemma[blockIdx.y] + k * kd;
    myainvu_col = ainvu[blockIdx.y] + k * stride;
  }
  __syncthreads();
  int i  = blockIdx.x * BLOCKSIZE + tid;
  int l  = i - kstart;
  T prod = 0.0;
  if (i < N)
  {
    for (int j = 0; j < N; j++) // multiply i-th row of ainv with U
      prod += mynewrow[j] * myainv[i + j * stride];
    // now we can calculate -A^-1 * dU
    // dU = A_k-U_k => A^-1*dU = A^-1*A_k - A^-1*U_k
    myainvu_col[i] = -prod;
    // A^-1*A_k is only different from zero when we multiply same row of A^-1 with the respective column from A
    if (i == k + kstart)
      myainvu_col[i] += 1.0;
  }
  // Use the portions needed for the lemma matrix
  if ((l >= 0) && (l < kd))
    mylemma_col[l] = prod;
}

void calc_lemma_column(float* ainv[],
                       float* newrow[],
                       float* lemma[],
                       float* ainvu[],
                       int k,
                       int kd,
                       int kstart,
                       int N,
                       int stride,
                       int num)
{
  int len = stride;
  dim3 dimBlock(BLOCKSIZE);
  dim3 dimGrid((len + BLOCKSIZE - 1) / BLOCKSIZE, num);
  calc_lemma_column<float><<<dimGrid, dimBlock>>>(ainv, newrow, lemma, ainvu, k, kd, kstart, N, stride);
}

void calc_lemma_column(double* ainv[],
                       double* newrow[],
                       double* lemma[],
                       double* ainvu[],
                       int k,
                       int kd,
                       int kstart,
                       int N,
                       int stride,
                       int num)
{
  int len = stride;
  dim3 dimBlock(BLOCKSIZE);
  dim3 dimGrid((len + BLOCKSIZE - 1) / BLOCKSIZE, num);
  calc_lemma_column<double><<<dimGrid, dimBlock>>>(ainv, newrow, lemma, ainvu, k, kd, kstart, N, stride);
}


template<typename T>
__global__ void copy_delayed(T** lemma_lu, T** lemma, T** ainv_row, T** ainv_kblock, int k, int kd, int stride)
{
  __shared__ T *mylemma_lu, *mylemma, *myainv_kblock, *myainv_row;
  if (threadIdx.x == 0)
  {
    mylemma_lu    = lemma_lu[blockIdx.y];
    mylemma       = lemma[blockIdx.y];
    myainv_kblock = ainv_kblock[blockIdx.y];
    myainv_row    = ainv_row[blockIdx.y];
  }
  __syncthreads();
  int i = blockIdx.x * BLOCKSIZE + threadIdx.x;
  if (i < stride)
    myainv_row[i] = myainv_kblock[i * stride + k];
  if (i < (k + 1) * kd)
    mylemma_lu[i] = mylemma[i];
}

void copy_delayed(float* lemma_lu[],
                  float* lemma[],
                  float* ainv_row[],
                  float* ainv_kblock[],
                  int k,
                  int kd,
                  int stride,
                  int num)
{
  int len = stride;
  if ((k + 1) * kd > len)
    len = (k + 1) * kd;
  dim3 dimBlock(BLOCKSIZE);
  dim3 dimGrid((len + BLOCKSIZE - 1) / BLOCKSIZE, num);
  copy_delayed<float><<<dimGrid, dimBlock>>>(lemma_lu, lemma, ainv_row, ainv_kblock, k, kd, stride);
}

void copy_delayed(double* lemma_lu[],
                  double* lemma[],
                  double* ainv_row[],
                  double* ainv_kblock[],
                  int k,
                  int kd,
                  int stride,
                  int num)
{
  int len = stride;
  if ((k + 1) * kd > len)
    len = (k + 1) * kd;
  dim3 dimBlock(BLOCKSIZE);
  dim3 dimGrid((len + BLOCKSIZE - 1) / BLOCKSIZE, num);
  copy_delayed<double><<<dimGrid, dimBlock>>>(lemma_lu, lemma, ainv_row, ainv_kblock, k, kd, stride);
}


template<typename T>
__global__ void copy_update_block(T** lemma_lu, T** lemma, T** ainv_work, T** ainv_kblock, int k1, int kd, int stride)
{
  __shared__ T *mylemma_lu, *mylemma, *myainv_kblock, *myainv_work;
  int tid = threadIdx.x;
  if (tid == 0)
  {
    mylemma_lu    = lemma_lu[blockIdx.y];
    mylemma       = lemma[blockIdx.y];
    myainv_kblock = ainv_kblock[blockIdx.y];
    myainv_work   = ainv_work[blockIdx.y];
  }
  __syncthreads();
  int i = blockIdx.x * BLOCKSIZE + tid;
  int j = i % k1;
  int k = i / k1;
  if (i < k1 * stride)
    myainv_work[i] = myainv_kblock[j + k * stride];
  if (i < k1 * kd)
    mylemma_lu[i] = mylemma[i];
}

void copy_update_block(float* lemma_lu[],
                       float* lemma[],
                       float* ainv_work[],
                       float* ainv_kblock[],
                       int k,
                       int kd,
                       int stride,
                       int num)
{
  int len = (k + 1) * stride;
  dim3 dimBlock(BLOCKSIZE);
  dim3 dimGrid((len + BLOCKSIZE - 1) / BLOCKSIZE, num);
  copy_update_block<float><<<dimGrid, dimBlock>>>(lemma_lu, lemma, ainv_work, ainv_kblock, k + 1, kd, stride);
}

void copy_update_block(double* lemma_lu[],
                       double* lemma[],
                       double* ainv_work[],
                       double* ainv_kblock[],
                       int k,
                       int kd,
                       int stride,
                       int num)
{
  int len = (k + 1) * stride;
  dim3 dimBlock(BLOCKSIZE);
  dim3 dimGrid((len + BLOCKSIZE - 1) / BLOCKSIZE, num);
  copy_update_block<double><<<dimGrid, dimBlock>>>(lemma_lu, lemma, ainv_work, ainv_kblock, k + 1, kd, stride);
}


template<typename T>
__global__ void calc_gradlapl_and_collect(T** lemma_lu,
                                          T** Ainv_row,
                                          T** GL_col,
                                          T* ratios,
                                          int k,
                                          int kdelay,
                                          int N,
                                          int rowstride)
{
  __shared__ T *my_row, *my_col;
  int curr = 5 * blockIdx.x;
  if (threadIdx.x == 0)
  {
    my_row       = Ainv_row[blockIdx.x];
    my_col       = GL_col[blockIdx.x];
    ratios[curr] = lemma_lu[blockIdx.x][k + k * kdelay];
  }
  __syncthreads();
  T prod = 0.0;
  if (threadIdx.x < 4)
  {
    for (int j = 0; j < N; j++)
      prod += my_row[j] * my_col[j + threadIdx.x * rowstride];
    ratios[curr + threadIdx.x + 1] = prod;
  }
}

void calc_gradlapl_and_collect(float* lemma_lu[],
                               float* Ainv_row[],
                               float* GL_col[],
                               float ratios[],
                               int k,
                               int kdelay,
                               int N,
                               int rowstride,
                               int num)
{
  dim3 dimBlock(4);
  dim3 dimGrid(num);
  calc_gradlapl_and_collect<float><<<dimGrid, dimBlock>>>(lemma_lu, Ainv_row, GL_col, ratios, k, kdelay, N, rowstride);
}

void calc_gradlapl_and_collect(double* lemma_lu[],
                               double* Ainv_row[],
                               double* GL_col[],
                               double ratios[],
                               int k,
                               int kdelay,
                               int N,
                               int rowstride,
                               int num)
{
  dim3 dimBlock(4);
  dim3 dimGrid(num);
  calc_gradlapl_and_collect<double><<<dimGrid, dimBlock>>>(lemma_lu, Ainv_row, GL_col, ratios, k, kdelay, N, rowstride);
}


template<typename T>
__global__ void calc_gradient_delayed(T** Ainv_row, T** GL_col, T* ratios, int N, int rowstride)
{
  __shared__ T *my_row, *my_col;
  if (threadIdx.x == 0)
  {
    my_row = Ainv_row[blockIdx.x];
    my_col = GL_col[blockIdx.x];
  }
  __syncthreads();
  T prod = 0.0;
  if (threadIdx.x < 3)
  {
    for (int j = 0; j < N; j++)
      prod += my_row[j] * my_col[j + threadIdx.x * rowstride];
    ratios[3 * blockIdx.x + threadIdx.x] = prod;
  }
}

void calc_gradient_delayed(float* Ainv_row[], float* GL_col[], float ratios[], int N, int rowstride, int num)
{
  dim3 dimBlock(3);
  dim3 dimGrid(num);
  calc_gradient_delayed<float><<<dimGrid, dimBlock>>>(Ainv_row, GL_col, ratios, N, rowstride);
}

void calc_gradient_delayed(double* Ainv_row[], double* GL_col[], double ratios[], int N, int rowstride, int num)
{
  dim3 dimBlock(3);
  dim3 dimGrid(num);
  calc_gradient_delayed<double><<<dimGrid, dimBlock>>>(Ainv_row, GL_col, ratios, N, rowstride);
}

#ifdef QMC_COMPLEX
void update_onemove(std::complex<float>* buff[],
                    int newrow_off,
                    int row_off,
                    int newgl_off,
                    int gl_off,
                    int ainvu_off,
                    int lemma_off,
                    int lemmainv_off,
                    int awork_off,
                    int accepted,
                    int k,
                    int kstart,
                    int kdelay,
                    int rowstride,
                    int num)
{
  dim3 dimBlock(BLOCKSIZE);
  dim3 dimGrid((rowstride + BLOCKSIZE - 1) / BLOCKSIZE, num);
  update_onemove<thrust::complex<float>><<<dimGrid, dimBlock>>>((thrust::complex<float>**)buff,
                                                                newrow_off,
                                                                row_off,
                                                                newgl_off,
                                                                gl_off,
                                                                ainvu_off,
                                                                lemma_off,
                                                                lemmainv_off,
                                                                awork_off,
                                                                accepted,
                                                                k,
                                                                kstart,
                                                                kdelay,
                                                                rowstride);
}

void update_onemove(std::complex<double>* buff[],
                    int newrow_off,
                    int row_off,
                    int newgl_off,
                    int gl_off,
                    int ainvu_off,
                    int lemma_off,
                    int lemmainv_off,
                    int awork_off,
                    int accepted,
                    int k,
                    int kstart,
                    int kdelay,
                    int rowstride,
                    int num)
{
  dim3 dimBlock(BLOCKSIZE);
  dim3 dimGrid((rowstride + BLOCKSIZE - 1) / BLOCKSIZE, num);
  update_onemove<thrust::complex<double>><<<dimGrid, dimBlock>>>((thrust::complex<double>**)buff,
                                                                 newrow_off,
                                                                 row_off,
                                                                 newgl_off,
                                                                 gl_off,
                                                                 ainvu_off,
                                                                 lemma_off,
                                                                 lemmainv_off,
                                                                 awork_off,
                                                                 accepted,
                                                                 k,
                                                                 kstart,
                                                                 kdelay,
                                                                 rowstride);
}


void multi_row_copy(std::complex<float>* dest[],
                    std::complex<float>* src[],
                    int len,
                    int offset,
                    int rows,
                    int stride,
                    int num)
{
  dim3 dimBlock(BLOCKSIZE);
  dim3 dimGrid((len + BLOCKSIZE - 1) / BLOCKSIZE, num);
  multi_row_copy<thrust::complex<float>>
      <<<dimGrid, dimBlock>>>((thrust::complex<float>**)dest, (thrust::complex<float>**)src, len, offset, rows, stride);
}


void multi_row_copy(std::complex<double>* dest[],
                    std::complex<double>* src[],
                    int len,
                    int offset,
                    int rows,
                    int stride,
                    int num)
{
  dim3 dimBlock(BLOCKSIZE);
  dim3 dimGrid(len / BLOCKSIZE, num);
  if (len % BLOCKSIZE)
    dimGrid.x++;
  multi_row_copy<thrust::complex<double>><<<dimGrid, dimBlock>>>((thrust::complex<double>**)dest,
                                                                 (thrust::complex<double>**)src,
                                                                 len,
                                                                 offset,
                                                                 rows,
                                                                 stride);
}

void calc_lemma_column(std::complex<float>* ainv[],
                       std::complex<float>* newrow[],
                       std::complex<float>* lemma[],
                       std::complex<float>* ainvu[],
                       int k,
                       int kd,
                       int kstart,
                       int N,
                       int stride,
                       int num)
{
  int len = stride;
  dim3 dimBlock(BLOCKSIZE);
  dim3 dimGrid((len + BLOCKSIZE - 1) / BLOCKSIZE, num);
  calc_lemma_column<thrust::complex<float>><<<dimGrid, dimBlock>>>((thrust::complex<float>**)ainv,
                                                                   (thrust::complex<float>**)newrow,
                                                                   (thrust::complex<float>**)lemma,
                                                                   (thrust::complex<float>**)ainvu,
                                                                   k,
                                                                   kd,
                                                                   kstart,
                                                                   N,
                                                                   stride);
}

void calc_lemma_column(std::complex<double>* ainv[],
                       std::complex<double>* newrow[],
                       std::complex<double>* lemma[],
                       std::complex<double>* ainvu[],
                       int k,
                       int kd,
                       int kstart,
                       int N,
                       int stride,
                       int num)
{
  int len = stride;
  dim3 dimBlock(BLOCKSIZE);
  dim3 dimGrid((len + BLOCKSIZE - 1) / BLOCKSIZE, num);
  calc_lemma_column<thrust::complex<double>><<<dimGrid, dimBlock>>>((thrust::complex<double>**)ainv,
                                                                    (thrust::complex<double>**)newrow,
                                                                    (thrust::complex<double>**)lemma,
                                                                    (thrust::complex<double>**)ainvu,
                                                                    k,
                                                                    kd,
                                                                    kstart,
                                                                    N,
                                                                    stride);
}

void copy_delayed(std::complex<float>* lemma_lu[],
                  std::complex<float>* lemma[],
                  std::complex<float>* ainv_row[],
                  std::complex<float>* ainv_kblock[],
                  int k,
                  int kd,
                  int stride,
                  int num)
{
  int len = stride;
  if ((k + 1) * kd > len)
    len = (k + 1) * kd;
  dim3 dimBlock(BLOCKSIZE);
  dim3 dimGrid((len + BLOCKSIZE - 1) / BLOCKSIZE, num);
  copy_delayed<thrust::complex<float>><<<dimGrid, dimBlock>>>((thrust::complex<float>**)lemma_lu,
                                                              (thrust::complex<float>**)lemma,
                                                              (thrust::complex<float>**)ainv_row,
                                                              (thrust::complex<float>**)ainv_kblock,
                                                              k,
                                                              kd,
                                                              stride);
}

void copy_delayed(std::complex<double>* lemma_lu[],
                  std::complex<double>* lemma[],
                  std::complex<double>* ainv_row[],
                  std::complex<double>* ainv_kblock[],
                  int k,
                  int kd,
                  int stride,
                  int num)
{
  int len = stride;
  if ((k + 1) * kd > len)
    len = (k + 1) * kd;
  dim3 dimBlock(BLOCKSIZE);
  dim3 dimGrid((len + BLOCKSIZE - 1) / BLOCKSIZE, num);
  copy_delayed<thrust::complex<double>><<<dimGrid, dimBlock>>>((thrust::complex<double>**)lemma_lu,
                                                               (thrust::complex<double>**)lemma,
                                                               (thrust::complex<double>**)ainv_row,
                                                               (thrust::complex<double>**)ainv_kblock,
                                                               k,
                                                               kd,
                                                               stride);
}

void copy_update_block(std::complex<float>* lemma_lu[],
                       std::complex<float>* lemma[],
                       std::complex<float>* ainv_work[],
                       std::complex<float>* ainv_kblock[],
                       int k,
                       int kd,
                       int stride,
                       int num)
{
  int len = (k + 1) * stride;
  dim3 dimBlock(BLOCKSIZE);
  dim3 dimGrid((len + BLOCKSIZE - 1) / BLOCKSIZE, num);
  copy_update_block<thrust::complex<float>><<<dimGrid, dimBlock>>>((thrust::complex<float>**)lemma_lu,
                                                                   (thrust::complex<float>**)lemma,
                                                                   (thrust::complex<float>**)ainv_work,
                                                                   (thrust::complex<float>**)ainv_kblock,
                                                                   k + 1,
                                                                   kd,
                                                                   stride);
}

void copy_update_block(std::complex<double>* lemma_lu[],
                       std::complex<double>* lemma[],
                       std::complex<double>* ainv_work[],
                       std::complex<double>* ainv_kblock[],
                       int k,
                       int kd,
                       int stride,
                       int num)
{
  int len = (k + 1) * stride;
  dim3 dimBlock(BLOCKSIZE);
  dim3 dimGrid((len + BLOCKSIZE - 1) / BLOCKSIZE, num);
  copy_update_block<thrust::complex<double>><<<dimGrid, dimBlock>>>((thrust::complex<double>**)lemma_lu,
                                                                    (thrust::complex<double>**)lemma,
                                                                    (thrust::complex<double>**)ainv_work,
                                                                    (thrust::complex<double>**)ainv_kblock,
                                                                    k + 1,
                                                                    kd,
                                                                    stride);
}

void calc_gradlapl_and_collect(std::complex<float>* lemma_lu[],
                               std::complex<float>* Ainv_row[],
                               std::complex<float>* GL_col[],
                               std::complex<float> ratios[],
                               int k,
                               int kdelay,
                               int N,
                               int rowstride,
                               int num)
{
  dim3 dimBlock(4);
  dim3 dimGrid(num);
  calc_gradlapl_and_collect<thrust::complex<float>><<<dimGrid, dimBlock>>>((thrust::complex<float>**)lemma_lu,
                                                                           (thrust::complex<float>**)Ainv_row,
                                                                           (thrust::complex<float>**)GL_col,
                                                                           (thrust::complex<float>*)ratios,
                                                                           k,
                                                                           kdelay,
                                                                           N,
                                                                           rowstride);
  cudaDeviceSynchronize();
}

void calc_gradlapl_and_collect(std::complex<double>* lemma_lu[],
                               std::complex<double>* Ainv_row[],
                               std::complex<double>* GL_col[],
                               std::complex<double> ratios[],
                               int k,
                               int kdelay,
                               int N,
                               int rowstride,
                               int num)
{
  dim3 dimBlock(4);
  dim3 dimGrid(num);
  calc_gradlapl_and_collect<thrust::complex<double>><<<dimGrid, dimBlock>>>((thrust::complex<double>**)lemma_lu,
                                                                            (thrust::complex<double>**)Ainv_row,
                                                                            (thrust::complex<double>**)GL_col,
                                                                            (thrust::complex<double>*)ratios,
                                                                            k,
                                                                            kdelay,
                                                                            N,
                                                                            rowstride);
  cudaDeviceSynchronize();
}

void calc_gradient_delayed(std::complex<float>* Ainv_row[],
                           std::complex<float>* GL_col[],
                           std::complex<float> ratios[],
                           int N,
                           int rowstride,
                           int num)
{
  dim3 dimBlock(3);
  dim3 dimGrid(num);
  calc_gradient_delayed<thrust::complex<float>><<<dimGrid, dimBlock>>>((thrust::complex<float>**)Ainv_row,
                                                                       (thrust::complex<float>**)GL_col,
                                                                       (thrust::complex<float>*)ratios,
                                                                       N,
                                                                       rowstride);
}

void calc_gradient_delayed(std::complex<double>* Ainv_row[],
                           std::complex<double>* GL_col[],
                           std::complex<double> ratios[],
                           int N,
                           int rowstride,
                           int num)
{
  dim3 dimBlock(3);
  dim3 dimGrid(num);
  calc_gradient_delayed<thrust::complex<double>><<<dimGrid, dimBlock>>>((thrust::complex<double>**)Ainv_row,
                                                                        (thrust::complex<double>**)GL_col,
                                                                        (thrust::complex<double>*)ratios,
                                                                        N,
                                                                        rowstride);
}
#endif
