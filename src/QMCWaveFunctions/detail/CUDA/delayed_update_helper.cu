//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "delayed_update_helper.h"
#include "config.h"
#ifndef QMC_CUDA2HIP
#include <cuComplex.h>
#else
#include <hip/hip_complex.h>
#include "ROCm/cuda2hip.h"
#endif
#include "subtractOne.cuh"

/** helper kernel for delayed update algorithm
 * W matrix is applied and copy selected rows of Ainv into V
 */
template<typename T, int BS>
__global__ void applyW_stageV_kernel(const int* restrict delay_list_gpu, const int delay_count,
                                     T* restrict temp_gpu, const int numorbs, const int ndelay,
                                     T* restrict V_gpu, const T* restrict Ainv)
{
  int col = threadIdx.x + blockIdx.x * BS;

  // move rows of Ainv to V
  for(int row=0; row<delay_count; row++)
  {
    const T* Ainv_row = Ainv + numorbs * delay_list_gpu[row];
    T* V_row = V_gpu + numorbs * row;
    if( col<numorbs ) V_row[col] = Ainv_row[col];
  }

  // apply W to temp
  if( col<delay_count )
    temp_gpu[ndelay*delay_list_gpu[col] + col] = subtractOne<T>(temp_gpu[ndelay*delay_list_gpu[col] + col]);
}

void applyW_stageV_cuda(const int *delay_list_gpu, const int delay_count,
                        float* temp_gpu, const int numorbs, const int ndelay,
                        float* V_gpu, const float* Ainv,
                        cudaStream_t& hstream)
{
  const int BS = 128;
  const int NB = (numorbs+BS-1)/BS;
  dim3 dimBlock(BS);
  dim3 dimGrid(NB);
  applyW_stageV_kernel<float, BS><<<dimGrid, dimBlock, 0, hstream>>>
  (delay_list_gpu, delay_count, temp_gpu, numorbs, ndelay, V_gpu, Ainv);
}

void applyW_stageV_cuda(const int *delay_list_gpu, const int delay_count,
                        std::complex<float>* temp_gpu, const int numorbs, const int ndelay,
                        std::complex<float>* V_gpu, const std::complex<float>* Ainv,
                        cudaStream_t& hstream)
{
  const int BS = 128;
  const int NB = (numorbs+BS-1)/BS;
  dim3 dimBlock(BS);
  dim3 dimGrid(NB);
  applyW_stageV_kernel<cuComplex, BS><<<dimGrid, dimBlock, 0, hstream>>>
  (delay_list_gpu, delay_count, (cuComplex*)temp_gpu, numorbs, ndelay, (cuComplex*)V_gpu, (cuComplex*)Ainv);
}

void applyW_stageV_cuda(const int *delay_list_gpu, const int delay_count,
                        double* temp_gpu, const int numorbs, const int ndelay,
                        double* V_gpu, const double* Ainv,
                        cudaStream_t& hstream)
{
  const int BS = 128;
  const int NB = (numorbs+BS-1)/BS;
  dim3 dimBlock(BS);
  dim3 dimGrid(NB);
  applyW_stageV_kernel<double, BS><<<dimGrid, dimBlock, 0, hstream>>>
  (delay_list_gpu, delay_count, temp_gpu, numorbs, ndelay, V_gpu, Ainv);
}

void applyW_stageV_cuda(const int *delay_list_gpu, const int delay_count,
                        std::complex<double>* temp_gpu, const int numorbs, const int ndelay,
                        std::complex<double>* V_gpu, const std::complex<double>* Ainv,
                        cudaStream_t& hstream)
{
  const int BS = 128;
  const int NB = (numorbs+BS-1)/BS;
  dim3 dimBlock(BS);
  dim3 dimGrid(NB);
  applyW_stageV_kernel<cuDoubleComplex, BS><<<dimGrid, dimBlock, 0, hstream>>>
  (delay_list_gpu, delay_count, (cuDoubleComplex*)temp_gpu, numorbs, ndelay, (cuDoubleComplex*)V_gpu, (cuDoubleComplex*)Ainv);
}

/** create identity matrix on the device
 */
template<typename T>
__host__ __device__ __inline__ T makeZero()
{
  return T(0);
}

template<>
__host__ __device__ __inline__ cuDoubleComplex makeZero<cuDoubleComplex>()
{
  return make_cuDoubleComplex(0.0, 0.0);
}

template<typename T>
__host__ __device__ __inline__ T makeOne()
{
  return T(1);
}

template<>
__host__ __device__ __inline__ cuDoubleComplex makeOne<cuDoubleComplex>()
{
  return make_cuDoubleComplex(1.0, 0.0);
}

template<typename T, int BS>
__global__ void make_identity_matrix_kernel(const int nrows, T* restrict mat, const int lda)
{
  int col = threadIdx.x + blockIdx.x * BS;
  if(col<nrows)
  {
    for(int row = blockIdx.y * BS; row < min((blockIdx.y+1)*BS, nrows); row++)
      mat[row*lda+col] = makeZero<T>();
    if(blockIdx.x==blockIdx.y)
      mat[col*lda+col] = makeOne<T>();
  }
}

void make_identity_matrix_cuda(const int nrows, double* mat, const int lda, cudaStream_t& hstream)
{
  const int BS = 128;
  const int NB = (nrows+BS-1)/BS;
  dim3 dimBlock(BS);
  dim3 dimGrid(NB,NB);
  make_identity_matrix_kernel<double, BS><<<dimGrid, dimBlock, 0, hstream>>>
  (nrows, mat, lda);
}

void make_identity_matrix_cuda(const int nrows, std::complex<double>* mat, const int lda, cudaStream_t& hstream)
{
  const int BS = 128;
  const int NB = (nrows+BS-1)/BS;
  dim3 dimBlock(BS);
  dim3 dimGrid(NB,NB);
  make_identity_matrix_kernel<cuDoubleComplex, BS><<<dimGrid, dimBlock, 0, hstream>>>
  (nrows, (cuDoubleComplex*)mat, lda);
}

/** extract matrix diagonal
 */
template<typename T, int BS>
__global__ void extract_matrix_diagonal_kernel(const int nrows, const T* restrict mat, const int lda, T* restrict diag)
{
  int col = threadIdx.x + blockIdx.x * BS;
  if(col<nrows) diag[col] = mat[col*lda+col];
}

void extract_matrix_diagonal_cuda(const int nrows, const double* mat, const int lda, double* diag, cudaStream_t& hstream)
{
  const int BS = 128;
  const int NB = (nrows+BS-1)/BS;
  dim3 dimBlock(BS);
  dim3 dimGrid(NB);
  extract_matrix_diagonal_kernel<double, BS><<<dimGrid, dimBlock, 0, hstream>>>
  (nrows, mat, lda, diag);
}

void extract_matrix_diagonal_cuda(const int nrows, const std::complex<double>* mat, const int lda, std::complex<double>* diag, cudaStream_t& hstream)
{
  const int BS = 128;
  const int NB = (nrows+BS-1)/BS;
  dim3 dimBlock(BS);
  dim3 dimGrid(NB);
  extract_matrix_diagonal_kernel<cuDoubleComplex, BS><<<dimGrid, dimBlock, 0, hstream>>>
  (nrows, (cuDoubleComplex*)mat, lda, (cuDoubleComplex*)diag);
}

/** copy matrix with precision difference
 */

template<typename T_IN, typename T_OUT, int BS>
__global__ void copy_matrix_kernel(const int nrows, const int ncols, const T_IN* restrict mat_in, const int lda, T_OUT* restrict mat_out, const int ldb)
{
  int col = threadIdx.x + blockIdx.x * BS;
  if(col<ncols)
  {
    for(int row = blockIdx.y * BS; row < min((blockIdx.y+1)*BS, nrows); row++)
      mat_out[row*ldb+col] = (T_OUT) mat_in[row*lda+col];
  }
}

void copy_matrix_cuda(const int nrows, const int ncols, const double* mat_in, const int lda, float* mat_out, const int ldb, cudaStream_t& hstream)
{
  const int BS = 128;
  const int NB1 = (ncols+BS-1)/BS;
  const int NB2 = (nrows+BS-1)/BS;
  dim3 dimBlock(BS);
  dim3 dimGrid(NB1,NB2);

  copy_matrix_kernel<double, float, BS><<<dimGrid, dimBlock, 0, hstream>>>
  (nrows, ncols, mat_in, lda, mat_out, ldb);
}

void copy_matrix_cuda(const int nrows, const int ncols, const float* mat_in, const int lda, double* mat_out, const int ldb, cudaStream_t& hstream)
{
  const int BS = 128;
  const int NB1 = (ncols+BS-1)/BS;
  const int NB2 = (nrows+BS-1)/BS;
  dim3 dimBlock(BS);
  dim3 dimGrid(NB1,NB2);

  copy_matrix_kernel<float, double, BS><<<dimGrid, dimBlock, 0, hstream>>>
  (nrows, ncols, mat_in, lda, mat_out, ldb);
}

void copy_matrix_cuda(const int nrows, const int ncols, const std::complex<double>* mat_in, const int lda, std::complex<float>* mat_out, const int ldb, cudaStream_t& hstream)
{
  const int BS = 128;
  const int NB1 = (ncols*2+BS-1)/BS;
  const int NB2 = (nrows+BS-1)/BS;
  dim3 dimBlock(BS);
  dim3 dimGrid(NB1,NB2);

  copy_matrix_kernel<double, float, BS><<<dimGrid, dimBlock, 0, hstream>>>
  (nrows, ncols*2, (const double*)mat_in, lda*2, (float*)mat_out, ldb*2);
}

void copy_matrix_cuda(const int nrows, const int ncols, const std::complex<float>* mat_in, const int lda, std::complex<double>* mat_out, const int ldb, cudaStream_t& hstream)
{
  const int BS = 128;
  const int NB1 = (ncols*2+BS-1)/BS;
  const int NB2 = (nrows+BS-1)/BS;
  dim3 dimBlock(BS);
  dim3 dimGrid(NB1,NB2);

  copy_matrix_kernel<float, double, BS><<<dimGrid, dimBlock, 0, hstream>>>
  (nrows, ncols*2, (const float*)mat_in, lda*2, (double*)mat_out, ldb*2);
}
