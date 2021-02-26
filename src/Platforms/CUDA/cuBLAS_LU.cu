//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include <cublas_v2.h>
#include "cuBLAS_LU.hpp"
#include "cudaError.h"
#include <stdexcept>
#include <cuComplex.h>
#include "cuBLAS.hpp"
#include <thrust/complex.h>
#include <thrust/system/cuda/detail/core/util.h>

namespace qmcplusplus
{
namespace QMC_CUDA
{
  
template<typename T, typename COMPT, int COLBS>
__global__ void computeLogDet_kernel(const int n, const T* const LU_diags, const int* const pivots, COMPT* logdets)
{
  const int iw                     = blockIdx.x;
  const int block_num              = blockIdx.y;
  const T* __restrict__ LU_diag_iw = LU_diags + iw * n;
  const int* __restrict__ pivot_iw = pivots + iw * n;
  int n_index                      = threadIdx.x + block_num * COLBS;
  __shared__ COMPT logdet_vals[COLBS];
  logdet_vals[threadIdx.x] = 0.0;
  if (n_index < n)
    logdet_vals[threadIdx.x] = log(((pivot_iw[n_index] == n_index + 1) ? LU_diag_iw[n_index] : -LU_diag_iw[n_index]));
  // insure that when we reduce logdet_vals all the threads in the block are done.
  __syncthreads();
  if(threadIdx.x == 0)
    {
      COMPT block_sum_log_det = 0.0;
      for(int iv = 0; iv < COLBS; ++iv)
	block_sum_log_det += logdet_vals[iv];
      atomicAdd(logdets + iw, block_sum_log_det);
    }
}

/** Calculates logdets using LU_diags and pivots
 *  \param[in] LU_mat - the LU output from cublasXgetrfBatched
 *  \param[out] LU_diags - the LU_diags from the LU
 *  \param[in] batch_size - no a big deal here.
 */
template<typename T, typename COMPT>
cudaError_t computeLogDet_batched_impl(cudaStream_t& hstream,
                                           const int n,
                                           const T* LU_diags,
					   const int* pivots,
					   COMPT* logdets,
                                           const int batch_size)
{
  // Perhaps this should throw an exception. I can think of no good reason it should ever happen other than
  // developer error.
  if (batch_size == 0 || n == 0)
    return cudaSuccess;

  const int COLBS          = 256;
  const int num_col_blocks = (n + COLBS - 1) / COLBS;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_size, num_col_blocks);
  computeLogDet_kernel<T, COMPT, COLBS><<<dimGrid, dimBlock, 0, hstream>>>(n, LU_diags, pivots, logdets);

  return cudaPeekAtLastError();
}

  
template<typename T, int COLBS>
__global__ void computeLUDiag_kernel(const int n, const int lda, const T* const invA[], T* LU_diag)
{
  const int iw                  = blockIdx.x;
  const int block_num           = blockIdx.y;
  const T* __restrict__ invA_iw = invA[iw * (blockDim.y * COLBS)];
  T* __restrict__ LU_diag_iw    = LU_diag + iw * n;
  int n_index                   = threadIdx.x + block_num * COLBS;
  if (n_index < n)
    *(LU_diag_iw + n_index) = *(invA_iw + n_index * lda + n_index);
}

/** Extracts the LU_diags from the LU in invA.
 *  \param[in] LU_mat - the LU output from cublasXgetrfBatched
 *  \param[out] LU_diags - the LU_diags from the LU
 *  \param[in] batch_size - no a big deal here.
 */
template<typename T>
cudaError_t computeLUDiag_batched_impl(cudaStream_t hstream,
                                           const int n,
                                           const int lda,
                                           T** LU_mat,
					   T* LU_diags,
                                           const int batch_size)
{
  // Perhaps this should throw an exception. I can think of no good reason it should ever happen other than
  // developer error.
  if (batch_size == 0 || n == 0)
    return cudaSuccess;

  const int COLBS          = 256;
  const int num_col_blocks = (n + COLBS - 1) / COLBS;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_size, num_col_blocks);
  computeLUDiag_kernel<T, COLBS><<<dimGrid, dimBlock, 0, hstream>>>(n, lda, LU_mat, LU_diags);

  return cudaPeekAtLastError();
}

/** Takes the transpose of PsiM using LU factorization calculates the log determinant and invPsiM
 *
 *  \param[inout] Ms -       pointers to pointers to working memory for Ms that are used to return invMs
 *  \param[in]    pivots -   pointer to n * nw ints allocated in device memory for pivots array.
 *  \param[in]    infos -    pointer to nw ints allocated in device memory factorization infos
 *  \param[out]   log_dets - pointer device memory for nw log determinant values to be returned, 
 *                           maybe this is supposed to be just RealType
 *  \param[in]    batch_size - if this changes over run a huge performance hit will be taken as memory allocation syncs device.
 */
void computeInverseAndDetLog_batched(cublasHandle_t& h_cublas,
				     cudaStream_t& hstream,
                                                  const int n,
                                                  const int lda,
						  double* Ms[],
						  double* LU_diags,
						  int* pivots,
						  int* infos,
                                                  double* log_dets,
                                                  const int batch_size)
{
  //LU is returned in Ms
  cublasErrorCheck(cuBLAS::getrf_batched(h_cublas, n, Ms, lda, pivots, infos, batch_size), "cuBLAS::getrf_batched failed in computeInverseAndDetLog_batched");
  cudaErrorCheck(computeLUDiag_batched_impl(hstream, n, lda, Ms, LU_diags, batch_size), "failed to extract LU diag values at cuomputeLUDiag_batched_impl");
  cudaErrorCheck(computeLogDet_batched_impl(hstream, n, LU_diags, pivots, log_dets, batch_size), "failed to calculate log determinant values in computeLogDet_batched_impl");
}


} // namespace cuBLAS_MFs
} // namespace qmcplusplus
