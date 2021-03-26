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
#include "Platforms/CUDA/cudaError.h"
#include <stdexcept>
#include <type_traits>
#include <complex>
#include <cuComplex.h>
#include "Platforms/CUDA/cuBLAS.hpp"
#include <thrust/system/cuda/detail/core/util.h>

namespace qmcplusplus
{
namespace cuBLAS_LU
{
template<int COLBS>
__global__ void computeLogDet_kernel(const int n,
                                     const cuDoubleComplex* const LU_diags,
                                     const int* const pivots,
                                     cuDoubleComplex* logdets)
{
  const int iw                                   = blockIdx.x;
  const int block_num                            = blockIdx.y;
  const cuDoubleComplex* __restrict__ LU_diag_iw = LU_diags + iw * n;
  const int* __restrict__ pivots_iw              = pivots + iw * n;
  int n_index                                    = threadIdx.x + block_num * COLBS;
  __shared__ cuDoubleComplex logdet_vals[COLBS];
  logdet_vals[threadIdx.x] = {0.0, 0.0};
  if (n_index < n)
  {
    cuDoubleComplex diag = LU_diag_iw[n_index];
    diag.x = diag.x * (1 - 2* ((pivots_iw[n_index] - 1) == n_index));
    logdet_vals[threadIdx.x].x =
      log(sqrt(diag.x * diag.x + diag.y * diag.y));
    logdet_vals[threadIdx.x].y = atan2(diag.y, diag.x);
  }
  if (threadIdx.x == 0 && blockIdx.y == 0)
    logdets[iw] = {0.0, 0.0};
  // insure that when we reduce logdet_vals all the threads in the block are done.
  __syncthreads();
  if (threadIdx.x == 0)
  {
    __shared__ cuDoubleComplex block_sum_log_det;
    block_sum_log_det = {0.0, 0.0};
    for (int iv = 0; iv < COLBS; ++iv)
    {
      block_sum_log_det.x += logdet_vals[iv].x;
      block_sum_log_det.y += logdet_vals[iv].y;
    }
    // No atomicAdd in cuda 10 for cuComplex/cuDoubleComplex
    atomicAdd((double*)&(logdets[iw].x), block_sum_log_det.x);
    atomicAdd((double*)&(logdets[iw].y), block_sum_log_det.y);

    // // this will not work of n > COLBS
    // logdets[iw].x = block_sum_log_det.x;
    // logdets[iw].y = block_sum_log_det.y;
  }
}

template<int COLBS>
__global__ void computeLogDet_kernel(const int n,
                                     const double* const LU_diags,
                                     const int* const pivots,
                                     cuDoubleComplex* logdets)
{
  const int iw                          = blockIdx.x;
  const int block_num                   = blockIdx.y;
  const double* __restrict__ LU_diag_iw = LU_diags + iw * n;
  const int* __restrict__ pivots_iw     = pivots + iw * n;
  int n_index                           = threadIdx.x + block_num * COLBS;
  __shared__ cuDoubleComplex logdet_vals[COLBS];
  logdet_vals[threadIdx.x] = {0.0, 0.0};
  if (n_index < n)
  {
    logdet_vals[threadIdx.x].x = log(abs(LU_diag_iw[n_index]));
    logdet_vals[threadIdx.x].y = ((LU_diag_iw[n_index] < 0) != ((pivots_iw[n_index] - 1) == n_index)) * M_PI;
  }
  if (threadIdx.x == 0 && blockIdx.y == 0)
    logdets[iw] = {0.0, 0.0};
  // insure that when we reduce logdet_vals all the threads in the block are done.
  __syncthreads();
  if (threadIdx.x == 0)
  {
    __shared__ cuDoubleComplex block_sum_log_det;
    block_sum_log_det = {0.0, 0.0};
    for (int iv = 0; iv < COLBS; ++iv)
    {
      block_sum_log_det.x += logdet_vals[iv].x;
      block_sum_log_det.y += logdet_vals[iv].y;
    }
    // No atomicAdd in cuda 10 for cuComplex/cuDoubleComplex
    atomicAdd((double*)&(logdets[iw].x), block_sum_log_det.x);
    atomicAdd((double*)&(logdets[iw].y), block_sum_log_det.y);

    // // this will not work of n > COLBS
    // logdets[iw].x = block_sum_log_det.x;
    // logdets[iw].y = block_sum_log_det.y;
  }
}

/** Calculates logdets using LU_diags and pivots
 *  \param[in] LU_mat - the LU output from cublasXgetrfBatched
 *  \param[out] LU_diags - the LU_diags from the LU
 *  \param[in] batch_size - no a big deal here.
 */

cudaError_t computeLogDet_batched_impl(cudaStream_t& hstream,
                                       const int n,
                                       const double* LU_diags,
                                       const int* pivots,
                                       cuDoubleComplex* logdets,
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
  computeLogDet_kernel<COLBS><<<dimGrid, dimBlock, 0, hstream>>>(n, LU_diags, pivots, logdets);
  return cudaPeekAtLastError();
}

cudaError_t computeLogDet_batched_impl(cudaStream_t& hstream,
                                       const int n,
                                       const std::complex<double>* LU_diags,
                                       const int* pivots,
                                       cuDoubleComplex* logdets,
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
  computeLogDet_kernel<COLBS><<<dimGrid, dimBlock, 0, hstream>>>(n, CUDATYPECAST(LU_diags), pivots, logdets);
  return cudaPeekAtLastError();
}

template<typename T>
void computeLogDet_batched(cudaStream_t& hstream,
                           const int n,
                           const T* LU_diags,
                           const int* pivots,
                           std::complex<double>* log_dets,
                           const int batch_size)
{
  cudaErrorCheck(computeLogDet_batched_impl(hstream, n, LU_diags, pivots, CUDATYPECAST(log_dets),
                                            batch_size),
                 "failed to calculate log determinant values in computeLogDet_batched_impl");
}

template<typename T, int COLBS>
__global__ void computeLUDiag_kernel(const int n, const int lda, T** mat_lus, T* LU_diag)
{
  const int iw                = blockIdx.x;
  const int block_num         = blockIdx.y;
  const T* __restrict__ lu_iw = mat_lus[iw];
  T* __restrict__ LU_diag_iw  = LU_diag + iw * n;
  int n_index                 = threadIdx.x + block_num * COLBS;
  if (n_index < n)
    *(LU_diag_iw + n_index) = *(lu_iw + n_index * lda + n_index);
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
                                     double* Cs[],
                                     double* LU_diags,
                                     int* pivots,
                                     int* infos,
                                     std::complex<double>* log_dets,
                                     const int batch_size)
{
  //LU is returned in Ms
  cublasErrorCheck(cuBLAS::getrf_batched(h_cublas, n, Ms, lda, pivots, infos, batch_size),
                   "cuBLAS::getrf_batched failed in computeInverseAndDetLog_batched");
  cudaErrorCheck(computeLUDiag_batched_impl(hstream, n, lda, Ms, LU_diags, batch_size),
                 "failed to extract LU diag values at cuomputeLUDiag_batched_impl");
  cudaErrorCheck(computeLogDet_batched_impl(hstream, n, LU_diags, pivots, CUDATYPECAST(log_dets), batch_size),
                 "failed to calculate log determinant values in computeLogDet_batched_impl");
  cublasErrorCheck(cuBLAS::getri_batched(h_cublas, n, Ms, lda, pivots, Cs, lda, infos, batch_size),
                   "cuBLAS::getri_batched failed in computeInverseAndDetLog_batched");
}

template<typename T>
void computeGetrf_batched(cublasHandle_t& h_cublas,
                          const int n,
                          const int lda,
                          T* Ms[],
                          int* pivots,
                          int* infos,
                          const int batch_size)
{
  cublasErrorCheck(cuBLAS::getrf_batched(h_cublas, n, Ms, lda, pivots, infos, batch_size),
                   "cuBLAS::getrf_batched failed in computeInverseAndDetLog_batched");
}

void computeGetri_batched(cublasHandle_t& h_cublas,
                          const int n,
                          const int lda,
                          double* Ms[],
                          double* Cs[],
                          int* pivots,
                          int* infos,
                          const int batch_size)
{
  cublasErrorCheck(cuBLAS::getri_batched(h_cublas, n, Ms, lda, pivots, Cs, lda, infos, batch_size),
                   "cuBLAS::getri_batched failed in computeInverseAndDetLog_batched");
}

void computeLUDiag_batched(cudaStream_t& hstream,
                           const int n,
                           const int lda,
                           double** Ms,
                           double* LU_diags,
                           const int batch_size)
{
  cudaErrorCheck(computeLUDiag_batched_impl(hstream, n, lda, Ms, LU_diags, batch_size),
                 "failed to extract LU diag values at cuomputeLUDiag_batched_impl");
}


#ifndef NDEBUG
template<typename T, typename COMPLT, int COLBS>
__global__ void peekinvM_kernel(T** M, T** invM, int* pivots, int* infos, COMPLT* log_dets)
{
  const int iw        = blockIdx.x;
  const int block_num = blockIdx.y;
  const T* invM_iw    = invM[iw];
  const T* M_iw       = M[iw];
  COMPLT* log_dets_iw = log_dets + iw;
}

template<typename T, typename COMPLT>
cudaError_t peekinvM_batched_impl(cudaStream_t hstream,
                                  T** M,
                                  T** invM,
                                  int* pivots,
                                  int* infos,
                                  COMPLT* log_dets,
                                  const int batch_size)
{
  const int COLBS = 256;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_size, 1);
  peekinvM_kernel<T, COMPLT, COLBS><<<dimGrid, dimBlock, 0, hstream>>>(M, invM, pivots, infos, log_dets);
  return cudaPeekAtLastError();
}

void peekinvM_batched(cudaStream_t& hstream,
                      double** Ms,
                      double** invMs,
                      int* pivots,
                      int* infos,
                      std::complex<double>* log_dets,
                      const int batch_size)
{
  cudaErrorCheck(peekinvM_batched_impl(hstream, Ms, invMs, pivots, infos, CUDATYPECAST(log_dets), batch_size),
                 "failed to extract LU diag values at cuomputeLUDiag_batched_impl");
}

// explicit instantiation for debug peek function.
template cudaError_t peekinvM_batched_impl<double, cuDoubleComplex>(cudaStream_t hstream,
                                                                    double** M,
                                                                    double** invM,
                                                                    int* pivots,
                                                                    int* infos,
                                                                    cuDoubleComplex* log_dets,
                                                                    const int batch_size);
#endif

template void computeGetrf_batched<double>(cublasHandle_t& h_cublas,
                                           const int n,
                                           const int lda,
                                           double* Ms[],
                                           int* pivots,
                                           int* infos,
                                           const int batch_size);
template void computeGetrf_batched<std::complex<double>>(cublasHandle_t& h_cublas,
                                                         const int n,
                                                         const int lda,
                                                         std::complex<double>* Ms[],
                                                         int* pivots,
                                                         int* infos,
                                                         const int batch_size);

template void computeLogDet_batched<double>(cudaStream_t& hstream,
                                            const int n,
                                            const double* LU_diags,
                                            const int* pivots,
                                            std::complex<double>* log_dets,
                                            const int batch_size);

template void computeLogDet_batched<std::complex<double>>(cudaStream_t& hstream,
                                                          const int n,
                                                          const std::complex<double>* LU_diags,
                                                          const int* pivots,
                                                          std::complex<double>* log_dets,
                                                          const int batch_size);


} // namespace cuBLAS_LU
} // namespace qmcplusplus
