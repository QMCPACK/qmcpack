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

#include "cuBLAS_LU.hpp"
#include "Platforms/CUDA/CUDAruntime.hpp"
#include "Platforms/CUDA/cuBLAS.hpp"
#include "Platforms/CUDA/CUDATypeMapping.hpp"
#include "Platforms/CUDA/CUDAfill.hpp"
#include <cuComplex.h>

/** \file
 *
 */
namespace qmcplusplus
{
namespace cuBLAS_LU
{
/** Because the primary branch of atan2 for a complex number z = x + iy
 *  log(z).real = sqrt( x^2 + y^2 )
 *  log(z).imag = atan2(y,x) 
 *
 *  the expressions with the pivot are to avoid divergence (although the compiler might 
 *  figure it out) and  is just a factor we need to decide sign of each determinant term.
 */

/** because atan2(y, x) = pi for (0, -x)
 *  we use a short cut for the real valued matrices
 */
__device__ cuDoubleComplex complexDetLog(const double lu_diag, int n_index, const int pivot)
{
  cuDoubleComplex log_value;
  double lud  = lu_diag * (1 - 2 * (pivot != n_index + 1));
  log_value.x = log(abs(lud));
  log_value.y = (lud < 0) * M_PI;
  return log_value;
}

__device__ cuDoubleComplex complexDetLog(const cuDoubleComplex lu_diag, int n_index, const int pivot)
{
  cuDoubleComplex diag;
  double pivot_factor = 1 - 2 * (pivot != n_index + 1);
  diag.x              = lu_diag.x * pivot_factor;
  diag.y              = lu_diag.y * pivot_factor;
  cuDoubleComplex log_value;
  log_value.x = log(sqrt(diag.x * diag.x + diag.y * diag.y));
  log_value.y = atan2(diag.y, diag.x);
  return log_value;
}

/** computes the log determinant using the output of LU factorization.
 *
 *  logdets are assumed to be zeroed on launch of kernel.
 */
template<int COLBS, typename T>
__global__ void computeLogDet_kernel(const int n,
                                     const int lda,
                                     T** mat_lus,
                                     const int* const pivots,
                                     cuDoubleComplex* logdets)
{
  const int iw         = blockIdx.x;
  const int block_num  = blockIdx.y;
  const T* lu_iw       = mat_lus[iw];
  const int* pivots_iw = pivots + iw * n;
  __shared__ cuDoubleComplex logdet_vals[COLBS];
  logdet_vals[threadIdx.x].x = 0.0;
  logdet_vals[threadIdx.x].y = 0.0;
  T lu_diag;
  int n_index = threadIdx.x + block_num * COLBS;
  if (n_index < n)
  {
    lu_diag                  = *(lu_iw + n_index * lda + n_index);
    logdet_vals[threadIdx.x] = complexDetLog(lu_diag, n_index, *(pivots_iw + n_index));
  }
  // insure that when we reduce logdet_vals all the threads in the block are done.
  __syncthreads();
  if (threadIdx.x == 0)
  {
    cuDoubleComplex block_sum_log_det;
    block_sum_log_det.x = 0.0;
    block_sum_log_det.y = 0.0;
    for (int iv = 0; iv < COLBS; ++iv)
    {
      block_sum_log_det.x += logdet_vals[iv].x;
      block_sum_log_det.y += logdet_vals[iv].y;
    }
    // No atomicAdd in cuda 10 for cuComplex/cuDoubleComplex
    atomicAdd((double*)&(logdets[iw].x), block_sum_log_det.x);
    atomicAdd((double*)&(logdets[iw].y), block_sum_log_det.y);
  }
}

/** Calculates logdets using LU_diags and pivots
 *  \param[in] LU_mat - the LU output from cublasXgetrfBatched
 *  \param[out] LU_diags - the LU_diags from the LU
 *  \param[in] batch_size - no a big deal here.
 */
template<typename TMAT>
cudaError_t computeLogDet_batched_impl(cudaStream_t& hstream,
                                       const int n,
                                       const int lda,
                                       TMAT** LU_mat,
                                       const int* pivots,
                                       std::complex<double>* logdets,
                                       const int batch_size)
{
  // Perhaps this should throw an exception. I can think of no good reason it should ever happen other than
  // developer error.
  if (batch_size == 0 || n == 0)
    return cudaSuccess;

  CUDAfill_n(logdets, batch_size, {0.0, 0.0});

  const int COLBS          = 256;
  const int num_col_blocks = (n + COLBS - 1) / COLBS;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_size, num_col_blocks);
  computeLogDet_kernel<COLBS>
      <<<dimGrid, dimBlock, 0, hstream>>>(n, lda, castCUDAType(LU_mat), pivots, castCUDAType(logdets));
  return cudaPeekAtLastError();
}

template<typename T>
void computeLogDet_batched(cudaStream_t& hstream,
                           const int n,
                           const int lda,
                           T** LU_mat,
                           const int* pivots,
                           std::complex<double>* log_dets,
                           const int batch_size)
{
  cudaErrorCheck(computeLogDet_batched_impl(hstream, n, lda, LU_mat, pivots, log_dets, batch_size),
                 "failed to calculate log determinant values in computeLogDet_batched_impl");
}

template<typename T>
void computeGetrf_batched(cublasHandle_t& h_cublas,
                          cudaStream_t& hstream,
                          const int n,
                          const int lda,
                          T* Ms[],
                          int* pivots,
                          int* host_infos,
                          int* infos,
                          const int batch_size)
{
  //LU is returned in Ms
  cublasErrorCheck(cuBLAS::getrf_batched(h_cublas, n, Ms, lda, pivots, infos, batch_size),
                   "cuBLAS::getrf_batched failed in computeInverseAndDetLog_batched");
  cudaErrorCheck(cudaMemcpyAsync(host_infos, infos, sizeof(int) * batch_size, cudaMemcpyDeviceToHost, hstream),
                 "cudaMemcpyAsync failed copying cuBLAS::getrf_batched infos from device");
  cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!");

  for (int iw = 0; iw < batch_size; ++iw)
  {
    if (*(host_infos + iw) != 0)
    {
      std::ostringstream err_msg;
      err_msg << "cuBLAS::getrf_batched failed with return code " << *(host_infos + iw);
      throw std::runtime_error(err_msg.str());
    }
  }
}

template<typename T>
void computeInverseAndDetLog_batched(cublasHandle_t& h_cublas,
                                     cudaStream_t& hstream,
                                     const int n,
                                     const int lda,
                                     T* Ms[],
                                     T* Cs[],
                                     T* LU_diags,
                                     int* pivots,
                                     int* host_infos,
                                     int* infos,
                                     std::complex<double>* log_dets,
                                     const int batch_size)
{
  computeGetrf_batched(h_cublas, hstream, n, lda, Ms, pivots, host_infos, infos, batch_size);
  cudaErrorCheck(computeLogDet_batched_impl(hstream, n, lda, Ms, pivots, log_dets, batch_size),
                 "failed to calculate log determinant values in computeLogDet_batched_impl");
  cublasErrorCheck(cuBLAS::getri_batched(h_cublas, n, Ms, lda, pivots, Cs, lda, infos, batch_size),
                   "cuBLAS::getri_batched failed in computeInverseAndDetLog_batched");
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

template void computeGetrf_batched<double>(cublasHandle_t& h_cublas,
                                           cudaStream_t& hstream,
                                           const int n,
                                           const int lda,
                                           double* Ms[],
                                           int* pivots,
                                           int* host_infos,
                                           int* infos,
                                           const int batch_size);

template void computeGetrf_batched<std::complex<double>>(cublasHandle_t& h_cublas,
                                                         cudaStream_t& hstream,
                                                         const int n,
                                                         const int lda,
                                                         std::complex<double>* Ms[],
                                                         int* pivots,
                                                         int* host_infos,
                                                         int* infos,
                                                         const int batch_size);


template void computeLogDet_batched<std::complex<double>>(cudaStream_t& hstream,
                                                          const int n,
                                                          const int lda,
                                                          std::complex<double>** LU_mat,
                                                          const int* pivots,
                                                          std::complex<double>* log_dets,
                                                          const int batch_size);

template void computeLogDet_batched<double>(cudaStream_t& hstream,
                                            const int n,
                                            const int lda,
                                            double** LU_mat,
                                            const int* pivots,
                                            std::complex<double>* log_dets,
                                            const int batch_size);

template void computeInverseAndDetLog_batched<double>(cublasHandle_t& h_cublas,
                                     cudaStream_t& hstream,
                                     const int n,
                                     const int lda,
                                     double* Ms[],
                                     double* Cs[],
                                     double* LU_diags,
                                     int* pivots,
                                     int* host_infos,
                                     int* infos,
                                     std::complex<double>* log_dets,
                                     const int batch_size);

template void computeInverseAndDetLog_batched<std::complex<double>>(cublasHandle_t& h_cublas,
                                     cudaStream_t& hstream,
                                     const int n,
                                     const int lda,
                                     std::complex<double>* Ms[],
                                     std::complex<double>* Cs[],
                                     std::complex<double>* LU_diags,
                                     int* pivots,
                                     int* host_infos,
                                     int* infos,
                                     std::complex<double>* log_dets,
                                     const int batch_size);


} // namespace cuBLAS_LU
} // namespace qmcplusplus
