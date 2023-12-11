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
#include <algorithm>
#include "Platforms/CUDA/CUDAruntime.hpp"
#include "Platforms/CUDA/cuBLAS.hpp"
#include "Platforms/CUDA/CUDATypeMapping.hpp"

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
__device__ cuDoubleComplex complexDetLog(const double lu_diag, const int n_index, const int pivot)
{
  cuDoubleComplex log_value;
  double lud  = lu_diag * (1 - 2 * (pivot != n_index + 1));
  log_value.x = log(abs(lud));
  log_value.y = (lud < 0) * M_PI;
  return log_value;
}

__device__ cuDoubleComplex complexDetLog(const cuDoubleComplex lu_diag, const int n_index, const int pivot)
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
  static_assert(COLBS != 0 && (COLBS & (COLBS - 1)) == 0, "COLBS is not a power of 2");
  const int iw         = blockIdx.x;
  const T* lu_iw       = mat_lus[iw];
  const int* pivots_iw = pivots + iw * n;

  __shared__ cuDoubleComplex logdet_vals[COLBS];
  const int tid      = threadIdx.x;
  logdet_vals[tid].x = 0.0;
  logdet_vals[tid].y = 0.0;

  const int num_col_blocks = (n + COLBS - 1) / COLBS;
  for (int block_num = 0; block_num < num_col_blocks; block_num++)
  {
    const int n_index = tid + block_num * COLBS;
    if (n_index < n)
    {
      const T lu_diag                   = *(lu_iw + n_index * lda + n_index);
      const cuDoubleComplex lu_diag_log = complexDetLog(lu_diag, n_index, *(pivots_iw + n_index));
      logdet_vals[tid].x += lu_diag_log.x;
      logdet_vals[tid].y += lu_diag_log.y;
    }
  }

  for (int iend = COLBS / 2; iend > 0; iend /= 2)
  {
    __syncthreads();
    if (tid < iend)
    {
      logdet_vals[tid].x += logdet_vals[tid + iend].x;
      logdet_vals[tid].y += logdet_vals[tid + iend].y;
    }
  }
  if (tid == 0)
    logdets[iw] = logdet_vals[0];
}

/** Calculates logdets using LU_diags and pivots
 *  \param[in] LU_mat - the LU output from cublasXgetrfBatched
 *  \param[out] LU_diags - the LU_diags from the LU
 *  \param[in] batch_size - no a big deal here.
 */
template<typename TMAT, typename T>
cudaError_t computeLogDet_batched_impl(cudaStream_t& hstream,
                                       const int n,
                                       const int lda,
                                       TMAT** LU_mat,
                                       const int* pivots,
                                       T* logdets,
                                       const int batch_size)
{
  // Perhaps this should throw an exception. I can think of no good reason it should ever happen other than
  // developer error.
  if (batch_size == 0 || n == 0)
    return cudaSuccess;

  const int COLBS = 128;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_size);
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

  if (std::any_of(host_infos, host_infos + batch_size, [](int i) { return i != 0; }))
  {
    std::ostringstream err_msg;
    err_msg << "cuBLAS::getrf_batched failed! Non-zero infos:" << std::endl;
    for (int iw = 0; iw < batch_size; ++iw)
      if (*(host_infos + iw) != 0)
        err_msg << "infos[" << iw << "] = " << *(host_infos + iw) << std::endl;
    throw std::runtime_error(err_msg.str());
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
  computeGetri_batched(h_cublas, hstream, n, lda, Ms, Cs, pivots, host_infos, infos, batch_size);
}


template<typename T>
void computeGetri_batched(cublasHandle_t& h_cublas,
                          cudaStream_t& hstream,
                          const int n,
                          const int lda,
                          T* Ms[],
                          T* Cs[],
                          int* pivots,
                          int* host_infos,
                          int* infos,
                          const int batch_size)
{
  cublasErrorCheck(cuBLAS::getri_batched(h_cublas, n, Ms, lda, pivots, Cs, lda, infos, batch_size),
                   "cuBLAS::getri_batched failed in computeInverseAndDetLog_batched");
  cudaErrorCheck(cudaMemcpyAsync(host_infos, infos, sizeof(int) * batch_size, cudaMemcpyDeviceToHost, hstream),
                 "cudaMemcpyAsync failed copying cuBLAS::getri_batched infos from device");
  cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!");

  if (std::any_of(host_infos, host_infos + batch_size, [](int i) { return i != 0; }))
  {
    std::ostringstream err_msg;
    err_msg << "cuBLAS::getri_batched failed! Non-zero infos:" << std::endl;
    for (int iw = 0; iw < batch_size; ++iw)
      if (*(host_infos + iw) != 0)
        err_msg << "infos[" << iw << "] = " << *(host_infos + iw) << std::endl;
    throw std::runtime_error(err_msg.str());
  }
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

template void computeGetri_batched<double>(cublasHandle_t& h_cublas,
                                           cudaStream_t& hstream,
                                           const int n,
                                           const int lda,
                                           double* Ms[],
                                           double* Cs[],
                                           int* pivots,
                                           int* host_infos,
                                           int* infos,
                                           const int batch_size);

template void computeGetri_batched<std::complex<double>>(cublasHandle_t& h_cublas,
                                                         cudaStream_t& hstream,
                                                         const int n,
                                                         const int lda,
                                                         std::complex<double>* Ms[],
                                                         std::complex<double>* Cs[],
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

template void computeLogDet_batched<float>(cudaStream_t& hstream,
                                           const int n,
                                           const int lda,
                                           float** LU_mat,
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
