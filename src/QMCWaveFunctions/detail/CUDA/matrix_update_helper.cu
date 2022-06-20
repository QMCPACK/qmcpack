//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "matrix_update_helper.hpp"
#include "config.h"
#ifndef QMC_CUDA2HIP
#include <cuComplex.h>
#include <thrust/system/cuda/detail/core/util.h>
namespace qmcplusplus
{
namespace CUDA
{
using namespace thrust::cuda_cub::core;
}
} // namespace qmcplusplus
#else
#include <hip/hip_complex.h>
#include "ROCm/cuda2hip.h"
#include "uninitialized_array.hpp"
#endif
#include "subtractOne.cuh"
#include <thrust/complex.h>

namespace qmcplusplus
{
/** interface to cuBLAS_inhouse calls for different data types S/C/D/Z
 */
namespace CUDA
{

template<typename T, int COLBS>
__global__ void copyAinvRow_saveGL_kernel(const int rowchanged,
                                          const int n,
                                          const T* const Ainv[],
                                          const int lda,
                                          T* const temp[],
                                          T* const rcopy[],
                                          const T* const phi_vgl_in[],
                                          const size_t phi_vgl_stride,
                                          T* const dphi_out[],
                                          T* const d2phi_out[])
{
  const int iw                    = blockIdx.x;
  const T* __restrict__ Ainv_iw   = Ainv[iw];
  T* __restrict__ temp_iw         = temp[iw];
  T* __restrict__ rcopy_iw        = rcopy[iw];
  const T* __restrict__ phi_in_iw = phi_vgl_in[iw];
  T* __restrict__ dphi_out_iw     = dphi_out[iw];
  T* __restrict__ d2phi_out_iw    = d2phi_out[iw];

  const int tid = threadIdx.x;
  if (tid == 0)
    temp_iw[rowchanged] = subtractOne<T>(temp_iw[rowchanged]);

  const int num_col_blocks = (n + COLBS - 1) / COLBS;
  for (int ib = 0; ib < num_col_blocks; ib++)
  {
    const int col_id = ib * COLBS + threadIdx.x;
    if (col_id < n)
    {
      rcopy_iw[col_id] = Ainv_iw[rowchanged * lda + col_id];

      // the following copying data on the device is not part of SM-1
      // it is intended to copy dphiV and d2phiV from temporary to final without a separate kernel.
      dphi_out_iw[col_id * 3]     = phi_in_iw[col_id + phi_vgl_stride];
      dphi_out_iw[col_id * 3 + 1] = phi_in_iw[col_id + phi_vgl_stride * 2];
      dphi_out_iw[col_id * 3 + 2] = phi_in_iw[col_id + phi_vgl_stride * 3];
      d2phi_out_iw[col_id]        = phi_in_iw[col_id + phi_vgl_stride * 4];
    }
  }
}

cudaError_t copyAinvRow_saveGL_cuda(cudaStream_t& hstream,
                                    const int rowchanged,
                                    const int n,
                                    const float* const Ainv[],
                                    const int lda,
                                    float* const temp[],
                                    float* const rcopy[],
                                    const float* const phi_vgl_in[],
                                    const size_t phi_vgl_stride,
                                    float* const dphi_out[],
                                    float* const d2phi_out[],
                                    const int batch_count)
{
  if (batch_count == 0)
    return cudaSuccess;

  const int COLBS = 64;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count);
  copyAinvRow_saveGL_kernel<float, COLBS><<<dimGrid, dimBlock, 0, hstream>>>(rowchanged, n, Ainv, lda, temp, rcopy,
                                                                             phi_vgl_in, phi_vgl_stride, dphi_out,
                                                                             d2phi_out);
  return cudaPeekAtLastError();
}

cudaError_t copyAinvRow_saveGL_cuda(cudaStream_t& hstream,
                                    const int rowchanged,
                                    const int n,
                                    const double* const Ainv[],
                                    const int lda,
                                    double* const temp[],
                                    double* const rcopy[],
                                    const double* const phi_vgl_in[],
                                    const size_t phi_vgl_stride,
                                    double* const dphi_out[],
                                    double* const d2phi_out[],
                                    const int batch_count)
{
  if (batch_count == 0)
    return cudaSuccess;

  const int COLBS = 64;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count);
  copyAinvRow_saveGL_kernel<double, COLBS><<<dimGrid, dimBlock, 0, hstream>>>(rowchanged, n, Ainv, lda, temp, rcopy,
                                                                              phi_vgl_in, phi_vgl_stride, dphi_out,
                                                                              d2phi_out);
  return cudaPeekAtLastError();
}

cudaError_t copyAinvRow_saveGL_cuda(cudaStream_t& hstream,
                                    const int rowchanged,
                                    const int n,
                                    const std::complex<float>* const Ainv[],
                                    const int lda,
                                    std::complex<float>* const temp[],
                                    std::complex<float>* const rcopy[],
                                    const std::complex<float>* const phi_vgl_in[],
                                    const size_t phi_vgl_stride,
                                    std::complex<float>* const dphi_out[],
                                    std::complex<float>* const d2phi_out[],
                                    const int batch_count)
{
  if (batch_count == 0)
    return cudaSuccess;

  const int COLBS = 64;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count);
  copyAinvRow_saveGL_kernel<cuComplex, COLBS>
      <<<dimGrid, dimBlock, 0, hstream>>>(rowchanged, n, (const cuComplex**)Ainv, lda, (cuComplex**)temp,
                                          (cuComplex**)rcopy, (const cuComplex**)phi_vgl_in, phi_vgl_stride,
                                          (cuComplex**)dphi_out, (cuComplex**)d2phi_out);
  return cudaPeekAtLastError();
}

cudaError_t copyAinvRow_saveGL_cuda(cudaStream_t& hstream,
                                    const int rowchanged,
                                    const int n,
                                    const std::complex<double>* const Ainv[],
                                    const int lda,
                                    std::complex<double>* const temp[],
                                    std::complex<double>* const rcopy[],
                                    const std::complex<double>* const phi_vgl_in[],
                                    const size_t phi_vgl_stride,
                                    std::complex<double>* const dphi_out[],
                                    std::complex<double>* const d2phi_out[],
                                    const int batch_count)
{
  if (batch_count == 0)
    return cudaSuccess;

  const int COLBS = 64;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count);
  copyAinvRow_saveGL_kernel<cuDoubleComplex, COLBS>
      <<<dimGrid, dimBlock, 0, hstream>>>(rowchanged, n, (const cuDoubleComplex**)Ainv, lda, (cuDoubleComplex**)temp,
                                          (cuDoubleComplex**)rcopy, (const cuDoubleComplex**)phi_vgl_in, phi_vgl_stride,
                                          (cuDoubleComplex**)dphi_out, (cuDoubleComplex**)d2phi_out);
  return cudaPeekAtLastError();
}

template<typename T, int COLBS, int DIM = 3>
__global__ void calcGradients_kernel(const int n,
                                     const T* const Ainvrow[],
                                     const T* const dpsiMrow[],
                                     T* const grads_now)
{
  const int iw                    = blockIdx.x;
  const T* __restrict__ invRow    = Ainvrow[iw];
  const T* __restrict__ dpsiM_row = dpsiMrow[iw];

  constexpr int SUM_SIZE = DIM * COLBS;
  __shared__ uninitialized_array<T, SUM_SIZE> sum;
  const int tid = threadIdx.x;
  for (int idim = 0; idim < DIM; idim++)
    sum[idim * COLBS + tid] = T(0);

  const int num_col_blocks = (n + COLBS - 1) / COLBS;
  for (int ib = 0; ib < num_col_blocks; ib++)
  {
    const int col_id = ib * COLBS + tid;
    for (int idim = 0; idim < DIM; idim++)
      if (col_id < n)
        sum[idim * COLBS + tid] += invRow[col_id] * dpsiM_row[col_id * DIM + idim];
  }

  for (int iend = COLBS / 2; iend > 0; iend /= 2)
  {
    __syncthreads();
    for (int idim = 0; idim < DIM; idim++)
      if (tid < iend)
        sum[idim * COLBS + tid] += sum[idim * COLBS + tid + iend];
  }

  if (tid == 0)
    for (int idim = 0; idim < DIM; idim++)
      grads_now[iw * DIM + idim] = sum[idim * COLBS];
}

cudaError_t calcGradients_cuda(cudaStream_t& hstream,
                               const int n,
                               const float* const Ainvrow[],
                               const float* const dpsiMrow[],
                               float* const grads_now,
                               const int batch_count)
{
  if (batch_count == 0)
    return cudaSuccess;

  const int COLBS = 64;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count);
  calcGradients_kernel<float, COLBS><<<dimGrid, dimBlock, 0, hstream>>>(n, Ainvrow, dpsiMrow, grads_now);
  return cudaPeekAtLastError();
}

cudaError_t calcGradients_cuda(cudaStream_t& hstream,
                               const int n,
                               const double* const Ainvrow[],
                               const double* const dpsiMrow[],
                               double* const grads_now,
                               const int batch_count)
{
  if (batch_count == 0)
    return cudaSuccess;

  const int COLBS = 64;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count);
  calcGradients_kernel<double, COLBS><<<dimGrid, dimBlock, 0, hstream>>>(n, Ainvrow, dpsiMrow, grads_now);
  return cudaPeekAtLastError();
}

cudaError_t calcGradients_cuda(cudaStream_t& hstream,
                               const int n,
                               const std::complex<float>* const Ainvrow[],
                               const std::complex<float>* const dpsiMrow[],
                               std::complex<float>* const grads_now,
                               const int batch_count)
{
  if (batch_count == 0)
    return cudaSuccess;

  const int COLBS = 64;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count);
  calcGradients_kernel<thrust::complex<float>, COLBS>
      <<<dimGrid, dimBlock, 0, hstream>>>(n, (const thrust::complex<float>**)Ainvrow,
                                          (const thrust::complex<float>**)dpsiMrow, (thrust::complex<float>*)grads_now);
  return cudaPeekAtLastError();
}

cudaError_t calcGradients_cuda(cudaStream_t& hstream,
                               const int n,
                               const std::complex<double>* const Ainvrow[],
                               const std::complex<double>* const dpsiMrow[],
                               std::complex<double>* const grads_now,
                               const int batch_count)
{
  if (batch_count == 0)
    return cudaSuccess;

  const int COLBS = 64;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count);
  calcGradients_kernel<thrust::complex<double>, COLBS>
      <<<dimGrid, dimBlock, 0, hstream>>>(n, (const thrust::complex<double>**)Ainvrow,
                                          (const thrust::complex<double>**)dpsiMrow,
                                          (thrust::complex<double>*)grads_now);
  return cudaPeekAtLastError();
}

template<typename T, int COLBS>
__global__ void add_delay_list_save_sigma_VGL_kernel(int* const delay_list[],
                                                     const int rowchanged,
                                                     const int delay_count,
                                                     T* const binv[],
                                                     const int binv_lda,
                                                     const T* const ratio_inv,
                                                     const T* const phi_vgl_in[],
                                                     const size_t phi_vgl_stride,
                                                     T* const phi_out[],
                                                     T* const dphi_out[],
                                                     T* const d2phi_out[],
                                                     const int norb,
                                                     const int n_accepted)
{
  const int tid = threadIdx.x;
  const int iw  = blockIdx.x;

  if (iw < n_accepted)
  {
    // real accept, settle y and Z
    int* __restrict__ delay_list_iw = delay_list[iw];
    T* __restrict__ binvrow_iw      = binv[iw] + delay_count * binv_lda;
    const T* __restrict__ phi_in_iw = phi_vgl_in[iw];
    T* __restrict__ phi_out_iw      = phi_out[iw];
    T* __restrict__ dphi_out_iw     = dphi_out[iw];
    T* __restrict__ d2phi_out_iw    = d2phi_out[iw];

    if (tid == 0)
    {
      delay_list_iw[delay_count] = rowchanged;
      binvrow_iw[delay_count]    = ratio_inv[iw];
    }

    const int num_delay_count_col_blocks = (delay_count + COLBS - 1) / COLBS;
    for (int ib = 0; ib < num_delay_count_col_blocks; ib++)
    {
      const int col_id = ib * COLBS + tid;
      if (col_id < delay_count)
        binvrow_iw[col_id] *= ratio_inv[iw];
    }

    const int num_col_blocks = (norb + COLBS - 1) / COLBS;
    for (int ib = 0; ib < num_col_blocks; ib++)
    {
      const int col_id = ib * COLBS + tid;
      if (col_id < norb)
      {
        // copy phiV, dphiV and d2phiV from temporary to final without a separate kernel.
        phi_out_iw[col_id]          = phi_in_iw[col_id];
        dphi_out_iw[col_id * 3]     = phi_in_iw[col_id + phi_vgl_stride];
        dphi_out_iw[col_id * 3 + 1] = phi_in_iw[col_id + phi_vgl_stride * 2];
        dphi_out_iw[col_id * 3 + 2] = phi_in_iw[col_id + phi_vgl_stride * 3];
        d2phi_out_iw[col_id]        = phi_in_iw[col_id + phi_vgl_stride * 4];
      }
    }
  }
  else
  {
    // fake accept. Set Y, Z with zero and x with 1
    T* __restrict__ Urow_iw   = phi_out[iw];
    const int num_blocks_norb = (norb + COLBS - 1) / COLBS;
    for (int ib = 0; ib < num_blocks_norb; ib++)
    {
      const int col_id = ib * COLBS + tid;
      if (col_id < norb)
        Urow_iw[col_id] = T(0);
    }

    T* __restrict__ binv_iw          = binv[iw];
    const int num_blocks_delay_count = (delay_count + COLBS - 1) / COLBS;
    for (int ib = 0; ib < num_blocks_delay_count; ib++)
    {
      const int col_id = ib * COLBS + tid;
      if (col_id < delay_count)
        binv_iw[delay_count * binv_lda + col_id] = binv_iw[delay_count + binv_lda * col_id] = T(0);
    }

    int* __restrict__ delay_list_iw = delay_list[iw];
    if (tid == 0)
    {
      binv_iw[delay_count * binv_lda + delay_count] = T(1);
      delay_list_iw[delay_count]                    = -1;
    }
  }
}

cudaError_t add_delay_list_save_sigma_VGL_batched(cudaStream_t& hstream,
                                                  int* const delay_list[],
                                                  const int rowchanged,
                                                  const int delay_count,
                                                  float* const binv[],
                                                  const int binv_lda,
                                                  const float* const ratio_inv,
                                                  const float* const phi_vgl_in[],
                                                  const size_t phi_vgl_stride,
                                                  float* const phi_out[],
                                                  float* const dphi_out[],
                                                  float* const d2phi_out[],
                                                  const int norb,
                                                  const int n_accepted,
                                                  const int batch_count)
{
  if (batch_count == 0)
    return cudaSuccess;

  const int COLBS = 64;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count);
  add_delay_list_save_sigma_VGL_kernel<float, COLBS>
      <<<dimGrid, dimBlock, 0, hstream>>>(delay_list, rowchanged, delay_count, binv, binv_lda, ratio_inv, phi_vgl_in,
                                          phi_vgl_stride, phi_out, dphi_out, d2phi_out, norb, n_accepted);
  return cudaPeekAtLastError();
}

cudaError_t add_delay_list_save_sigma_VGL_batched(cudaStream_t& hstream,
                                                  int* const delay_list[],
                                                  const int rowchanged,
                                                  const int delay_count,
                                                  double* const binv[],
                                                  const int binv_lda,
                                                  const double* const ratio_inv,
                                                  const double* const phi_vgl_in[],
                                                  const size_t phi_vgl_stride,
                                                  double* const phi_out[],
                                                  double* const dphi_out[],
                                                  double* const d2phi_out[],
                                                  const int norb,
                                                  const int n_accepted,
                                                  const int batch_count)
{
  if (batch_count == 0)
    return cudaSuccess;

  const int COLBS = 64;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count);
  add_delay_list_save_sigma_VGL_kernel<double, COLBS>
      <<<dimGrid, dimBlock, 0, hstream>>>(delay_list, rowchanged, delay_count, binv, binv_lda, ratio_inv, phi_vgl_in,
                                          phi_vgl_stride, phi_out, dphi_out, d2phi_out, norb, n_accepted);
  return cudaPeekAtLastError();
}

cudaError_t add_delay_list_save_sigma_VGL_batched(cudaStream_t& hstream,
                                                  int* const delay_list[],
                                                  const int rowchanged,
                                                  const int delay_count,
                                                  std::complex<float>* const binv[],
                                                  const int binv_lda,
                                                  const std::complex<float>* const ratio_inv,
                                                  const std::complex<float>* const phi_vgl_in[],
                                                  const size_t phi_vgl_stride,
                                                  std::complex<float>* const phi_out[],
                                                  std::complex<float>* const dphi_out[],
                                                  std::complex<float>* const d2phi_out[],
                                                  const int norb,
                                                  const int n_accepted,
                                                  const int batch_count)
{
  if (batch_count == 0)
    return cudaSuccess;

  const int COLBS = 64;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count);
  add_delay_list_save_sigma_VGL_kernel<thrust::complex<float>, COLBS>
      <<<dimGrid, dimBlock, 0, hstream>>>(delay_list, rowchanged, delay_count, (thrust::complex<float>**)binv, binv_lda,
                                          (const thrust::complex<float>*)ratio_inv,
                                          (const thrust::complex<float>**)phi_vgl_in, phi_vgl_stride,
                                          (thrust::complex<float>**)phi_out, (thrust::complex<float>**)dphi_out,
                                          (thrust::complex<float>**)d2phi_out, norb, n_accepted);
  return cudaPeekAtLastError();
}

cudaError_t add_delay_list_save_sigma_VGL_batched(cudaStream_t& hstream,
                                                  int* const delay_list[],
                                                  const int rowchanged,
                                                  const int delay_count,
                                                  std::complex<double>* const binv[],
                                                  const int binv_lda,
                                                  const std::complex<double>* const ratio_inv,
                                                  const std::complex<double>* const phi_vgl_in[],
                                                  const size_t phi_vgl_stride,
                                                  std::complex<double>* const phi_out[],
                                                  std::complex<double>* const dphi_out[],
                                                  std::complex<double>* const d2phi_out[],
                                                  const int norb,
                                                  const int n_accepted,
                                                  const int batch_count)
{
  if (batch_count == 0)
    return cudaSuccess;

  const int COLBS = 64;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count);
  add_delay_list_save_sigma_VGL_kernel<thrust::complex<double>, COLBS>
      <<<dimGrid, dimBlock, 0, hstream>>>(delay_list, rowchanged, delay_count, (thrust::complex<double>**)binv,
                                          binv_lda, (const thrust::complex<double>*)ratio_inv,
                                          (const thrust::complex<double>**)phi_vgl_in, phi_vgl_stride,
                                          (thrust::complex<double>**)phi_out, (thrust::complex<double>**)dphi_out,
                                          (thrust::complex<double>**)d2phi_out, norb, n_accepted);
  return cudaPeekAtLastError();
}

template<typename T, int COLBS>
__global__ void applyW_kernel(const int* const delay_list[], const int delay_count, T* const tempMat[], const int lda)
{
  const int iw                          = blockIdx.x;
  const int* __restrict__ delay_list_iw = delay_list[iw];
  T* __restrict__ tempMat_iw            = tempMat[iw];

  const int tid        = threadIdx.x;
  const int num_blocks = (delay_count + COLBS - 1) / COLBS;
  for (int ib = 0; ib < num_blocks; ib++)
  {
    const int col_id = ib * COLBS + tid;
    if (col_id < delay_count)
    {
      const int row_id = delay_list_iw[col_id];
      if (row_id >= 0)
        tempMat_iw[row_id * lda + col_id] = subtractOne<T>(tempMat_iw[row_id * lda + col_id]);
    }
  }
}

cudaError_t applyW_batched(cudaStream_t& hstream,
                           const int* const delay_list[],
                           const int delay_count,
                           float* const tempMat[],
                           const int lda,
                           const int batch_count)
{
  if (batch_count == 0)
    return cudaSuccess;

  const int COLBS = 32;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count);
  applyW_kernel<float, COLBS><<<dimGrid, dimBlock, 0, hstream>>>(delay_list, delay_count, tempMat, lda);
  return cudaPeekAtLastError();
}

cudaError_t applyW_batched(cudaStream_t& hstream,
                           const int* const delay_list[],
                           const int delay_count,
                           double* const tempMat[],
                           const int lda,
                           const int batch_count)
{
  if (batch_count == 0)
    return cudaSuccess;

  const int COLBS = 32;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count);
  applyW_kernel<double, COLBS><<<dimGrid, dimBlock, 0, hstream>>>(delay_list, delay_count, tempMat, lda);
  return cudaPeekAtLastError();
}

cudaError_t applyW_batched(cudaStream_t& hstream,
                           const int* const delay_list[],
                           const int delay_count,
                           std::complex<float>* const tempMat[],
                           const int lda,
                           const int batch_count)
{
  if (batch_count == 0)
    return cudaSuccess;

  const int COLBS = 32;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count);
  applyW_kernel<cuComplex, COLBS>
      <<<dimGrid, dimBlock, 0, hstream>>>(delay_list, delay_count, (cuComplex**)tempMat, lda);
  return cudaPeekAtLastError();
}

cudaError_t applyW_batched(cudaStream_t& hstream,
                           const int* const delay_list[],
                           const int delay_count,
                           std::complex<double>* const tempMat[],
                           const int lda,
                           const int batch_count)
{
  if (batch_count == 0)
    return cudaSuccess;

  const int COLBS = 32;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count);
  applyW_kernel<cuDoubleComplex, COLBS>
      <<<dimGrid, dimBlock, 0, hstream>>>(delay_list, delay_count, (cuDoubleComplex**)tempMat, lda);
  return cudaPeekAtLastError();
}

__global__ void print_delay_list_kernel(int* const delay_list[], const int delay_count)
{
  const int tid                   = threadIdx.x;
  const int iw                    = blockIdx.x;
  int* __restrict__ delay_list_iw = delay_list[iw];

  if (tid == 0)
    printf("check delay_list %p %d last %d\n", delay_list_iw, delay_list_iw[delay_count],
           delay_list_iw[(delay_count - 1 < 0) ? 0 : (delay_count - 1)]);
}

cudaError_t print_delay_list_batched(cudaStream_t& hstream,
                                     int* const delay_list[],
                                     const int delay_count,
                                     const int batch_count)
{
  if (batch_count == 0)
    return cudaSuccess;

  const int COLBS = 64;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count);
  print_delay_list_kernel<<<dimGrid, dimBlock, 0, hstream>>>(delay_list, delay_count);
  return cudaPeekAtLastError();
}

} // namespace CUDA
} // namespace qmcplusplus
