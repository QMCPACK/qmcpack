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


#include "QMCWaveFunctions/Fermion/matrix_update_helper.hpp"

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
                                          const T* const dphi_in[],
                                          const T* const d2phi_in[],
                                          T* const dphi_out[],
                                          T* const d2phi_out[])
{
  const int iw                      = blockIdx.x;
  const T* __restrict__ Ainv_iw     = Ainv[iw];
  T* __restrict__ temp_iw           = temp[iw];
  T* __restrict__ rcopy_iw          = rcopy[iw];
  const T* __restrict__ dphi_in_iw  = dphi_in[iw];
  const T* __restrict__ d2phi_in_iw = d2phi_in[iw];
  T* __restrict__ dphi_out_iw       = dphi_out[iw];
  T* __restrict__ d2phi_out_iw      = d2phi_out[iw];

  const int tid = threadIdx.x;
  if (tid == 0)
    temp_iw[rowchanged] -= T(1);

  const int num_col_blocks = (n + COLBS - 1) / COLBS;
  for (int ib = 0; ib < num_col_blocks; ib++)
  {
    const int col_id = ib * COLBS + threadIdx.x;
    if (col_id < n)
    {
      rcopy_iw[col_id] = Ainv_iw[rowchanged * lda + col_id];

      // the following copying data on the device is not part of SM-1
      // it is intended to copy dphiV and d2phiV from temporary to final without a separate kernel.
      dphi_out_iw[col_id * 3]     = dphi_in_iw[col_id * 3];
      dphi_out_iw[col_id * 3 + 1] = dphi_in_iw[col_id * 3 + 1];
      dphi_out_iw[col_id * 3 + 2] = dphi_in_iw[col_id * 3 + 2];
      d2phi_out_iw[col_id]        = d2phi_in_iw[col_id];
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
                                    const float* const dphi_in[],
                                    const float* const d2phi_in[],
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
                                                                             dphi_in, d2phi_in, dphi_out, d2phi_out);
  return cudaPeekAtLastError();
}

cudaError_t copyAinvRow_saveGL_cuda(cudaStream_t& hstream,
                                    const int rowchanged,
                                    const int n,
                                    const double* const Ainv[],
                                    const int lda,
                                    double* const temp[],
                                    double* const rcopy[],
                                    const double* const dphi_in[],
                                    const double* const d2phi_in[],
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
                                                                              dphi_in, d2phi_in, dphi_out, d2phi_out);
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

  __shared__ T sum[DIM][COLBS];
  const int tid = threadIdx.x;
  for (int idim = 0; idim < DIM; idim++)
    sum[idim][tid] = T(0);

  const int num_col_blocks = (n + COLBS - 1) / COLBS;
  for (int ib = 0; ib < num_col_blocks; ib++)
  {
    const int col_id = ib * COLBS + tid;
    for (int idim = 0; idim < DIM; idim++)
      if (col_id < n)
        sum[idim][tid] += invRow[col_id] * dpsiM_row[col_id * DIM + idim];
  }

  for (int iend = COLBS / 2; iend > 0; iend /= 2)
  {
    __syncthreads();
    for (int idim = 0; idim < DIM; idim++)
      if (tid < iend)
        sum[idim][tid] += sum[idim][tid + iend];
  }

  if (tid == 0)
    for (int idim = 0; idim < DIM; idim++)
      grads_now[iw * DIM + idim] = sum[idim][0];
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

template<typename T, int COLBS>
__global__ void add_delay_list_save_y_VGL_kernel(int* const delay_list[],
                                                const int rowchanged,
                                                const int delay_count,
                                                T* const binv[],
                                                const int binv_lda,
                                                const T* const ratio_inv,
                                                const T* const phi_in[],
                                                const T* const dphi_in[],
                                                const T* const d2phi_in[],
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
    int* __restrict__ delay_list_iw   = delay_list[iw];
    T* __restrict__ binvrow_iw        = binv[iw] + delay_count * binv_lda;
    const T* __restrict__ phi_in_iw   = phi_in[iw];
    const T* __restrict__ dphi_in_iw  = dphi_in[iw];
    const T* __restrict__ d2phi_in_iw = d2phi_in[iw];
    T* __restrict__ phi_out_iw        = phi_out[iw];
    T* __restrict__ dphi_out_iw       = dphi_out[iw];
    T* __restrict__ d2phi_out_iw      = d2phi_out[iw];

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
        dphi_out_iw[col_id * 3]     = dphi_in_iw[col_id * 3];
        dphi_out_iw[col_id * 3 + 1] = dphi_in_iw[col_id * 3 + 1];
        dphi_out_iw[col_id * 3 + 2] = dphi_in_iw[col_id * 3 + 2];
        d2phi_out_iw[col_id]        = d2phi_in_iw[col_id];
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

cudaError_t add_delay_list_save_y_VGL_batched(cudaStream_t& hstream,
                                              int* const delay_list[],
                                              const int rowchanged,
                                              const int delay_count,
                                              float* const binv[],
                                              const int binv_lda,
                                              const float* const ratio_inv,
                                              const float* const phi_in[],
                                              const float* const dphi_in[],
                                              const float* const d2phi_in[],
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
  add_delay_list_save_y_VGL_kernel<float, COLBS>
      <<<dimGrid, dimBlock, 0, hstream>>>(delay_list, rowchanged, delay_count, binv, binv_lda, ratio_inv, phi_in,
                                          dphi_in, d2phi_in, phi_out, dphi_out, d2phi_out, norb, n_accepted);
  return cudaPeekAtLastError();
}

cudaError_t add_delay_list_save_y_VGL_batched(cudaStream_t& hstream,
                                              int* const delay_list[],
                                              const int rowchanged,
                                              const int delay_count,
                                              double* const binv[],
                                              const int binv_lda,
                                              const double* const ratio_inv,
                                              const double* const phi_in[],
                                              const double* const dphi_in[],
                                              const double* const d2phi_in[],
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
  add_delay_list_save_y_VGL_kernel<double, COLBS>
      <<<dimGrid, dimBlock, 0, hstream>>>(delay_list, rowchanged, delay_count, binv, binv_lda, ratio_inv, phi_in,
                                          dphi_in, d2phi_in, phi_out, dphi_out, d2phi_out, norb, n_accepted);
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
      //printf("check applyW %d %d\n", col_id, delay_list_iw[col_id]);
      if (row_id >= 0)
        tempMat_iw[row_id * lda + col_id] -= T(1);
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

  //printf("YYdebug batch_count %d\n", batch_count);

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
