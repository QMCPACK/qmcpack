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
__global__ void add_delay_list_compute_y_kernel(int* const delay_list[],
                                             const int rowchanged,
                                             const int delay_count,
                                             T* const binvrow[],
                                             const T* const p[],
                                             const T* const ratio,
                                             T* const y)
{
  // original CPU code
  //delay_list[delay_count] = rowchanged;
  // x
  //T y = c_ratio;
  //for (int i = 0; i < delay_count; i++)
  //  y += Binv[delay_count][i] * p[i];
  //Binv[delay_count][delay_count] = y = T(1) / y;

  const int tid = threadIdx.x;
  const int iw                    = blockIdx.x;
  int* __restrict__ delay_list_iw    = delay_list[iw];

  if (tid == 0)
    delay_list_iw[delay_count] = rowchanged;

  T* __restrict__ binvrow_iw = binvrow[iw];
  const T* __restrict__ p_iw = p[iw];

  __shared__ T sum[COLBS];
  sum[tid] = T(0);

  const int num_col_blocks = (delay_count + COLBS - 1) / COLBS;
  for (int ib = 0; ib < num_col_blocks; ib++)
  {
    const int col_id = ib * COLBS + tid;
    if (col_id < delay_count)
      sum[tid] += binvrow_iw[col_id] * p_iw[col_id];
  }

  for (int iend = COLBS / 2; iend > 0; iend /= 2)
  {
    __syncthreads();
    if (tid < iend)
      sum[tid] += sum[tid + iend];
  }

  if (tid == 0)
    binvrow_iw[delay_count] = y[iw] = T(1) / (ratio[iw] + sum[0]);
}

cudaError_t add_delay_list_compute_y_batched(cudaStream_t& hstream,
                                             int* const delay_list[],
                                             const int rowchanged,
                                             const int delay_count,
                                             float* const binvrow[],
                                             const float* const p[],
                                             const float* const ratio,
                                             float* const y,
                                             const int batch_count)
{
  if (batch_count == 0)
    return cudaSuccess;

  const int COLBS = 64;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count);
  add_delay_list_compute_y_kernel<float, COLBS><<<dimGrid, dimBlock, 0, hstream>>>(delay_list, rowchanged, delay_count, binvrow, p, ratio, y);
  return cudaPeekAtLastError();
}

cudaError_t add_delay_list_compute_y_batched(cudaStream_t& hstream,
                                             int* const delay_list[],
                                             const int rowchanged,
                                             const int delay_count,
                                             double* const binvrow[],
                                             const double* const p[],
                                             const double* const ratio,
                                             double* const y,
                                             const int batch_count)
{
  if (batch_count == 0)
    return cudaSuccess;

  const int COLBS = 64;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count);
  add_delay_list_compute_y_kernel<double, COLBS><<<dimGrid, dimBlock, 0, hstream>>>(delay_list, rowchanged, delay_count, binvrow, p, ratio, y);
  return cudaPeekAtLastError();
}

} // namespace CUDA
} // namespace qmcplusplus
