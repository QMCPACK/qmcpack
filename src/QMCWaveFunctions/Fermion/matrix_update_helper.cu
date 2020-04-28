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
  const int COLBS = 32;
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
  const int COLBS = 32;
  dim3 dimBlock(COLBS);
  dim3 dimGrid(batch_count);
  copyAinvRow_saveGL_kernel<double, COLBS><<<dimGrid, dimBlock, 0, hstream>>>(rowchanged, n, Ainv, lda, temp, rcopy,
                                                                              dphi_in, d2phi_in, dphi_out, d2phi_out);
  return cudaPeekAtLastError();
}

} // namespace CUDA
} // namespace qmcplusplus
