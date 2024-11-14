//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef CUDA_MATRIX_UPDATE_HELPER_H
#define CUDA_MATRIX_UPDATE_HELPER_H

#include "config.h"
#ifndef QMC_CUDA2HIP
#include <cuda_runtime_api.h>
#else
#include <hip/hip_runtime.h>
#include "ROCm/cuda2hip.h"
#endif

namespace qmcplusplus
{
/** interface to cuBLAS_inhouse calls for different data types S/C/D/Z
 */
namespace CUDA
{
/** helper function for SM-1 Fahy update
 * subtract one in temp
 * copy Ainv changed row to rcopy
 * save phi G and L as accept.
 */
template<typename T>
cudaError_t copyAinvRow_saveGL_batched(cudaStream_t hstream,
                                       const int rowchanged,
                                       const int n,
                                       const T* const Ainv[],
                                       const int lda,
                                       T* const temp[],
                                       T* const rcopy[],
                                       const T* const phi_vgl_in[],
                                       const size_t phi_vgl_stride,
                                       T* const dphi_out[],
                                       T* const d2phi_out[],
                                       const int batch_count);

/** calculate gradients
 */
template<typename T>
cudaError_t calcGradients_batched(cudaStream_t hstream,
                                  const int n,
                                  const T* const Ainvrow[],
                                  const T* const dpsiMrow[],
                                  T* const grads_now,
                                  const int batch_count);

template<typename T>
cudaError_t add_delay_list_save_sigma_VGL_batched(cudaStream_t hstream,
                                                  int* const delay_list[],
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
                                                  const int n_accepted,
                                                  const int batch_count);

template<typename T>
cudaError_t applyW_batched(cudaStream_t hstream,
                           const int* const delay_list[],
                           const int delay_count,
                           T* const tempMat[],
                           const int lda,
                           const int batch_count);

cudaError_t print_delay_list_batched(cudaStream_t hstream,
                                     int* const delay_list[],
                                     const int delay_count,
                                     const int batch_count);

} // namespace CUDA
} // namespace qmcplusplus
#endif
