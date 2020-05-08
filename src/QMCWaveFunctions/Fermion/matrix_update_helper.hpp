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


#ifndef CUDA_MATRIX_UPDATE_HELPER_H
#define CUDA_MATRIX_UPDATE_HELPER_H

#include <complex>
#include <cuda_runtime_api.h>

namespace qmcplusplus
{
/** interface to cuBLAS_inhouse calls for different data types S/C/D/Z
 */
namespace CUDA
{
/** helper function for SM-1 Fahy update
 * substract one in temp
 * copy Ainv changed row to rcopy
 * save phi G and L as accept.
 */
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
                                    const int batch_count);

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
                                    const int batch_count);
/** calculate gradients
 */
cudaError_t calcGradients_cuda(cudaStream_t& hstream,
                               const int n,
                               const float* const Ainvrow[],
                               const float* const dpsiMrow[],
                               float* const grads_now,
                               const int batch_count);

cudaError_t calcGradients_cuda(cudaStream_t& hstream,
                               const int n,
                               const double* const Ainvrow[],
                               const double* const dpsiMrow[],
                               double* const grads_now,
                               const int batch_count);

cudaError_t add_delay_list_compute_y_batched(cudaStream_t& hstream,
                                             int* const delay_list[],
                                             const int rowchanged,
                                             const int delay_count,
                                             float* const binvrow[],
                                             const float* const p[],
                                             const float* const ratio,
                                             float* const y,
                                             const int batch_count);

cudaError_t add_delay_list_compute_y_batched(cudaStream_t& hstream,
                                             int* const delay_list[],
                                             const int rowchanged,
                                             const int delay_count,
                                             double* const binvrow[],
                                             const double* const p[],
                                             const double* const ratio,
                                             double* const y,
                                             const int batch_count);

cudaError_t fake_accept_add_delay_list_update_Binv_U_batched(cudaStream_t& hstream,
                                                             int* const delay_list[],
                                                             const int delay_count,
                                                             float* const binv[],
                                                             const int lda,
                                                             float* const Urow[],
                                                             const int norb,
                                                             const int batch_count);

cudaError_t fake_accept_add_delay_list_update_Binv_U_batched(cudaStream_t& hstream,
                                                             int* const delay_list[],
                                                             const int delay_count,
                                                             double* const binv[],
                                                             const int lda,
                                                             double* const Urow[],
                                                             const int norb,
                                                             const int batch_count);

cudaError_t applyW_batched(cudaStream_t& hstream,
                           const int* const delay_list[],
                           const int delay_count,
                           float* const tempMat[],
                           const int lda,
                           const int batch_count);

cudaError_t applyW_batched(cudaStream_t& hstream,
                           const int* const delay_list[],
                           const int delay_count,
                           double* const tempMat[],
                           const int lda,
                           const int batch_count);

} // namespace CUDA
} // namespace qmcplusplus
#endif
