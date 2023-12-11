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
                                    const int batch_count);

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
                                    const int batch_count);

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
                                    const int batch_count);

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

cudaError_t calcGradients_cuda(cudaStream_t& hstream,
                               const int n,
                               const std::complex<float>* const Ainvrow[],
                               const std::complex<float>* const dpsiMrow[],
                               std::complex<float>* const grads_now,
                               const int batch_count);

cudaError_t calcGradients_cuda(cudaStream_t& hstream,
                               const int n,
                               const std::complex<double>* const Ainvrow[],
                               const std::complex<double>* const dpsiMrow[],
                               std::complex<double>* const grads_now,
                               const int batch_count);

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
                                                  const int batch_count);

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
                                                  const int batch_count);

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
                                                  const int batch_count);

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

cudaError_t applyW_batched(cudaStream_t& hstream,
                           const int* const delay_list[],
                           const int delay_count,
                           std::complex<float>* const tempMat[],
                           const int lda,
                           const int batch_count);

cudaError_t applyW_batched(cudaStream_t& hstream,
                           const int* const delay_list[],
                           const int delay_count,
                           std::complex<double>* const tempMat[],
                           const int lda,
                           const int batch_count);

cudaError_t print_delay_list_batched(cudaStream_t& hstream,
                                     int* const delay_list[],
                                     const int delay_count,
                                     const int batch_count);

} // namespace CUDA
} // namespace qmcplusplus
#endif
