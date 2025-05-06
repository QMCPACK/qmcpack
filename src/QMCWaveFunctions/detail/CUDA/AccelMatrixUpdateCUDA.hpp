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


#ifndef QMCPLUSPLUS_COMPUTE_MATRIX_UPDATE_CUDA_H
#define QMCPLUSPLUS_COMPUTE_MATRIX_UPDATE_CUDA_H

#include <QueueAliases.hpp>
#include "matrix_update_helper.hpp"
#include "delayed_update_helper.h"

namespace qmcplusplus
{

namespace compute
{

template<typename T>
void copyAinvRow_saveGL_batched(Queue<PlatformKind::CUDA>& queue,
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
                                const int batch_count)
{
  cudaErrorCheck(CUDA::copyAinvRow_saveGL_batched(queue.getNative(), rowchanged, n, Ainv, lda, temp, rcopy, phi_vgl_in,
                                                  phi_vgl_stride, dphi_out, d2phi_out, batch_count),
                 "CUDA::copyAinvRow_saveGL_cuda failed!");
}

template<typename T>
void calcGradients_batched(Queue<PlatformKind::CUDA>& queue,
                           const int n,
                           const T* const Ainvrow[],
                           const T* const dpsiMrow[],
                           T* const grads_now,
                           const int batch_count)
{
  cudaErrorCheck(CUDA::calcGradients_batched(queue.getNative(), n, Ainvrow, dpsiMrow, grads_now, batch_count),
                 "CUDA::calcGradients_cuda failed!");
}

template<typename T>
void add_delay_list_save_sigma_VGL_batched(Queue<PlatformKind::CUDA>& queue,
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
                                           const int batch_count)
{
  cudaErrorCheck(CUDA::add_delay_list_save_sigma_VGL_batched(queue.getNative(), delay_list, rowchanged, delay_count,
                                                             binv, binv_lda, ratio_inv, phi_vgl_in, phi_vgl_stride,
                                                             phi_out, dphi_out, d2phi_out, norb, n_accepted,
                                                             batch_count),
                 "CUDA::add_delay_list_save_y_VGL_batched failed!");
}


template<typename T>
void applyW_batched(Queue<PlatformKind::CUDA>& queue,
                    const int* const delay_list[],
                    const int delay_count,
                    T* const tempMat[],
                    const int lda,
                    const int batch_count)
{
  cudaErrorCheck(CUDA::applyW_batched(queue.getNative(), delay_list, delay_count, tempMat, lda, batch_count),
                 "CUDA::applyW_batched failed!");
}

template<typename T>
void applyW_stageV(Queue<PlatformKind::CUDA>& queue,
                   const int* delay_list_gpu,
                   const int delay_count,
                   T* temp_gpu,
                   const int numorbs,
                   const int ndelay,
                   T* V_gpu,
                   const T* Ainv)
{
  applyW_stageV_cuda(delay_list_gpu, delay_count, temp_gpu, numorbs, ndelay, V_gpu, Ainv, queue.getNative());
}

} // namespace compute
} // namespace qmcplusplus
#endif
