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


#ifndef QMCPLUSPLUS_COMPUTE_MATRIX_UPDATE_SYCL_H
#define QMCPLUSPLUS_COMPUTE_MATRIX_UPDATE_SYCL_H

#include <QueueAliases.hpp>
#include "matrix_update_helper.hpp"

namespace qmcplusplus
{

namespace compute
{

template<typename T>
void copyAinvRow_saveGL_batched(Queue<PlatformKind::SYCL>& queue,
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
  try
  {
    SYCL::copyAinvRow_saveGL_batched(queue.getNative(), rowchanged, n, Ainv, lda, temp, rcopy, phi_vgl_in,
                                     phi_vgl_stride, dphi_out, d2phi_out, batch_count);
  }
  catch (sycl::exception& e)
  {
    throw std::runtime_error(std::string("SYCL::copyAinvRow_saveGL_batched exception: ") + e.what());
  }
}

template<typename T>
void calcGradients_batched(Queue<PlatformKind::SYCL>& queue,
                           const int n,
                           const T* const Ainvrow[],
                           const T* const dpsiMrow[],
                           T* const grads_now,
                           const int batch_count)
{
  try
  {
    SYCL::calcGradients_batched(queue.getNative(), n, Ainvrow, dpsiMrow, grads_now, batch_count);
  }
  catch (sycl::exception& e)
  {
    throw std::runtime_error(std::string("SYCL::calcGradients_batched exception: ") + e.what());
  }
}

template<typename T>
void add_delay_list_save_sigma_VGL_batched(Queue<PlatformKind::SYCL>& queue,
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
  try
  {
    SYCL::add_delay_list_save_sigma_VGL_batched(queue.getNative(), delay_list, rowchanged, delay_count, binv, binv_lda,
                                                ratio_inv, phi_vgl_in, phi_vgl_stride, phi_out, dphi_out, d2phi_out,
                                                norb, n_accepted, batch_count);
  }
  catch (sycl::exception& e)
  {
    throw std::runtime_error(std::string("SYCL::add_delay_list_save_y_VGL_batched exception: ") + e.what());
  }
}


template<typename T>
void applyW_batched(Queue<PlatformKind::SYCL>& queue,
                    const int* const delay_list[],
                    const int delay_count,
                    T* const tempMat[],
                    const int lda,
                    const int batch_count)
{
  try
  {
    SYCL::applyW_batched(queue.getNative(), delay_list, delay_count, tempMat, lda, batch_count);
  }
  catch (sycl::exception& e)
  {
    throw std::runtime_error(std::string("SYCL::applyW_batched exception: ") + e.what());
  }
}


} // namespace compute
} // namespace qmcplusplus
#endif
