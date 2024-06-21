//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SYCL_MATRIX_UPDATE_HELPER_H
#define QMCPLUSPLUS_SYCL_MATRIX_UPDATE_HELPER_H

#include <sycl/sycl.hpp>

namespace qmcplusplus
{
namespace SYCL
{
template<typename T>
sycl::event copyAinvRow_saveGL_batched(sycl::queue& aq,
                                       const int rowchanged,
                                       const int n,
                                       const T* const Ainv[],
                                       const int lda,
                                       T* const temp[],
                                       T* const rcopy[],
                                       const T* const phi_vgl_in[],
                                       const int phi_vgl_stride,
                                       T* const dphi_out[],
                                       T* const d2phi_out[],
                                       int batch_count,
                                       const std::vector<sycl::event>& dependencies = {});

template<typename T, int DIM = 3>
sycl::event calcGradients_batched(sycl::queue& aq,
                                  const int n,
                                  const T* const Ainvrow[],
                                  const T* const dpsiMrow[],
                                  T* const grads_now,
                                  int batch_count,
                                  const std::vector<sycl::event>& dependencies = {});

template<typename T>
sycl::event add_delay_list_save_sigma_VGL_batched(sycl::queue& aq,
                                                  int* const delay_list[],
                                                  const int rowchanged,
                                                  const int delay_count,
                                                  T* const binv[],
                                                  const int binv_lda,
                                                  const T* const ratio_inv,
                                                  const T* const phi_vgl_in[],
                                                  const int phi_vgl_stride,
                                                  T* const phi_out[],
                                                  T* const dphi_out[],
                                                  T* const d2phi_out[],
                                                  const int norb,
                                                  const int n_accepted,
                                                  const int batch_count,
                                                  const std::vector<sycl::event>& dependencies = {});

template<typename T>
sycl::event applyW_batched(sycl::queue& aq,
                           const int* const delay_list[],
                           const int delay_count,
                           T* const tempMat[],
                           const int lda,
                           const int batch_cout,
                           const std::vector<sycl::event>& dependencies = {});

} // namespace SYCL
} // namespace qmcplusplus
#endif
