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

#include "matrix_update_helper.hpp"
#include <complex>

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
                                       const int batch_count,
                                       const std::vector<sycl::event>& dependencies)
{
  constexpr int COLBS = 128;

  return aq
      .parallel_for(sycl::nd_range<1>{{static_cast<size_t>(batch_count * COLBS)}, {static_cast<size_t>(COLBS)}},
                    dependencies, [=](sycl::nd_item<1> item) {
                      const int iw                    = item.get_group(0); //blockIdx.x;
                      const T* __restrict__ Ainv_iw   = Ainv[iw];
                      T* __restrict__ temp_iw         = temp[iw];
                      T* __restrict__ rcopy_iw        = rcopy[iw];
                      const T* __restrict__ phi_in_iw = phi_vgl_in[iw];
                      T* __restrict__ dphi_out_iw     = dphi_out[iw];
                      T* __restrict__ d2phi_out_iw    = d2phi_out[iw];

                      const int tid = item.get_local_id(0); //threadIdx.x;
                      if (tid == 0)
                        temp_iw[rowchanged] -= T(1); //temp_iw[rowchanged] = subtractOne<T>(temp_iw[rowchanged]);

                      const int num_col_blocks = (n + COLBS - 1) / COLBS;
                      for (int ib = 0; ib < num_col_blocks; ib++)
                      {
                        const int col_id = ib * COLBS + tid; //threadIdx.x;
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
                    });
}

template sycl::event copyAinvRow_saveGL_batched(sycl::queue& aq,
                                                const int rowchanged,
                                                const int n,
                                                const float* const Ainv[],
                                                const int lda,
                                                float* const temp[],
                                                float* const rcopy[],
                                                const float* const phi_vgl_in[],
                                                const int phi_vgl_stride,
                                                float* const dphi_out[],
                                                float* const d2phi_out[],
                                                const int batch_count,
                                                const std::vector<sycl::event>& dependencies);

template sycl::event copyAinvRow_saveGL_batched(sycl::queue& aq,
                                                const int rowchanged,
                                                const int n,
                                                const double* const Ainv[],
                                                const int lda,
                                                double* const temp[],
                                                double* const rcopy[],
                                                const double* const phi_vgl_in[],
                                                const int phi_vgl_stride,
                                                double* const dphi_out[],
                                                double* const d2phi_out[],
                                                const int batch_count,
                                                const std::vector<sycl::event>& dependencies);

template sycl::event copyAinvRow_saveGL_batched(sycl::queue& aq,
                                                const int rowchanged,
                                                const int n,
                                                const std::complex<float>* const Ainv[],
                                                const int lda,
                                                std::complex<float>* const temp[],
                                                std::complex<float>* const rcopy[],
                                                const std::complex<float>* const phi_vgl_in[],
                                                const int phi_vgl_stride,
                                                std::complex<float>* const dphi_out[],
                                                std::complex<float>* const d2phi_out[],
                                                const int batch_count,
                                                const std::vector<sycl::event>& dependencies);

template sycl::event copyAinvRow_saveGL_batched(sycl::queue& aq,
                                                const int rowchanged,
                                                const int n,
                                                const std::complex<double>* const Ainv[],
                                                const int lda,
                                                std::complex<double>* const temp[],
                                                std::complex<double>* const rcopy[],
                                                const std::complex<double>* const phi_vgl_in[],
                                                const int phi_vgl_stride,
                                                std::complex<double>* const dphi_out[],
                                                std::complex<double>* const d2phi_out[],
                                                const int batch_count,
                                                const std::vector<sycl::event>& dependencies);

template<typename T, int DIM>
sycl::event calcGradients_batched(sycl::queue& aq,
                                  const int n,
                                  const T* const Ainvrow[],
                                  const T* const dpsiMrow[],
                                  T* const grads_now,
                                  const int batch_count,
                                  const std::vector<sycl::event>& dependencies)
{
  constexpr int COLBS = 128;

  return aq.submit([&](sycl::handler& cgh) {
    cgh.depends_on(dependencies);

    sycl::local_accessor<T, 1> sum((static_cast<size_t>(DIM * COLBS)), cgh);
    cgh.parallel_for(sycl::nd_range<1>{{static_cast<size_t>(batch_count * COLBS)}, {static_cast<size_t>(COLBS)}},
                     [=](sycl::nd_item<1> item) {
                       const int iw                    = item.get_group(0); //blockIdx.x;
                       const T* __restrict__ invRow    = Ainvrow[iw];
                       const T* __restrict__ dpsiM_row = dpsiMrow[iw];

                       const int tid = item.get_local_id(0); //threadIdx.x;
                       for (int idim = 0; idim < DIM; idim++)
                         sum[idim * COLBS + tid] = T{};

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
                         item.barrier(sycl::access::fence_space::local_space);
                         for (int idim = 0; idim < DIM; idim++)
                           if (tid < iend)
                             sum[idim * COLBS + tid] += sum[idim * COLBS + tid + iend];
                       }

                       if (tid == 0)
                         for (int idim = 0; idim < DIM; idim++)
                           grads_now[iw * DIM + idim] = sum[idim * COLBS];
                     });
  });
}

template sycl::event calcGradients_batched(sycl::queue& aq,
                                           const int n,
                                           const float* const Ainvrow[],
                                           const float* const dpsiMrow[],
                                           float* const grads_now,
                                           const int batch_count,
                                           const std::vector<sycl::event>& dependencies);

template sycl::event calcGradients_batched(sycl::queue& aq,
                                           const int n,
                                           const double* const Ainvrow[],
                                           const double* const dpsiMrow[],
                                           double* const grads_now,
                                           const int batch_count,
                                           const std::vector<sycl::event>& dependencies);

template sycl::event calcGradients_batched(sycl::queue& aq,
                                           const int n,
                                           const std::complex<float>* const Ainvrow[],
                                           const std::complex<float>* const dpsiMrow[],
                                           std::complex<float>* const grads_now,
                                           const int batch_count,
                                           const std::vector<sycl::event>& dependencies);

template sycl::event calcGradients_batched(sycl::queue& aq,
                                           const int n,
                                           const std::complex<double>* const Ainvrow[],
                                           const std::complex<double>* const dpsiMrow[],
                                           std::complex<double>* const grads_now,
                                           const int batch_count,
                                           const std::vector<sycl::event>& dependencies);


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
                                                  const std::vector<sycl::event>& dependencies)
{
  constexpr int COLBS = 64;

  return aq.parallel_for(sycl::nd_range<1>{{static_cast<size_t>(batch_count * COLBS)}, {static_cast<size_t>(COLBS)}},
                         dependencies, [=](sycl::nd_item<1> item) {
                           const int tid = item.get_local_id(0); //threadIdx.x;
                           const int iw  = item.get_group(0);    //blockIdx.x;

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
                                 Urow_iw[col_id] = T{};
                             }

                             T* __restrict__ binv_iw          = binv[iw];
                             const int num_blocks_delay_count = (delay_count + COLBS - 1) / COLBS;
                             for (int ib = 0; ib < num_blocks_delay_count; ib++)
                             {
                               const int col_id = ib * COLBS + tid;
                               if (col_id < delay_count)
                                 binv_iw[delay_count * binv_lda + col_id] = binv_iw[delay_count + binv_lda * col_id] =
                                     T(0);
                             }

                             int* __restrict__ delay_list_iw = delay_list[iw];
                             if (tid == 0)
                             {
                               binv_iw[delay_count * binv_lda + delay_count] = T(1);
                               delay_list_iw[delay_count]                    = -1;
                             }
                           }
                         });
}

template sycl::event add_delay_list_save_sigma_VGL_batched(sycl::queue& aq,
                                                           int* const delay_list[],
                                                           const int rowchanged,
                                                           const int delay_count,
                                                           float* const binv[],
                                                           const int binv_lda,
                                                           const float* const ratio_inv,
                                                           const float* const phi_vgl_in[],
                                                           const int phi_vgl_stride,
                                                           float* const phi_out[],
                                                           float* const dphi_out[],
                                                           float* const d2phi_out[],
                                                           const int norb,
                                                           const int n_accepted,
                                                           const int batch_count,
                                                           const std::vector<sycl::event>& dependencies);

template sycl::event add_delay_list_save_sigma_VGL_batched(sycl::queue& aq,
                                                           int* const delay_list[],
                                                           const int rowchanged,
                                                           const int delay_count,
                                                           double* const binv[],
                                                           const int binv_lda,
                                                           const double* const ratio_inv,
                                                           const double* const phi_vgl_in[],
                                                           const int phi_vgl_stride,
                                                           double* const phi_out[],
                                                           double* const dphi_out[],
                                                           double* const d2phi_out[],
                                                           const int norb,
                                                           const int n_accepted,
                                                           const int batch_count,
                                                           const std::vector<sycl::event>& dependencies);

template sycl::event add_delay_list_save_sigma_VGL_batched(sycl::queue& aq,
                                                           int* const delay_list[],
                                                           const int rowchanged,
                                                           const int delay_count,
                                                           std::complex<float>* const binv[],
                                                           const int binv_lda,
                                                           const std::complex<float>* const ratio_inv,
                                                           const std::complex<float>* const phi_vgl_in[],
                                                           const int phi_vgl_stride,
                                                           std::complex<float>* const phi_out[],
                                                           std::complex<float>* const dphi_out[],
                                                           std::complex<float>* const d2phi_out[],
                                                           const int norb,
                                                           const int n_accepted,
                                                           const int batch_count,
                                                           const std::vector<sycl::event>& dependencies);

template sycl::event add_delay_list_save_sigma_VGL_batched(sycl::queue& aq,
                                                           int* const delay_list[],
                                                           const int rowchanged,
                                                           const int delay_count,
                                                           std::complex<double>* const binv[],
                                                           const int binv_lda,
                                                           const std::complex<double>* const ratio_inv,
                                                           const std::complex<double>* const phi_vgl_in[],
                                                           const int phi_vgl_stride,
                                                           std::complex<double>* const phi_out[],
                                                           std::complex<double>* const dphi_out[],
                                                           std::complex<double>* const d2phi_out[],
                                                           const int norb,
                                                           const int n_accepted,
                                                           const int batch_count,
                                                           const std::vector<sycl::event>& dependencies);

template<typename T>
sycl::event applyW_batched(sycl::queue& aq,
                           const int* const delay_list[],
                           const int delay_count,
                           T* const tempMat[],
                           const int lda,
                           const int batch_count,
                           const std::vector<sycl::event>& dependencies)
{
  constexpr int COLBS = 32;

  return aq.parallel_for(sycl::nd_range<1>{{static_cast<size_t>(batch_count * COLBS)}, {static_cast<size_t>(COLBS)}},
                         dependencies, [=](sycl::nd_item<1> item) {
                           const int iw                          = item.get_group(0); //blockIdx.x;
                           const int* __restrict__ delay_list_iw = delay_list[iw];
                           T* __restrict__ tempMat_iw            = tempMat[iw];

                           const T mone(-1);
                           const int tid        = item.get_local_id(0); //threadIdx.x;
                           const int num_blocks = (delay_count + COLBS - 1) / COLBS;
                           for (int ib = 0; ib < num_blocks; ib++)
                           {
                             const int col_id = ib * COLBS + tid;
                             if (col_id < delay_count)
                             {
                               const int row_id = delay_list_iw[col_id];
                               if (row_id >= 0)
                                 tempMat_iw[row_id * lda + col_id] +=
                                     mone; //subtractOne<T>(tempMat_iw[row_id * lda + col_id]);
                             }
                           }
                         });
}

template sycl::event applyW_batched(sycl::queue& aq,
                                    const int* const delay_list[],
                                    const int delay_count,
                                    float* const tempMat[],
                                    const int lda,
                                    const int batch_count,
                                    const std::vector<sycl::event>& dependencies);

template sycl::event applyW_batched(sycl::queue& aq,
                                    const int* const delay_list[],
                                    const int delay_count,
                                    double* const tempMat[],
                                    const int lda,
                                    const int batch_count,
                                    const std::vector<sycl::event>& dependencies);

template sycl::event applyW_batched(sycl::queue& aq,
                                    const int* const delay_list[],
                                    const int delay_count,
                                    std::complex<float>* const tempMat[],
                                    const int lda,
                                    const int batch_count,
                                    const std::vector<sycl::event>& dependencies);

template sycl::event applyW_batched(sycl::queue& aq,
                                    const int* const delay_list[],
                                    const int delay_count,
                                    std::complex<double>* const tempMat[],
                                    const int lda,
                                    const int batch_count,
                                    const std::vector<sycl::event>& dependencies);

} // namespace SYCL
} // namespace qmcplusplus
