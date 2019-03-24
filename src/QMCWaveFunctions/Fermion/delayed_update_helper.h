//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef CUDA_DELAYED_UPDATE_HELPER_H
#define CUDA_DELAYED_UPDATE_HELPER_H

#include <complex>

/** helper function for delayed update algorithm
 * W matrix is applied and copy selected rows of Ainv into V
 */

void applyW_stageV_cuda(const int *delay_list_gpu, const int delay_count,
                        float* temp_gpu, const int numorbs, const int ndelay,
                        float* V_gpu, const float* Ainv,
                        cudaStream_t& hstream);

void applyW_stageV_cuda(const int *delay_list_gpu, const int delay_count,
                        std::complex<float>* temp_gpu, const int numorbs, const int ndelay,
                        std::complex<float>* V_gpu, const std::complex<float>* Ainv,
                        cudaStream_t& hstream);

void applyW_stageV_cuda(const int *delay_list_gpu, const int delay_count,
                        double* temp_gpu, const int numorbs, const int ndelay,
                        double* V_gpu, const double* Ainv,
                        cudaStream_t& hstream);

void applyW_stageV_cuda(const int *delay_list_gpu, const int delay_count,
                        std::complex<double>* temp_gpu, const int numorbs, const int ndelay,
                        std::complex<double>* V_gpu, const std::complex<double>* Ainv,
                        cudaStream_t& hstream);

/** create identity matrix on the device
 */
void make_identity_matrix_cuda(const int nrows, double* mat, const int lda, cudaStream_t& hstream);

void make_identity_matrix_cuda(const int nrows, std::complex<double>* mat, const int lda, cudaStream_t& hstream);

#endif
