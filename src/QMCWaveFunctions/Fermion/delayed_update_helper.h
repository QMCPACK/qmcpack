//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef CUDA_DELAYED_UPDATE_HELPER_H
#define CUDA_DELAYED_UPDATE_HELPER_H

#include <complex>

void applyW_stageV_cuda(const int *delay_list_gpu, const int delay_count,
                        float* temp_gpu, const int numorbs, const int ndelay,
                        float* V_gpu, const float* Ainv,
                        cudaStream_t& hstream);

void applyW_stageV_cuda(const int *delay_list_gpu, const int delay_count,
                        complex<float>* temp_gpu, const int numorbs, const int ndelay,
                        complex<float>* V_gpu, const complex<float>* Ainv,
                        cudaStream_t& hstream);

void applyW_stageV_cuda(const int *delay_list_gpu, const int delay_count,
                        double* temp_gpu, const int numorbs, const int ndelay,
                        double* V_gpu, const double* Ainv,
                        cudaStream_t& hstream);

void applyW_stageV_cuda(const int *delay_list_gpu, const int delay_count,
                        complex<double>* temp_gpu, const int numorbs, const int ndelay,
                        complex<double>* V_gpu, const complex<double>* Ainv,
                        cudaStream_t& hstream);

#endif
