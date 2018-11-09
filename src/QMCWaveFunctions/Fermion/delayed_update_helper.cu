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



#include "determinant_update.h"

template<typename T, int BS>
__global__ void applyW_stageV_kernel(const int *delay_list_gpu, const int delay_count,
                                     T* temp_gpu, const int numorbs, const int ndelay,
                                     T* V_gpu, const T* Ainv)
{
  int col = threadIdx.x + blockIdx.x * BS;

  // move rows of Ainv to V
  for(int row=0; row<delay_count; row++)
  {
    const T* Ainv_row = Ainv + numorbs * delay_list_gpu[row];
    T* V_row = V_gpu + numorbs * row;
    if( col<numorbs ) V_row[col] = Ainv_row[col];
  }

  // apply W to temp
  if( col<delay_count ) temp_gpu[ndelay*delay_list_gpu[col] + col] += T(-1);
}

void applyW_stageV_cuda(const int *delay_list_gpu, const int delay_count,
                        float* temp_gpu, const int numorbs, const int ndelay,
                        float* V_gpu, const float* Ainv,
                        cudaStream_t& hstream)
{
  const int BS = 128;
  const int NB = (numorbs+BS-1)/BS;
  dim3 dimBlock(BS);
  dim3 dimGrid(NB);
  applyW_stageV_kernel<float, BS><<<dimGrid, dimBlock, 0, hstream>>>
  (delay_list_gpu, delay_count, temp_gpu, numorbs, ndelay, V_gpu, Ainv);
}

void applyW_stageV_cuda(const int *delay_list_gpu, const int delay_count,
                        double* temp_gpu, const int numorbs, const int ndelay,
                        double* V_gpu, const double* Ainv,
                        cudaStream_t& hstream)
{
  const int BS = 128;
  const int NB = (numorbs+BS-1)/BS;
  dim3 dimBlock(BS);
  dim3 dimGrid(NB);
  applyW_stageV_kernel<double, BS><<<dimGrid, dimBlock, 0, hstream>>>
  (delay_list_gpu, delay_count, temp_gpu, numorbs, ndelay, V_gpu, Ainv);
}
