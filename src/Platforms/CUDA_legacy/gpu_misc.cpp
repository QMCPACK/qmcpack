//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ying Wai Li, yingwaili@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "config.h"
#include "gpu_misc.h"

namespace gpu
{
cudaStream_t kernelStream;
cudaStream_t memoryStream;

cudaEvent_t syncEvent;

cudaEvent_t gradientSyncDiracEvent;
cudaEvent_t gradientSyncOneBodyEvent;
cudaEvent_t gradientSyncTwoBodyEvent;
cudaEvent_t ratioSyncDiracEvent;
cudaEvent_t ratioSyncOneBodyEvent;
cudaEvent_t ratioSyncTwoBodyEvent;
cublasHandle_t cublasHandle;

size_t MaxGPUSpineSizeMB;
int rank;
int relative_rank;                     // relative rank number on the node the rank is on, counting starts at zero
int device_group_size;                 // size of the lists below
bool cudamps;                          // is set to true if Cuda MPS service is running
std::vector<int> device_group_numbers; // on node list of GPU device numbers with respect to relative rank number
std::vector<int>
    device_rank_numbers; // on node list of MPI rank numbers (absolute) with respect to relative rank number

void initCUDAStreams()
{
  cudaCheck(cudaStreamCreate(&kernelStream));
  cudaCheck(cudaStreamCreate(&memoryStream));
}

void initCUDAEvents()
{
  cudaCheck(cudaEventCreateWithFlags(&syncEvent, cudaEventDisableTiming));
  cudaCheck(cudaEventCreateWithFlags(&gradientSyncDiracEvent, cudaEventDisableTiming));
  cudaCheck(cudaEventCreateWithFlags(&gradientSyncOneBodyEvent, cudaEventDisableTiming));
  cudaCheck(cudaEventCreateWithFlags(&gradientSyncTwoBodyEvent, cudaEventDisableTiming));
  cudaCheck(cudaEventCreateWithFlags(&ratioSyncDiracEvent, cudaEventDisableTiming));
  cudaCheck(cudaEventCreateWithFlags(&ratioSyncOneBodyEvent, cudaEventDisableTiming));
  cudaCheck(cudaEventCreateWithFlags(&ratioSyncTwoBodyEvent, cudaEventDisableTiming));
}

void initCublas() { cublasCreate(&cublasHandle); }

void finalizeCUDAStreams()
{
  cudaCheck(cudaStreamDestroy(kernelStream));
  cudaCheck(cudaStreamDestroy(memoryStream));
}

void finalizeCUDAEvents()
{
  cudaCheck(cudaEventDestroy(syncEvent));
  cudaCheck(cudaEventDestroy(gradientSyncDiracEvent));
  cudaCheck(cudaEventDestroy(gradientSyncOneBodyEvent));
  cudaCheck(cudaEventDestroy(gradientSyncTwoBodyEvent));
  cudaCheck(cudaEventDestroy(ratioSyncDiracEvent));
  cudaCheck(cudaEventDestroy(ratioSyncOneBodyEvent));
  cudaCheck(cudaEventDestroy(ratioSyncTwoBodyEvent));
}

void finalizeCublas() { cublasDestroy(cublasHandle); }

void synchronize() { cudaCheck(cudaDeviceSynchronize()); }

void streamsSynchronize() { cudaCheck(cudaEventRecord(syncEvent, 0)); }

} // namespace gpu
