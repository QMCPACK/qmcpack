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

void
initCUDAStreams()
{
  cudaStreamCreate(&kernelStream);
  cudaStreamCreate(&memoryStream);
}

void
initCUDAEvents()
{
  cudaEventCreateWithFlags(&syncEvent, cudaEventDisableTiming);
  cudaEventCreateWithFlags(&gradientSyncDiracEvent, cudaEventDisableTiming);
  cudaEventCreateWithFlags(&gradientSyncOneBodyEvent, cudaEventDisableTiming);
  cudaEventCreateWithFlags(&gradientSyncTwoBodyEvent, cudaEventDisableTiming);
  cudaEventCreateWithFlags(&ratioSyncDiracEvent, cudaEventDisableTiming);
  cudaEventCreateWithFlags(&ratioSyncOneBodyEvent, cudaEventDisableTiming);
  cudaEventCreateWithFlags(&ratioSyncTwoBodyEvent, cudaEventDisableTiming);
}

void
initCublas()
{
  cublasCreate(&cublasHandle);
}

void
finalizeCUDAStreams()
{
  cudaStreamDestroy(kernelStream);
  cudaStreamDestroy(memoryStream);
}

void
finalizeCUDAEvents()
{
  cudaEventDestroy(syncEvent);
  cudaEventDestroy(gradientSyncDiracEvent);
  cudaEventDestroy(gradientSyncOneBodyEvent);
  cudaEventDestroy(gradientSyncTwoBodyEvent);
  cudaEventDestroy(ratioSyncDiracEvent);
  cudaEventDestroy(ratioSyncOneBodyEvent);
  cudaEventDestroy(ratioSyncTwoBodyEvent);
}

void
finalizeCublas()
{
  cublasDestroy(cublasHandle);
}

void
synchronize()
{
  cudaDeviceSynchronize();
}

void
streamsSynchronize()
{
  cudaEventRecord(syncEvent, 0);
}

}
