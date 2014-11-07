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
