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
