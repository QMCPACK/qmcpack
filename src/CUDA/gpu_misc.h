#ifndef GPU_MISC_H
#define GPU_MISC_H

#include <cstdlib>
#include <cstdio>

#include <cuda_runtime_api.h>
#include <vector>

namespace gpu
{

extern cudaStream_t kernelStream;
extern cudaStream_t memoryStream;

extern cudaEvent_t syncEvent;

extern cudaEvent_t gradientSyncDiracEvent;
extern cudaEvent_t gradientSyncOneBodyEvent;
extern cudaEvent_t gradientSyncTwoBodyEvent;

extern cudaEvent_t ratioSyncDiracEvent;
extern cudaEvent_t ratioSyncOneBodyEvent;
extern cudaEvent_t ratioSyncTwoBodyEvent;



void initCUDAStreams();
void initCUDAEvents();

void synchronize();

void streamsSynchronize();

}
#endif

