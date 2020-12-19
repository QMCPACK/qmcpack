#ifndef QMCPLUSPLUS_SCOPEDPROFILER_H
#define QMCPLUSPLUS_SCOPEDPROFILER_H

#include "config.h"
#ifdef ENABLE_CUDA
#include "CUDA/cudaError.h"
#include <cuda_profiler_api.h>
#endif
#ifdef USE_VTUNE_API
#include <ittnotify.h>
#endif

namespace qmcplusplus
{

//non-nested
class ScopedProfiler
{
public:
  ScopedProfiler(bool active) : isActive(active)
  {
    if (isActive)
    {
#ifdef ENABLE_CUDA
      cudaErrorCheck(cudaProfilerStart(), "cudaProfilerStart failed!");
#endif
#ifdef USE_VTUNE_API
      __itt_resume();
#endif
    }
  }

  ~ScopedProfiler()
  {
    if (isActive)
    {
#ifdef ENABLE_CUDA
      cudaErrorCheck(cudaProfilerStop(), "cudaProfilerStop failed!");
#endif
#ifdef USE_VTUNE_API
      __itt_pause();
#endif
    }
  }
private:
  const bool isActive;
};
}
#endif
