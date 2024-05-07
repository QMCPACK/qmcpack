//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_QUEUE_CUDA_H
#define QMCPLUSPLUS_QUEUE_CUDA_H

#include "Queue.hpp"
#include "CUDAruntime.hpp"

namespace qmcplusplus
{

namespace compute
{

template<>
class Queue<PlatformKind::CUDA>
{
public:
  Queue() { cudaErrorCheck(cudaStreamCreate(&hstream_), "cudaStreamCreate failed!"); }

  ~Queue() { cudaErrorCheck(cudaStreamDestroy(hstream_), "cudaStreamDestroy failed!"); }

  // dualspace container
  template<class DSC>
  void enqueueH2D(DSC& dataset, typename DSC::size_type size = 0, typename DSC::size_type offset = 0)
  {
    if (dataset.data() == dataset.device_data()) return;

    if (size == 0)
    {
      cudaErrorCheck(cudaMemcpyAsync(dataset.device_data() + offset, dataset.data() + offset,
                                     dataset.size() * sizeof(typename DSC::value_type), cudaMemcpyHostToDevice, hstream_),
                     "Queue<PlatformKind::CUDA>::enqueueH2D cudaMemcpyAsync failed!");
    }
    else
    {
      cudaErrorCheck(cudaMemcpyAsync(dataset.device_data() + offset, dataset.data() + offset,
                                     size * sizeof(typename DSC::value_type), cudaMemcpyHostToDevice, hstream_),
                     "Queue<PlatformKind::CUDA>::enqueueH2D cudaMemcpyAsync failed!");
    }
  }

  template<class DSC>
  void enqueueD2H(DSC& dataset, typename DSC::size_type size = 0, typename DSC::size_type offset = 0)
  {
    if (dataset.data() == dataset.device_data()) return;

    if (size == 0)
    {
      cudaErrorCheck(cudaMemcpyAsync(dataset.data() + offset, dataset.device_data() + offset,
                                     dataset.size() * sizeof(typename DSC::value_type), cudaMemcpyDeviceToHost, hstream_),
                     "Queue<PlatformKind::CUDA>::enqueueD2H cudaMemcpyAsync failed!");
    }
    else
    {
      cudaErrorCheck(cudaMemcpyAsync(dataset.data() + offset, dataset.device_data() + offset,
                                     size * sizeof(typename DSC::value_type), cudaMemcpyDeviceToHost, hstream_),
                     "Queue<PlatformKind::CUDA>::enqueueD2H cudaMemcpyAsync failed!");
    }
  }

  void sync() { cudaErrorCheck(cudaStreamSynchronize(hstream_), "cudaStreamSynchronize failed!"); }

  cudaStream_t getNative() { return hstream_; }

private:
  cudaStream_t hstream_;
};

} // namespace compute

} // namespace qmcplusplus

#endif
