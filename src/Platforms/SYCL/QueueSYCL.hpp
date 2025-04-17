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

#ifndef QMCPLUSPLUS_QUEUE_SYCL_H
#define QMCPLUSPLUS_QUEUE_SYCL_H

#include "Common/Queue.hpp"
#include "SYCLruntime.hpp"

namespace qmcplusplus
{

namespace compute
{

template<>
class Queue<PlatformKind::SYCL> : public QueueBase
{
public:
  Queue() : queue_(createSYCLInOrderQueueOnDefaultDevice()) {}

  template<class T>
  inline void memcpy(T* dst, const T* src, size_t size)
  {
    queue_.memcpy(dst, src, sizeof(T) * size);
  }

  // dualspace container
  template<class DSC>
  void enqueueH2D(DSC& dataset, typename DSC::size_type size = 0, typename DSC::size_type offset = 0)
  {
    if (dataset.data() == dataset.device_data())
      return;
    memcpy(dataset.device_data() + offset, dataset.data() + offset, size ? size : dataset.size());
  }

  template<class DSC>
  void enqueueD2H(DSC& dataset, typename DSC::size_type size = 0, typename DSC::size_type offset = 0)
  {
    if (dataset.data() == dataset.device_data())
      return;

    memcpy(dataset.data() + offset, dataset.device_data() + offset, size ? size : dataset.size());
  }

  void sync() { queue_.wait(); }

  sycl::queue& getNative() { return queue_; }

private:
  sycl::queue queue_;
};

} // namespace compute

} // namespace qmcplusplus

#endif
