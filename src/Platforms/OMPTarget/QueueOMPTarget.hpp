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

#ifndef QMCPLUSPLUS_QUEUE_OMPTARGET_H
#define QMCPLUSPLUS_QUEUE_OMPTARGET_H

#include "Queue.hpp"

namespace qmcplusplus
{

namespace compute
{

template<>
class Queue<PlatformKind::OMPTARGET> : public QueueBase
{
public:
  // dualspace container
  template<class DSC>
  void enqueueH2D(DSC& dataset, typename DSC::size_type size = 0, typename DSC::size_type offset = 0)
  {
    if (dataset.data() == dataset.device_data())
      return;

    auto host_ptr = dataset.data();
    if (size == 0)
    {
      PRAGMA_OFFLOAD("omp target update to(host_ptr[offset:dataset.size()])")
    }
    else
    {
      PRAGMA_OFFLOAD("omp target update to(host_ptr[offset:size])")
    }
  }

  template<class DSC>
  void enqueueD2H(DSC& dataset, typename DSC::size_type size = 0, typename DSC::size_type offset = 0)
  {
    if (dataset.data() == dataset.device_data())
      return;

    auto host_ptr = dataset.data();
    if (size == 0)
    {
      PRAGMA_OFFLOAD("omp target update from(host_ptr[offset:dataset.size()])")
    }
    else
    {
      PRAGMA_OFFLOAD("omp target update from(host_ptr[offset:size])")
    }
  }

  void sync() {}
};

} // namespace compute

} // namespace qmcplusplus

#endif
