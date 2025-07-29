//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SYCL_RUNTIME_H
#define QMCPLUSPLUS_SYCL_RUNTIME_H

#include <sycl/sycl.hpp>

namespace qmcplusplus
{
/// return a reference to the per-device default queue
sycl::queue& getSYCLDefaultDeviceDefaultQueue();
/// create an in-order queue using the default device
sycl::queue createSYCLInOrderQueueOnDefaultDevice();
/// create a out-of-order queue using the default device
sycl::queue createSYCLQueueOnDefaultDevice();
} // namespace qmcplusplus

#endif
