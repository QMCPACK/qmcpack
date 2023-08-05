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

#include <sycl/sycl.hpp>
#include "SYCLDeviceManager.h"
#include "SYCLruntime.hpp"

namespace qmcplusplus
{
sycl::queue& getSYCLDefaultDeviceDefaultQueue() { return SYCLDeviceManager::getDefaultDeviceDefaultQueue(); }
size_t getSYCLdeviceFreeMem()
{
  auto device = getSYCLDefaultDeviceDefaultQueue().get_device();
  if (device.has(sycl::aspect::ext_intel_free_memory))
    return getSYCLDefaultDeviceDefaultQueue().get_device().get_info<sycl::ext::intel::info::device::free_memory>();
  else
    return 0;
}
} // namespace qmcplusplus
