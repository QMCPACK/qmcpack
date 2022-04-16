//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Thomas Applencourt, apl@anl.gov, Argonne National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Thomas Applencourt, apl@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SYCLDEVICEMANAGER_H
#define QMCPLUSPLUS_SYCLDEVICEMANAGER_H

#include <CL/sycl.hpp>
#include <vector>
#include "config.h"

namespace qmcplusplus
{

struct syclDeviceInfo
{
  sycl::context context;
  sycl::device device;
};

/** SYCL device manager
 */
class SYCLDeviceManager
{
  int sycl_default_device_num;
  std::vector<syclDeviceInfo> visible_devices;
  sycl::queue default_device_queue;

public:
  SYCLDeviceManager(int& default_device_num, int& num_devices, int local_rank, int local_size);

  /** access the the DeviceManager owned default queue.
   * Restrict the use of it to performance non-critical operations.
   * Note: CUDA has a default queue but all the SYCL queues are explicit.
   * Right now we return a copy of the default queue. Queues hold contexts and devices by referece.
   * So making a copy is expected to be cheap. If this is not the case, we will find a cheap solution.
   */
  sycl::queue getDefaultDeviceQueue() const { return default_device_queue; }
};
} // namespace qmcplusplus

#endif
