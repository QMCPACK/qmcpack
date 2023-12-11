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

#include <vector>
#include <memory>
#include <sycl/sycl.hpp>
#include "config.h"
#if defined(ENABLE_OFFLOAD)
#include <omp.h>
#endif

namespace qmcplusplus
{
class syclDeviceInfo
{
public:
#if defined(ENABLE_OFFLOAD)
  syclDeviceInfo(const sycl::context& context, const sycl::device& device, const omp_interop_t& interop);
#else
  syclDeviceInfo(const sycl::context& context, const sycl::device& device);
#endif
  ~syclDeviceInfo();
  const sycl::context& get_context() const { return context_; }
  const sycl::device& get_device() const { return device_; }

private:
  sycl::context context_;
  sycl::device device_;
#if defined(ENABLE_OFFLOAD)
  omp_interop_t interop_;
#endif
};

/** SYCL device manager
 */
class SYCLDeviceManager
{
  int sycl_default_device_num;
  std::vector<syclDeviceInfo> visible_devices;

  /// the global singleton which can be used to access the default queue of the default device.
  static std::unique_ptr<sycl::queue> default_device_queue;

public:
  SYCLDeviceManager(int& default_device_num, int& num_devices, int local_rank, int local_size);

  /** access the the DeviceManager owned default queue.
   * Restrict the use of it to performance non-critical operations.
   * Note: CUDA has a default queue but all the SYCL queues are explicit.
   */
  static sycl::queue& getDefaultDeviceDefaultQueue();
  sycl::queue createQueueDefaultDevice() const;
};
} // namespace qmcplusplus

#endif
