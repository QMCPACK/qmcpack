//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Thomas Applencourt, apl@anl.gov, Argonne National Laboratory
//
// File created by: Thomas Applencourt, apl@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SYCLDEVICEMANAGER_H
#define QMCPLUSPLUS_SYCLDEVICEMANAGER_H

#include <stdexcept>
#include <CL/sycl.hpp>
#include "Host/OutputManager.h"
#include "determineDefaultDeviceNum.h"

namespace qmcplusplus
{

struct syclDeviceInfo {
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
  SYCLDeviceManager(int& default_device_num, int& num_devices, int local_rank, int local_size)
      : sycl_default_device_num(-1)
  {
    // Potentially multiple GPU platform.
    std::vector<sycl::platform> platforms = sycl::platform::get_platforms();
    if (platforms.empty())
      throw std::runtime_error("Cannot find SYCL platforms!");

    // find out devices from the first platform with GPUs.
    std::vector<sycl::device> devices;
    app_log() << "Visible SYCL platforms are :" << std::endl;
    for (auto& platform : platforms)
    {
      std::vector<sycl::device> gpu_devices = platform.get_devices(sycl::info::device_type::gpu);
      const auto gpu_count = gpu_devices.size();
      bool selected = false;
      if (devices.empty() && gpu_count > 0)
      {
        selected = true;
        devices = std::move(gpu_devices);
      }
      app_log() << (selected ? " ** " : "    ") << platform.get_info<sycl::info::platform::name>()
                << " with " << gpu_count << " GPUs." << std::endl;
    }
    app_log() << std::endl;

    const size_t sycl_device_count = devices.size();
    if (num_devices == 0)
      num_devices = sycl_device_count;
    else if (num_devices != sycl_device_count)
      throw std::runtime_error("Inconsistent number of SYCL devices with the previous record!");
    if (sycl_device_count > local_size)
      app_warning() << "More SYCL devices than the number of MPI ranks. "
                    << "Some devices will be left idle.\n"
                    << "There is potential performance issue with the GPU affinity.\n";
    if (num_devices > 0)
    {
      sycl_default_device_num = determineDefaultDeviceNum(sycl_device_count, local_rank, local_size);
      if (default_device_num < 0)
        default_device_num = sycl_default_device_num;
      else if (default_device_num != sycl_default_device_num)
        throw std::runtime_error("Inconsistent assigned SYCL devices with the previous record!");

      visible_devices.reserve(num_devices);
      for (int id = 0; id < num_devices; id++)
        visible_devices.push_back({sycl::context(devices[sycl_default_device_num]), devices[sycl_default_device_num]});
      default_device_queue = sycl::queue(visible_devices[sycl_default_device_num].context, visible_devices[sycl_default_device_num].device);
    }
  }

  sycl::queue getDefaultDeviceQueue() const { return default_device_queue; }
};
} // namespace qmcplusplus

#endif
