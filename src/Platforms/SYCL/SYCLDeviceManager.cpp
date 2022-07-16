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


#include "SYCLDeviceManager.h"
#include <stdexcept>
#include <string>
#include <algorithm>
#include "config.h"
#include "Platforms/Host/OutputManager.h"
#include "Platforms/Host/determineDefaultDeviceNum.h"
namespace qmcplusplus
{

std::unique_ptr<sycl::queue> SYCLDeviceManager::default_device_queue;

SYCLDeviceManager::SYCLDeviceManager(int& default_device_num, int& num_devices, int local_rank, int local_size)
    : sycl_default_device_num(-1)
{
#if defined(ENABLE_OFFLOAD)
  const int sycl_device_count=omp_get_num_devices();
  sycl_default_device_num = determineDefaultDeviceNum(sycl_device_count, local_rank, local_size);
  if (default_device_num < 0)
    default_device_num = sycl_default_device_num;

  visible_devices.resize(sycl_device_count);
  for (int id = 0; id < sycl_device_count; id++)
  {
    omp_interop_t interop;
#pragma omp interop device(id) init(prefer_type("sycl"),targetsync: interop)
    if(id == sycl_default_device_num)
    {
      int result;
      sycl::queue* omp_queue = static_cast<sycl::queue *>(omp_get_interop_ptr(interop, omp_ipr_targetsync, &result));
      if(result != omp_irc_success)
        throw std::runtime_error("SYCLDeviceManager::SYCLDeviceManager fail to obtain sycl::queue by interop");

      default_device_queue.reset(new sycl::queue(*omp_queue));
      visible_devices[id].queue_ = omp_queue;
    }
    visible_devices[id].interop_ = interop;
  } 

#else
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
    const auto gpu_count                  = gpu_devices.size();
    bool selected                         = false;
    if (devices.empty() && gpu_count > 0)
    {
      selected = true;
      devices  = std::move(gpu_devices);
    }
    app_log() << (selected ? " ** " : "    ") << platform.get_info<sycl::info::platform::name>() << " with "
              << gpu_count << " GPUs." << std::endl;
  }
  app_log() << std::endl;

  visible_devices.reserve(devices.size());
  for (int id = 0; id < devices.size(); id++)
    visible_devices.push_back({sycl::context(devices[id]), devices[id]});

  const size_t sycl_device_count = visible_devices.size();
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
    default_device_queue = std::make_unique<sycl::queue>(visible_devices[sycl_default_device_num].context_,
                                                         visible_devices[sycl_default_device_num].device_);
  }
#endif
}

sycl::queue& SYCLDeviceManager::getDefaultDeviceQueue()
{
  if (!default_device_queue)
    throw std::runtime_error("SYCLDeviceManager::getDefaultDeviceQueue() the global instance not initialized.");
  return *default_device_queue;
}

} // namespace qmcplusplus
