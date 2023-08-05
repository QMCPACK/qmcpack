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
#include "OutputManager.h"
#include "determineDefaultDeviceNum.h"

namespace qmcplusplus
{
#if defined(ENABLE_OFFLOAD)
syclDeviceInfo::syclDeviceInfo(const sycl::context& context, const sycl::device& device, const omp_interop_t& interop)
    : context_(context), device_(device), interop_(interop)
#else
syclDeviceInfo::syclDeviceInfo(const sycl::context& context, const sycl::device& device)
    : context_(context), device_(device)
#endif
{}

syclDeviceInfo::~syclDeviceInfo()
{
#if defined(ENABLE_OFFLOAD)
#pragma omp interop destroy(interop_)
#endif
}

std::unique_ptr<sycl::queue> SYCLDeviceManager::default_device_queue;

SYCLDeviceManager::SYCLDeviceManager(int& default_device_num, int& num_devices, int local_rank, int local_size)
    : sycl_default_device_num(-1)
{
#if defined(ENABLE_OFFLOAD)
  const size_t omp_num_devices = omp_get_num_devices();
  visible_devices.reserve(omp_num_devices);
  for (int id = 0; id < omp_num_devices; id++)
  {
    omp_interop_t interop;
#pragma omp interop device(id) init(prefer_type("level_zero"), targetsync : interop)

    int err = -1;
    const std::string omp_backend_name(omp_get_interop_str(interop, omp_ipr_fr_name, &err));
    if (err != omp_irc_success)
      throw std::runtime_error("omp_get_interop_str(omp_ipr_fr_name) failed!");

    if (omp_backend_name.find("level_zero") != 0)
      throw std::runtime_error("Interop between OpenMP and SYCL is only supported when both implementations are built "
                               "on top of Level Zero API.");

    auto hPlatform = omp_get_interop_ptr(interop, omp_ipr_platform, &err);
    if (err != omp_irc_success)
      throw std::runtime_error("omp_get_interop_ptr(omp_ipr_platform) failed!");
    auto hContext = omp_get_interop_ptr(interop, omp_ipr_device_context, &err);
    if (err != omp_irc_success)
      throw std::runtime_error("omp_get_interop_ptr(omp_ipr_device_context) failed!");
    auto hDevice = omp_get_interop_ptr(interop, omp_ipr_device, &err);
    if (err != omp_irc_success)
      throw std::runtime_error("omp_get_interop_ptr(omp_ipr_device) failed!");

    const sycl::platform sycl_platform =
        sycl::ext::oneapi::level_zero::make_platform(reinterpret_cast<pi_native_handle>(hPlatform));

    const sycl::device sycl_device =
        sycl::ext::oneapi::level_zero::make_device(sycl_platform, reinterpret_cast<pi_native_handle>(hDevice));

    visible_devices
        .emplace_back(sycl::ext::oneapi::level_zero::make_context({sycl_device},
                                                                  reinterpret_cast<pi_native_handle>(hContext),
                                                                  true /* keep the ownership, no transfer */),
                      sycl_device, interop);
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
    visible_devices.emplace_back(sycl::context(devices[id]), devices[id]);
#endif

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
    default_device_queue = std::make_unique<sycl::queue>(visible_devices[sycl_default_device_num].get_context(),
                                                         visible_devices[sycl_default_device_num].get_device(),
                                                         sycl::property::queue::in_order());
    if (!visible_devices[sycl_default_device_num].get_device().has(sycl::aspect::ext_intel_free_memory))
      app_warning()
          << "Free memory queries always return 0 due to inactive 'oneAPI' System Resource Management (sysman). "
          << "Set environment variable ZES_ENABLE_SYSMAN to 1 to activate the query feature." << std::endl;
  }
}

sycl::queue& SYCLDeviceManager::getDefaultDeviceDefaultQueue()
{
  if (!default_device_queue)
    throw std::runtime_error("SYCLDeviceManager::getDefaultDeviceQueue() the global instance not initialized.");
  return *default_device_queue;
}

sycl::queue SYCLDeviceManager::createQueueDefaultDevice() const
{
  return sycl::queue(visible_devices[sycl_default_device_num].get_context(),
                     visible_devices[sycl_default_device_num].get_device(), sycl::property::queue::in_order());
}

} // namespace qmcplusplus
