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
#include "OutputManager.h"
#include "determineDefaultDeviceNum.h"
#if defined(_OPENMP)
#include <omp.h>
#endif

namespace qmcplusplus
{
#if defined(_OPENMP)
/** create SYCL device/contexts from OpenMP owned ones to ensure interoperability.
 * CUDA has the notion of primary context while SYCL requires explicitly sharing context.
 */
static std::vector<struct syclDeviceInfo> xomp_get_sycl_devices();
#endif

SYCLDeviceManager::SYCLDeviceManager(int& default_device_num, int& num_devices, int local_rank, int local_size)
    : sycl_default_device_num(-1)
{
  // Potentially multiple GPU platform.
  std::vector<sycl::platform> platforms = sycl::platform::get_platforms();
  if (platforms.empty())
    throw std::runtime_error("Cannot find SYCL platforms!");

    // find out devices from the first platform with GPUs.
#if defined(ENABLE_OFFLOAD)
  visible_devices = xomp_get_sycl_devices();
#else
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

    default_device_queue = std::make_unique<sycl::queue>(visible_devices[sycl_default_device_num].context,
                                                         visible_devices[sycl_default_device_num].device);
  }
}

std::unique_ptr<sycl::queue> SYCLDeviceManager::default_device_queue;

sycl::queue& SYCLDeviceManager::getDefaultDeviceQueue()
{
  if (!default_device_queue)
    throw std::runtime_error("SYCLDeviceManager::getDefaultDeviceQueue() the global instance not initialized.");
  return *default_device_queue;
}

#if defined(_OPENMP)
static std::vector<struct syclDeviceInfo> xomp_get_sycl_devices()
{
  enum class Backend
  {
    UNKNOWN,
    LEVEL_ZERO,
    OPENCL
  };

  const auto num_omp_devices = omp_get_num_devices();
  Backend selected_backend   = Backend::UNKNOWN;
  std::vector<struct syclDeviceInfo> devices(num_omp_devices);
  for (int id = 0; id < num_omp_devices; id++)
  {
    omp_interop_t o    = 0;
    Backend my_backend = Backend::UNKNOWN;
#pragma omp interop init(prefer_type("sycl"), targetsync : o) device(id)
    int err = -1;

    const std::string omp_backend(omp_get_interop_str(o, omp_ipr_fr_name, &err));
    assert(err >= 0 && "omp_get_interop_str(omp_ipr_fr_name)");

    if (omp_backend.find("level_zero") == 0)
    {
      my_backend = Backend::LEVEL_ZERO;

      auto hPlatform = omp_get_interop_ptr(o, omp_ipr_platform, &err);
      assert(err >= 0 && "omp_get_interop_ptr(omp_ipr_platform)");
      auto hContext = omp_get_interop_ptr(o, omp_ipr_device_context, &err);
      assert(err >= 0 && "omp_get_interop_ptr(omp_ipr_device_context)");
      auto hDevice = omp_get_interop_ptr(o, omp_ipr_device, &err);
      assert(err >= 0 && "omp_get_interop_ptr(omp_ipr_device)");

      const sycl::platform sycl_platform =
          sycl::ext::oneapi::level_zero::make_platform(reinterpret_cast<pi_native_handle>(hPlatform));
      devices[id].device =
          sycl::ext::oneapi::level_zero::make_device(sycl_platform, reinterpret_cast<pi_native_handle>(hDevice));

      devices[id].context = sycl::ext::oneapi::level_zero::make_context({devices[id].device},
                                                                        reinterpret_cast<pi_native_handle>(hContext),
                                                                        true /* keep the ownership, no transfer */);
    }
    else if (omp_backend.find("opencl") == 0)
    {
      my_backend = Backend::OPENCL;
      /*
                auto hContext = omp_get_interop_ptr(o, omp_ipr_device_context, &err);
                assert (err >= 0 && "omp_get_interop_ptr(omp_ipr_device_context)");
                auto hDevice =  omp_get_interop_ptr(o, omp_ipr_device, &err);
                assert (err >= 0 && "omp_get_interop_ptr(omp_ipr_device)");
          devices[id].device = sycl::make_device<sycl::backend::opencl>(static_cast<cl_device>(hDevice));
          devices[id].context = sycl::make_context<sycl::backend::opencl>(static_cast<cl_context>(hContext));
*/
    }

#pragma omp interop destroy(o)

    if (selected_backend == Backend::UNKNOWN)
      selected_backend = my_backend;
    else if (selected_backend != my_backend)
      throw std::runtime_error("Inconsistent backends detected among OpenMP devices.");
  }

  if (devices.size() > 0)
    switch (selected_backend)
    {
    case Backend::LEVEL_ZERO:
      app_log() << "SYCL adopts the Level Zero backend chosen by OpenMP." << std::endl;
      break;
    case Backend::OPENCL:
      throw std::runtime_error("OpenMP has chosen the OpenCL backend. "
                               "We have not yet worked out its interoperability with SYCL. "
                               "Please set the Level Zero backend in OpenMP!");
      break;
    default:
      throw std::runtime_error("Failed in extracting OpenMP backend supported by SYCL.");
    }

  return devices;
}
#endif

} // namespace qmcplusplus
