//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////


#include "DeviceManager.h"
#include <memory>
#include <stdexcept>
#include "Host/OutputManager.h"

namespace qmcplusplus
{

DeviceManager::DeviceManager(int local_rank, int local_size)
    : default_device_num(-1),
      num_devices(0)
#if defined(ENABLE_CUDA)
      , cuda_dm_(default_device_num, num_devices, local_rank, local_size)
#endif
#if defined(ENABLE_OFFLOAD)
      , omptarget_dm_(default_device_num, num_devices, local_rank, local_size)
#endif
#if defined(ENABLE_SYCL)
      , sycl_dm_(default_device_num, num_devices, local_rank, local_size)
#endif
{
  if (num_devices > 0)
  {
    if (local_size % num_devices != 0)
      app_warning() << "The number of MPI ranks on the node is not divisible by the number of accelerators. "
                    << "Imbalanced load may cause performance loss.\n";
  }
}

DeviceManager::~DeviceManager() = default;

std::unique_ptr<DeviceManager> DeviceManager::global;

void DeviceManager::initializeGlobalDeviceManager(int local_rank, int local_size)
{
  // throw error on subsequent calls to initializeGlobalDeviceManager
  // if the desired behavior is no-op. std::call_once can be used.
  if (global)
    throw std::runtime_error(
        "DeviceManager::initializeGlobalDeviceManager the global instance cannot be initialized again.");
  global = std::make_unique<DeviceManager>(local_rank, local_size);
}

const DeviceManager& DeviceManager::getGlobal()
{
  if (!global)
    throw std::runtime_error("DeviceManager::getGlobal the global instance was not initialized.");
  return *global;
}
} // namespace qmcplusplus
