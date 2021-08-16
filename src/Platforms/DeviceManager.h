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


#ifndef QMCPLUSPLUS_DEVICEMANAGER_H
#define QMCPLUSPLUS_DEVICEMANAGER_H

#include <memory>
#include <config.h>
#if defined(ENABLE_CUDA)
#include "CUDA/CUDADeviceManager.h"
#endif
#if defined(ENABLE_OFFLOAD)
#include "OMPTarget/OMPDeviceManager.h"
#endif

namespace qmcplusplus
{

class DeviceManager
{
public:
  DeviceManager(int local_rank, int local_size);

  ~DeviceManager();

  int getDefaultDeviceNum() const { return default_device_num; }
  int getNumDevices() const { return num_devices; }

#if defined(ENABLE_CUDA)
  const auto& getCUDADM() const { return cuda_dm_; }
#endif
#if defined(ENABLE_OFFLOAD)
  const auto& getOMPDM() const { return omptarget_dm_; }
#endif

  static std::unique_ptr<DeviceManager> global;

private:
  int default_device_num;
  int num_devices;
#if defined(ENABLE_CUDA)
  CUDADeviceManager cuda_dm_;
#endif
#if defined(ENABLE_OFFLOAD)
  OMPDeviceManager omptarget_dm_;
#endif
};

} // namespace qmcplusplus

#endif
