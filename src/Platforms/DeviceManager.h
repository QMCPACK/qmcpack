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
#if defined(ENABLE_SYCL)
#include "SYCL/SYCLDeviceManager.h"
#endif

namespace qmcplusplus
{

/** orchestrate platforms (programming models) and set up the default device on each platform
 *
 * Currently support CUDA, OpenMP, SYCL
 *
 * Each platform provides its own device manager class. Each individual device manager
 * initializes all the devices visible to the process and also set the default device.
 * The first initialized platform selects the default device for the current process and all the rest
 * platforms follow that choice.
 *
 * The device manager on each platform may be constomized to store a list of known device numbers or
 * a list of context device handles depending on the need.
 *
 * DeviceManager assumes there is only one type of accelerators although they may support multiple
 * platforms. Under this assumption, the numbers of devices captured on all the platforms must agree.
 *
 * DeviceManager::global is intended to the per-process global instance and should be initialized
 * after MPI initialization in main().
 */
class DeviceManager
{
public:
  /** constructor
   * @param local_rank the rank id of the node local MPI communicator
   * @param local_size the size of the node local MPI communicator
   */
  DeviceManager(int local_rank, int local_size);

  ~DeviceManager();

  // accessors
  int getDefaultDeviceNum() const { return default_device_num; }
  int getNumDevices() const { return num_devices; }

#if defined(ENABLE_CUDA)
  const auto& getCUDADM() const { return cuda_dm_; }
#endif
#if defined(ENABLE_OFFLOAD)
  const auto& getOMPDM() const { return omptarget_dm_; }
#endif
#if defined(ENABLE_SYCL)
  const auto& getSYCLDM() const { return sycl_dm_; }
#endif

  /** initialize the global instance of DeviceManager
   * arguments are the same as the constructor
   */
  static void initializeGlobalDeviceManager(int local_rank, int local_size);
  /// global instance accessor
  static const DeviceManager& getGlobal();

private:
  /// the global singleton which can be used to query default devices and all the captured devices.
  static std::unique_ptr<DeviceManager> global;
  /// the id of default device. Must be defined before platform device manager objects
  int default_device_num;
  /// the number of devices. Must be defined before platform device manager objects
  int num_devices;
#if defined(ENABLE_CUDA)
  /// CUDA device manager object
  CUDADeviceManager cuda_dm_;
#endif
#if defined(ENABLE_OFFLOAD)
  /// OpenMP device manager object
  OMPDeviceManager omptarget_dm_;
#endif
#if defined(ENABLE_SYCL)
  /// SYCL device manager object
  SYCLDeviceManager sycl_dm_;
#endif
};

} // namespace qmcplusplus

#endif
