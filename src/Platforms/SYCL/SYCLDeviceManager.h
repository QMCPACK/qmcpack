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

#if defined(_OPENMP)
std::vector<struct syclDeviceInfo> xomp_get_sycl_devices();
#endif

/** SYCL device manager
 */
class SYCLDeviceManager
{
  int sycl_default_device_num;
  std::vector<syclDeviceInfo> visible_devices;
  sycl::queue default_device_queue;

public:
  SYCLDeviceManager(int& default_device_num, int& num_devices, int local_rank, int local_size);

  sycl::queue getDefaultDeviceQueue() const { return default_device_queue; }
};
} // namespace qmcplusplus

#endif
