//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////


#include "OMPDeviceManager.h"
#include <stdexcept>
#include <omp.h>
#include "OutputManager.h"
#include "determineDefaultDeviceNum.h"

namespace qmcplusplus
{

OMPDeviceManager::OMPDeviceManager(int& default_device_num, int& num_devices, int local_rank, int local_size)
    : omp_default_device_num(-1), omp_device_count(omp_get_num_devices())
{
  if (num_devices == 0)
    num_devices = omp_device_count;
  else if (num_devices != omp_device_count)
    throw std::runtime_error("Inconsistent number of OpenMP devices with the previous record!");
  if (omp_device_count > local_size)
    app_warning() << "More OpenMP devices than the number of MPI ranks. "
                  << "Some devices will be left idle.\n"
                  << "There is potential performance issue with the GPU affinity.\n";
  if (num_devices > 0)
  {
    omp_default_device_num = determineDefaultDeviceNum(omp_device_count, local_rank, local_size);
    if (default_device_num < 0)
      default_device_num = omp_default_device_num;
    else if (default_device_num != omp_default_device_num)
      throw std::runtime_error("Inconsistent assigned OpenMP devices with the previous record!");
    omp_set_default_device(omp_default_device_num);
  }
}
} // namespace qmcplusplus
