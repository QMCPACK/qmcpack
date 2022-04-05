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


#ifndef QMCPLUSPLUS_CUDADEVICEMANAGER_H
#define QMCPLUSPLUS_CUDADEVICEMANAGER_H

#include <stdexcept>
#include "CUDAruntime.hpp"
#include "Host/OutputManager.h"
#include "determineDefaultDeviceNum.h"

namespace qmcplusplus
{

/** CUDA device manager
 */
class CUDADeviceManager
{
  int cuda_default_device_num;
  int cuda_device_count;

public:
  CUDADeviceManager(int& default_device_num, int& num_devices, int local_rank, int local_size)
      : cuda_default_device_num(-1), cuda_device_count(0)
  {
    cudaErrorCheck(cudaGetDeviceCount(&cuda_device_count), "cudaGetDeviceCount failed!");
    if (num_devices == 0)
      num_devices = cuda_device_count;
    else if (num_devices != cuda_device_count)
      throw std::runtime_error("Inconsistent number of CUDA devices with the previous record!");
    if (cuda_device_count > local_size)
      app_warning() << "More CUDA devices than the number of MPI ranks. "
                    << "Some devices will be left idle.\n"
                    << "There is potential performance issue with the GPU affinity. "
                    << "Use CUDA_VISIBLE_DEVICE or MPI launcher to expose desired devices.\n";
    if (num_devices > 0)
    {
      cuda_default_device_num = determineDefaultDeviceNum(cuda_device_count, local_rank, local_size);
      if (default_device_num < 0)
        default_device_num = cuda_default_device_num;
      else if (default_device_num != cuda_default_device_num)
        throw std::runtime_error("Inconsistent assigned CUDA devices with the previous record!");
#pragma omp parallel
      {
        cudaErrorCheck(cudaSetDevice(cuda_default_device_num), "cudaSetDevice failed!");
        cudaErrorCheck(cudaFree(0), "cudaFree failed!");
      }
    }
  }
};
} // namespace qmcplusplus

#endif
