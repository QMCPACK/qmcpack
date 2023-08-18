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

namespace qmcplusplus
{

/** CUDA device manager
 */
class CUDADeviceManager
{
  int cuda_default_device_num;
  int cuda_device_count;

public:
  CUDADeviceManager(int& default_device_num, int& num_devices, int local_rank, int local_size);
};
} // namespace qmcplusplus

#endif
