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


#ifndef QMCPLUSPLUS_OMPDEVICEMANAGER_H
#define QMCPLUSPLUS_OMPDEVICEMANAGER_H

namespace qmcplusplus
{

/** OpenMP device manager
 */
class OMPDeviceManager
{
  int omp_default_device_num;
  const int omp_device_count;

public:
  OMPDeviceManager(int& default_device_num, int& num_devices, int local_rank, int local_size);
};
} // namespace qmcplusplus

#endif
