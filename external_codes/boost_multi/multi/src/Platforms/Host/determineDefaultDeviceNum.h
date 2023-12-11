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


#ifndef QMCPLUSPLUS_DETERMINEDEFAULTDEVICENUM_H
#define QMCPLUSPLUS_DETERMINEDEFAULTDEVICENUM_H

namespace qmcplusplus
{
/** distribute MPI ranks among devices
 *
 * the amount of MPI ranks for each device differs by 1 at maximum.
 * larger id has more MPI ranks.
 */
inline int determineDefaultDeviceNum(int num_devices, int rank_id, int num_ranks)
{
  if (num_ranks < num_devices)
    num_devices = num_ranks;
  // ranks are equally distributed among devices
  int min_ranks_per_device = num_ranks / num_devices;
  int residual             = num_ranks % num_devices;
  int assigned_device_id;
  if (rank_id < min_ranks_per_device * (num_devices - residual))
    assigned_device_id = rank_id / min_ranks_per_device;
  else
    assigned_device_id = (rank_id + num_devices - residual) / (min_ranks_per_device + 1);
  return assigned_device_id;
}
}

#endif
