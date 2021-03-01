//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "config.h"
#ifdef ENABLE_CUDA
#include <cuda_runtime_api.h>
#include "CUDA/cudaError.h"
#endif
#include <omp.h>
#include "Utilities/qmc_common.h"
#include "accelerators.hpp"
#include "Platforms/Host/OutputManager.h"

namespace qmcplusplus
{
/** distribute MPI ranks among devices
 *
 * the amount of MPI ranks for each device differs by 1 at maximum.
 * larger id has more MPI ranks.
 */
int getDeviceID(int rank_id, int num_ranks, int num_devices)
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

int assignAccelerators(Communicate& NodeComm)
{
  int num_accelerators(0);
  int assigned_accelerators_id(-1);
#ifdef ENABLE_CUDA
  int cudaDeviceCount;
  int cudaDeviceID;
  cudaErrorCheck(cudaGetDeviceCount(&cudaDeviceCount), "cudaGetDeviceCount failed!");
  if (num_accelerators == 0)
    num_accelerators = cudaDeviceCount;
  else if (num_accelerators != cudaDeviceCount)
    throw std::runtime_error("Inconsistent number of CUDA devices with the previous record!");
  if (cudaDeviceCount > NodeComm.size())
    app_warning() << "More CUDA devices than the number of MPI ranks. "
                  << "Some devices will be left idle.\n"
                  << "There is potential performance issue with the GPU affinity. "
                  << "Use CUDA_VISIBLE_DEVICE or MPI launcher to expose desired devices.\n";
  if (num_accelerators > 0)
  {
    cudaDeviceID = getDeviceID(NodeComm.rank(), NodeComm.size(), cudaDeviceCount);
    if (assigned_accelerators_id < 0)
      assigned_accelerators_id = cudaDeviceID;
    else if (assigned_accelerators_id != cudaDeviceID)
      throw std::runtime_error("Inconsistent assigned CUDA devices with the previous record!");
#pragma omp parallel
    {
      cudaErrorCheck(cudaSetDevice(cudaDeviceID), "cudaSetDevice failed!");
      cudaErrorCheck(cudaFree(0), "cudaFree failed!");
    }
  }
#endif
#ifdef ENABLE_OFFLOAD
  int ompDeviceCount = omp_get_num_devices();
  int ompDeviceID;
  if (num_accelerators == 0)
    num_accelerators = ompDeviceCount;
  else if (num_accelerators != ompDeviceCount)
    throw std::runtime_error("Inconsistent number of OpenMP devices with the previous record!");
  if (ompDeviceCount > NodeComm.size())
    app_warning() << "More OpenMP devices than the number of MPI ranks. "
                  << "Some devices will be left idle.\n"
                  << "There is potential performance issue with the GPU affinity.\n";
  if (num_accelerators > 0)
  {
    ompDeviceID = getDeviceID(NodeComm.rank(), NodeComm.size(), ompDeviceCount);
    if (assigned_accelerators_id < 0)
      assigned_accelerators_id = ompDeviceID;
    else if (assigned_accelerators_id != ompDeviceID)
      throw std::runtime_error("Inconsistent assigned OpenMP devices with the previous record!");
    omp_set_default_device(ompDeviceID);
  }
#endif

  if (num_accelerators > 0)
  {
    if (NodeComm.size() % num_accelerators != 0)
      app_warning() << "The number of MPI ranks on the node is not divisible by the number of accelerators. "
                    << "Imbalanced load may cause performance loss.\n";
    qmc_common.default_accelerator_id = assigned_accelerators_id;
  }

  return num_accelerators;
}

} // namespace qmcplusplus
