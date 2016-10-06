//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ying Wai Li, yingwaili@ornl.gov, Oak Ridge National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_CUDA_INIT_H
#define QMCPLUSPLUS_CUDA_INIT_H

#ifdef QMC_CUDA
#include "Configuration.h"
#include "Message/Communicate.h"
#include "Utilities/OhmmsInfo.h"
#include "Message/CommOperators.h"
#include <cuda_runtime_api.h>
#include <unistd.h>
#include <CUDA/gpu_misc.h>

#define MAX_GPU_SPLINE_SIZE_MB 81920

inline int get_device_num()
{
  const int MAX_LEN = 200;
  int size = OHMMS::Controller->size();
  int rank = OHMMS::Controller->rank();
  std::vector<char> myname(MAX_LEN);
  gethostname(&myname[0], MAX_LEN);
  std::vector<char> host_list(MAX_LEN*size);
  for (int i=0; i<MAX_LEN; i++)
    host_list[rank*MAX_LEN+i] = myname[i];
  OHMMS::Controller->allgather(myname, host_list, MAX_LEN);
  std::vector<std::string> hostnames;
  for (int i=0; i<size; i++)
    hostnames.push_back(&(host_list[i*MAX_LEN]));
  std::string myhostname = &myname[0];
  int devnum = 0;
  for (int i=0; i<rank; i++)
    if (hostnames[i] == myhostname)
      devnum++;
  return devnum;
}

inline int get_num_appropriate_devices()
{
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  int num_appropriate=0;
  for (int device=0; device < deviceCount; ++device)
  {
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, device);
    if (((deviceProp.major >= 1) && (deviceProp.minor >= 3)) ||
        deviceProp.major >= 2)
      num_appropriate++;
  }
  return num_appropriate;
}

inline void set_appropriate_device_num(int num)
{
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  int num_appropriate=0, device=0;
  for (device = 0; device < deviceCount; ++device)
  {
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, device);
    if (((deviceProp.major >= 1) && (deviceProp.minor >= 3)) ||
        deviceProp.major >= 2)
    {
      num_appropriate++;
      if (num_appropriate == num+1)
        cudaSetDevice (device);
    }
  }
}

inline void Finalize_CUDA()
{
  gpu::finalizeCublas();
  gpu::finalizeCUDAEvents();
  gpu::finalizeCUDAStreams();
  cudaDeviceReset();
}

inline void Init_CUDA(int rank, int size)
{
  int devNum = get_device_num();
  std::ostringstream o;
  o << "Rank = " << rank
       << "  My device number = " << devNum << std::endl;
  std::cerr << o.str();
  int num_appropriate = get_num_appropriate_devices();
  if (devNum >= num_appropriate)
  {
    std::cerr << "Not enough double-precision capable GPUs for MPI rank "
         << rank << ".\n";
    abort();
  }
  set_appropriate_device_num (devNum);
  cudaDeviceSetLimit(cudaLimitPrintfFifoSize, 1024 * 1024 * 50);
  gpu::initCUDAStreams();
  gpu::initCUDAEvents();
  gpu::initCublas();
  gpu::MaxGPUSpineSizeMB = MAX_GPU_SPLINE_SIZE_MB;
  gpu::rank=rank;
  if(!rank) std::cerr << "Default MAX_GPU_SPLINE_SIZE_MB is " << gpu::MaxGPUSpineSizeMB << " MB." << std::endl;
  return;
  int numGPUs;
  cudaGetDeviceCount(&numGPUs);
  std::cerr << "There are " << numGPUs << " GPUs.";
  std::cerr << "Size = " << size << std::endl;
  int chunk = size/numGPUs;
  //  int device_num = rank % numGPUs;
  int device_num = rank * numGPUs / size;
  std::cerr << "My device number is " << device_num << std::endl;
  cudaSetDevice (device_num);
}
#else
inline void Init_CUDA(int rank, int size)
{
  std::cerr << "Flag \"--gpu\" was used, but QMCPACK was built without "
       << "GPU code.\nPlease use cmake -DQMC_CUDA=1.\n";
  APP_ABORT("GPU disabled");
}

inline void Finalize_CUDA()
{
  std::cerr << "Flag \"--gpu\" was used, but QMCPACK was built without "
       << "GPU code.\nPlease use cmake -DQMC_CUDA=1.\n";
  APP_ABORT("GPU disabled");
}
#endif
#endif
