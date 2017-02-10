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

/** Obtains the number of appropriate Cuda devices on the current MPI rank's node
 */
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

/** First determines the current MPI rank's relative number on its node,
 *  then returns which Cuda device to use.
 */
inline int get_device_num()
{
  const int MAX_LEN = 200;
  int size = OHMMS::Controller->size(); // how many MPI ranks
  int rank = OHMMS::Controller->rank(); // current MPI rank number
  std::vector<char> myname(MAX_LEN);
  gethostname(&myname[0], MAX_LEN);
  std::vector<char> host_list(MAX_LEN*size);
  // copy current MPI rank's hostname to (shared) host_list
  for (int i=0; i<MAX_LEN; i++)
    host_list[rank*MAX_LEN+i] = myname[i];
  OHMMS::Controller->allgather(myname, host_list, MAX_LEN); // wait till every rank finishes
  // copy host_list to hostnames vector
  std::vector<std::string> hostnames;
  for (int i=0; i<size; i++)
    hostnames.push_back(&(host_list[i*MAX_LEN]));
  std::string myhostname = &myname[0];
  // calculate how many MPI ranks are ahead of the current one on the current node
  int number_ranks = 0;
  int curr_ranknum = 0;
  for (int i=0; i<size; i++) // loop over all ranks
  {
    if (hostnames[i] == myhostname)
    {
      number_ranks++; // count all ranks with the same host name as the current one
      if (i<rank)
        curr_ranknum++; // count only the ones schedule before the current one
    }
  }
  int num_cuda_devices=get_num_appropriate_devices();
  if (curr_ranknum==0) // output some information if we're the first rank on a node
  {
    std::cerr << number_ranks << " MPI ranks on node " << myhostname << " with " << num_cuda_devices << " appropriate CUDA devices." << std::endl;
    if(number_ranks<num_cuda_devices)
    {
      std::cerr << "WARNING: Less MPI ranks than Cuda devices. Some Cuda devices (device # >= " << number_ranks << ") will not be used." << std::endl;
    }
    else
    {
      if(number_ranks%num_cuda_devices) // is only true (>0) when number of MPI ranks is not a multiple of Cuda device number
        std::cerr << "WARNING: Number of MPI ranks is not a multiple of the number of Cuda devices." << std::endl;
    }
  }
  // return Cuda device number based on how many appropriate ones exist on the current rank's node
  int devnum=curr_ranknum % num_cuda_devices;
  return devnum;
}

/** Sets the Cuda device of the current MPI rank from the pool of appropriate Cuda devices
 * @param num requested Cuda device number
 */
inline void set_appropriate_device_num(int num)
{
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  int num_appropriate=0, device=0;
  bool set_cuda_device=false;
  for (device = 0; device < deviceCount; ++device)
  {
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, device);
    if (((deviceProp.major >= 1) && (deviceProp.minor >= 3)) ||
        deviceProp.major >= 2)
    {
      num_appropriate++;
      if (num_appropriate == num+1)
      {
        cudaSetDevice (device);
        set_cuda_device=true;
        std::cerr << "<- Rank " << OHMMS::Controller->rank() << " has acquired CUDA device #" << device << std::endl;
        break; // the device is set, nothing more to do here
      }
    }
  }
  // This is a fail-safe that should never be triggered
  if(!set_cuda_device)
  {
    APP_ABORT("Failure to obtain requested CUDA device.");
  }
}

inline void Finalize_CUDA()
{
  gpu::finalizeCublas();
  gpu::finalizeCUDAEvents();
  gpu::finalizeCUDAStreams();
  cudaDeviceReset();
}

/** Initialize Cuda device on current MPI rank
 */
inline void Init_CUDA()
{
  int devNum = get_device_num();
  set_appropriate_device_num(devNum);
  cudaDeviceSetLimit(cudaLimitPrintfFifoSize, 1024 * 1024 * 50);
  gpu::rank=OHMMS::Controller->rank();
  gpu::initCUDAStreams();
  gpu::initCUDAEvents();
  gpu::initCublas();
  gpu::MaxGPUSpineSizeMB = MAX_GPU_SPLINE_SIZE_MB;
  // Output maximum spline buffer size for first MPI rank
  if(gpu::rank==0) std::cerr << "Default MAX_GPU_SPLINE_SIZE_MB is " << gpu::MaxGPUSpineSizeMB << " MB." << std::endl;
  return;
}
#else
inline void Init_CUDA()
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
