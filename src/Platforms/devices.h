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
  OHMMS::Controller->allgather(myname, host_list, MAX_LEN); // gathering data from every rank
  // copy host_list to hostnames vector
  std::vector<std::string> hostnames;
  for (int i=0; i<size; i++)
    hostnames.push_back(&(host_list[i*MAX_LEN]));
  std::string myhostname = &myname[0];

  // calculate how many MPI ranks there on the current node, how many are ahead of the current one, and how many nodes are used in total
  std::vector<int> number_ranks(1);
  std::vector<int> cuda_devices(1);
  std::vector<int> rank_node(1);

  number_ranks[0] = 0;
  cuda_devices[0]=get_num_appropriate_devices();

  int relative_ranknum = 0;
  int num_nodes=0;
  std::string curr_host;
  for (int i=0; i<size; i++)
  {
    // counts the number of different hostnames (in other words, how many nodes there are)
    if (hostnames[i] != curr_host)
    {
      curr_host = hostnames[i];
      num_nodes++;
    }
    if (hostnames[i] == myhostname)
    {
      number_ranks[0]++; // count all ranks with the same host name as the current one
      if (i<rank)
        relative_ranknum++; // count only the ones scheduled before the current one
    }
    if (i==rank) rank_node[0]=num_nodes; // node number of current rank (NOTE: node numbers start at 1)
  }

  // gather all the information
  std::vector<int> ranks_per_node(size+1); ranks_per_node[size]=-1;
  OHMMS::Controller->allgather(number_ranks,ranks_per_node,1);
  std::vector<int> num_cuda_devices(size+1); num_cuda_devices[size]=-1;
  OHMMS::Controller->allgather(cuda_devices,num_cuda_devices,1);
  std::vector<int> node_of_rank(size+1); node_of_rank[size]=num_nodes;
  OHMMS::Controller->allgather(rank_node,node_of_rank,1);

  // output information for every rank with a different configuration from the previous one (i.e. with all nodes equal this will only be the first rank)
  if ((ranks_per_node[rank] != ranks_per_node[(rank+size)%(size+1)]) || (num_cuda_devices[rank] != num_cuda_devices[(rank+size)%(size+1)]))
  {
    std::ostringstream out;
    out << ranks_per_node[rank] << " MPI ranks on node(s) " << hostnames[rank];
    // things get a tiny bit more complicated now that we want to output the other hostnames for which the current configuration information is true
    int r=rank;
    int curr_node=node_of_rank[rank];
    // loop over successive ranks with same node configuration
    while ((ranks_per_node[rank] == ranks_per_node[r]) && (num_cuda_devices[rank] == num_cuda_devices[r]))
    {
      if (node_of_rank[r] != curr_node) // when the node number changes, output the new node's host name
      {
        out << ", " << hostnames[r];
        curr_node=node_of_rank[r];
      }
      r++;
      if (r>=size) break; // safety first
    }
    out << " with " << num_cuda_devices[rank] << " appropriate CUDA devices." << std::endl;
    // Output sanity check information for the user
    if(ranks_per_node[rank]<num_cuda_devices[rank])
    {
      out << "WARNING: Fewer MPI ranks than CUDA devices (" << num_cuda_devices[rank] << "). Some CUDA devices (device # >= " << ranks_per_node[rank] << ") will not be used." << std::endl;
    }
    else
    {
      if(ranks_per_node[rank]%num_cuda_devices[rank]) // is only true (>0) when number of MPI ranks is not a multiple of Cuda device number
        out << "WARNING: Number of MPI ranks is not a multiple of the number of CUDA devices (" << num_cuda_devices[rank] << ")." << std::endl;
    }
    std::cerr << out.str();
    std::cerr.flush();
  }
  // return CUDA device number based on how many appropriate ones exist on the current rank's node and what the relative rank number is
  return relative_ranknum % num_cuda_devices[rank];
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
        std::ostringstream out;
        out << "<- Rank " << OHMMS::Controller->rank() << " will use CUDA device #" << device ;
        out << " (" << deviceProp.name << ")" << std::endl;
        std::cerr << out.str();
        std::cerr.flush();
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
  if(gpu::rank==0)
    std::cerr << "Default MAX_GPU_SPLINE_SIZE_MB is " << gpu::MaxGPUSpineSizeMB << " MB." << std::endl;
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
