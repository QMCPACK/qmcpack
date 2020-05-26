//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_MEMORY_UTILITIES_HPP
#define AFQMC_MEMORY_UTILITIES_HPP

#ifdef ENABLE_CUDA
#include <cuda_runtime.h>
#include "AFQMC/Memory/CUDA/cuda_utilities.h"
#elif ENABLE_HIP
#include <hip_runtime.h>
#include "AFQMC/Memory/HIP/cuda_utilities.h"
#endif
#include "Message/OpenMP.h"

namespace qmc_cuda {
extern bool afqmc_cuda_handles_init;
}

inline int number_of_devices()
{
  int num_devices=0;
#ifdef ENABLE_CUDA
  if(not qmc_cuda::afqmc_cuda_handles_init)
    throw std::runtime_error(" Error: Uninitialized CUDA environment.");
  cudaGetDeviceCount(&num_devices);
#elif ENABLE_HIP
  if(not qmc_hip::afqmc_hip_handles_init)
    throw std::runtime_error(" Error: Uninitialized HIP environment.");
  hipGetDeviceCount(&num_devices);
#endif
  return num_devices;
}

#endif
