//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include <cstddef>
#include <cuda_runtime_api.h>
#include "cudaError.h"

size_t getCUDAdeviceFreeMem()
{
  size_t free, total;
  cudaErrorCheck(cudaMemGetInfo(&free, &total), "cudaMemGetInfo failed!");
  return free;
}
