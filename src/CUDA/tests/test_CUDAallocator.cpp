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


#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <memory>
#include <iostream>
#include "CUDA/cudaError.h"
#include "CUDA/CUDAallocator.hpp"
#include "OhmmsPETE/OhmmsVector.h"

namespace qmcplusplus
{
TEST_CASE("CUDA_allocators", "[CUDA]")
{
  { // CUDAManagedAllocator
    Vector<double, CUDAManagedAllocator<double>> vec(1024);
    cudaPointerAttributes attr;
    cudaErrorCheck(cudaPointerGetAttributes(&attr, vec.data()) , "cudaPointerGetAttributes failed!");
    REQUIRE(attr.type == cudaMemoryTypeManaged);
  }
  { // CUDAAllocator
    Vector<double, CUDAAllocator<double>> vec(1024);
    cudaPointerAttributes attr;
    cudaErrorCheck(cudaPointerGetAttributes(&attr, vec.data()) , "cudaPointerGetAttributes failed!");
    REQUIRE(attr.type == cudaMemoryTypeDevice);
  }
  { // CUDAHostAllocator
    Vector<double, CUDAHostAllocator<double>> vec(1024);
    cudaPointerAttributes attr;
    cudaErrorCheck(cudaPointerGetAttributes(&attr, vec.data()) , "cudaPointerGetAttributes failed!");
    REQUIRE(attr.type == cudaMemoryTypeHost);
  }
  { // CUDALockedPageAllocator
    Vector<double, CUDALockedPageAllocator<double>> vec(1024);
    cudaPointerAttributes attr;
    cudaErrorCheck(cudaPointerGetAttributes(&attr, vec.data()) , "cudaPointerGetAttributes failed!");
    REQUIRE(attr.type == cudaMemoryTypeHost);
  }
}

} // namespace qmcplusplus
