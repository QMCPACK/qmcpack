//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////


#include <cassert>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <cuda_runtime.h>
#include "AFQMC/Memory/CUDA/cuda_init.h"
#include "cuda_arch.h"
#include "AFQMC/Memory/device_pointers.hpp"
#include "mpi3/communicator.hpp"
#include "mpi3/shared_communicator.hpp"
#include "cublas_v2.h"
//#include "cublasXt.h"
#include "cusparse.h"
#include "cusolverDn.h"
#include "curand.h"

namespace arch
{
cublasHandle_t afqmc_cublas_handle;
//  cublasXtHandle_t afqmc_cublasXt_handle;
cusparseHandle_t afqmc_cusparse_handle;
cusolverDnHandle_t afqmc_cusolverDn_handle;
curandGenerator_t afqmc_curand_generator;

cudaMemcpyKind tocudaMemcpyKind(MEMCOPYKIND v)
{
  switch (v)
  {
  case memcopyH2H: {
    return cudaMemcpyHostToHost;
  }
  case memcopyH2D: {
    return cudaMemcpyHostToDevice;
  }
  case memcopyD2H: {
    return cudaMemcpyDeviceToHost;
  }
  case memcopyD2D: {
    return cudaMemcpyDeviceToDevice;
  }
  case memcopyDefault: {
    return cudaMemcpyDefault;
  }
  }
  return cudaMemcpyDefault;
}

void INIT(boost::mpi3::shared_communicator& node, unsigned long long int iseed) { qmc_cuda::CUDA_INIT(node, iseed); }

void memcopy(void* dst, const void* src, size_t count, MEMCOPYKIND kind, const std::string& message)
{
  cudaError_t status = cudaMemcpy(dst, src, count, tocudaMemcpyKind(kind));
  if (status != cudaSuccess)
  {
    if (message != "")
    {
      std::cerr << "Error: " << message << std::endl;
    }
    std::cerr << " Error when calling cudaMemcpy: " << cudaGetErrorString(status) << std::endl;
    throw std::runtime_error("Error: cudaMemcpy returned error code.");
  }
}

void memcopy2D(void* dst,
               size_t dpitch,
               const void* src,
               size_t spitch,
               size_t width,
               size_t height,
               MEMCOPYKIND kind,
               const std::string& message)
{
  cudaError_t status = cudaMemcpy2D(dst, dpitch, src, spitch, width, height, tocudaMemcpyKind(kind));
  if (status != cudaSuccess)
  {
    if (message != "")
    {
      std::cerr << "Error: " << message << std::endl;
    }
    std::cerr << " Error when calling cudaMemcpy2D: " << cudaGetErrorString(status) << std::endl;
    throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
  }
}

void malloc(void** devPtr, size_t size, const std::string& message)
{
  cudaError_t status = cudaMalloc(devPtr, size);
  if (status != cudaSuccess)
  {
    std::cerr << " Error allocating " << size / 1024.0 / 1024.0 << " MBs on GPU." << std::endl;
    if (message != "")
    {
      std::cerr << " Error from: " << message << std::endl;
    }
    std::cerr << " Error when calling cudaMalloc: " << cudaGetErrorString(status) << std::endl;
    throw std::runtime_error("Error: cudaMalloc returned error code.");
  }
}

void free(void* p, const std::string& message)
{
  cudaError_t status = cudaFree(p);
  if (status != cudaSuccess)
  {
    if (message != "")
    {
      std::cerr << " Error from: " << message << std::endl;
    }
    std::cerr << " Error from calling cudaFree: " << cudaGetErrorString(status) << std::endl;
  }
}

} // namespace arch

namespace device
{
boost::multi::array<std::complex<double>, 1, device::device_allocator<std::complex<double>>>* cusparse_buffer(nullptr);

arch::device_handles base_device_pointer::handles{&arch::afqmc_cublas_handle,
                                                  //                                         &afqmc_cublasXt_handle,
                                                  &arch::afqmc_cusparse_handle, &arch::afqmc_cusolverDn_handle,
                                                  &arch::afqmc_curand_generator};

} // namespace device
