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
#include "AFQMC/Memory/CUDA/cuda_arch.h"
#include "AFQMC/Memory/device_pointers.hpp"
#include "AFQMC/Memory/buffer_allocators.h"
#include "mpi3/communicator.hpp"
#include "mpi3/shared_communicator.hpp"
#include "cublas_v2.h"
//#include "cublasXt.h"
#include "cusparse.h"
#include "cusolverDn.h"
#include "curand.h"

namespace qmcplusplus
{
namespace afqmc
{
extern std::shared_ptr<device_allocator_generator_type> device_buffer_generator;
extern std::shared_ptr<localTG_allocator_generator_type> localTG_buffer_generator;
} // namespace afqmc
} // namespace qmcplusplus

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

void INIT(boost::mpi3::shared_communicator& node, unsigned long long int iseed)
{
  qmc_cuda::CUDA_INIT(node, iseed);
  using qmcplusplus::afqmc::device_allocator_generator_type;
  using qmcplusplus::afqmc::device_buffer_generator;
  using qmcplusplus::afqmc::localTG_buffer_generator;
  if (device_buffer_generator == nullptr)
  {
    device_buffer_generator =
        std::make_shared<device_allocator_generator_type>(device::memory_resource{}, std::size_t(20 * 1024 * 1024),
                                                          device::constructor<char>{});
    // same memory space, so use same memory resource
    localTG_buffer_generator = device_buffer_generator;
  }
  else
  {
    std::cerr << " Warning: device_buffer_generator already initialized in arch::INIT." << std::endl;
  }
}

void memcopy(void* dst, const void* src, size_t count, MEMCOPYKIND kind, std::string message)
{
  if (cudaSuccess != cudaMemcpy(dst, src, count, tocudaMemcpyKind(kind)))
    throw std::runtime_error("Error: cudaMemcpy returned error code.");
}

void memcopy2D(void* dst, size_t dpitch, const void* src, size_t spitch, size_t width, size_t height, MEMCOPYKIND kind, std::string message)
{
  if (cudaSuccess != cudaMemcpy2D(dst, dpitch, src, spitch, width, height, tocudaMemcpyKind(kind)))
    throw std::runtime_error("Error: cudaMemcpy2D returned error code.");
}

void malloc(void** devPtr, size_t size, std::string message)
{
  if (cudaSuccess != cudaMalloc(devPtr, size))
  {
    std::cerr << " Error allocating " << size * 1024.0 / 1024.0 << " MBs on GPU." << std::endl;
    if (message != "")
    {
      std::cerr << " Error from : " << message << std::endl;
    }
    throw std::runtime_error("Error: cudaMalloc returned error code.");
  }
}

void free(void* p, std::string message)
{
  if (cudaSuccess != cudaFree(p))
  {
    if (message != "")
    {
      std::cerr << " Error from: " << message << std::endl;
    }
    throw std::runtime_error("Error: cudaFree returned error code.");
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
