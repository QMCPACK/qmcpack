///////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Fionn Malone, malone14@llnl.gov, LLNL
//
// File created by: Fionn Malone, malone14@llnl.gov, LLNL
////////////////////////////////////////////////////////////////////////////////


#include "hip_arch.h"
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <hip/hip_runtime.h>
#include "AFQMC/Memory/device_pointers.hpp"
#include "mpi3/communicator.hpp"
#include "mpi3/shared_communicator.hpp"
#include <hipblas/hipblas.h>
#include <hipsparse/hipsparse.h>
#include <rocsolver/rocsolver.h>
#include <rocrand/rocrand.h>

namespace arch
{
hipblasHandle_t afqmc_hipblas_handle;
//  cublasXtHandle_t afqmc_cublasXt_handle;
hipsparseHandle_t afqmc_hipsparse_handle;
rocsolver_handle afqmc_rocsolver_handle;
//hiprandGenerator_t afqmc_rand_generator;
rocrand_generator afqmc_rocrand_generator;

hipMemcpyKind tohipMemcpyKind(MEMCOPYKIND v)
{
  switch (v)
  {
  case memcopyH2H: {
    return hipMemcpyHostToHost;
  }
  case memcopyH2D: {
    return hipMemcpyHostToDevice;
  }
  case memcopyD2H: {
    return hipMemcpyDeviceToHost;
  }
  case memcopyD2D: {
    return hipMemcpyDeviceToDevice;
  }
  case memcopyDefault: {
    return hipMemcpyDefault;
  }
  }
  return hipMemcpyDefault;
}

void INIT(boost::mpi3::shared_communicator& node, unsigned long long int iseed) { qmc_hip::HIP_INIT(node, iseed); }

void memcopy(void* dst, const void* src, size_t count, MEMCOPYKIND kind, const std::string& message)
{
  hipError_t status = hipMemcpy(dst, src, count, tohipMemcpyKind(kind));
  if (status != hipSuccess)
  {
    if (message != "")
    {
      std::cerr << "Error: " << message << std::endl;
    }
    std::cerr << " Error when calling hipMemcpy: " << hipGetErrorString(status) << std::endl;
    throw std::runtime_error("Error: hipMemcpy returned error code.");
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
  hipError_t status = hipMemcpy2D(dst, dpitch, src, spitch, width, height, tohipMemcpyKind(kind));
  if (status != hipSuccess)
  {
    if (message != "")
    {
      std::cerr << "Error: " << message << std::endl;
    }
    std::cerr << " Error when calling hipMemcpy2D: " << hipGetErrorString(status) << std::endl;
    throw std::runtime_error("Error: hipMemcpy2D returned error code.");
  }
}

void malloc(void** devPtr, size_t size, const std::string& message)
{
  hipError_t status = hipMalloc(devPtr, size);
  if (status != hipSuccess)
  {
    std::cerr << " Error allocating " << size * 1024.0 / 1024.0 << " MBs on GPU." << std::endl;
    if (message != "")
    {
      std::cerr << " Error from : " << message << std::endl;
    }
    std::cerr << " Error when call hipMalloc " << hipGetErrorString(status) << std::endl;
    throw std::runtime_error("Error: hipMalloc returned error code.");
  }
}

void free(void* p, const std::string& message)
{
  hipError_t status = hipFree(p);
  if (status != hipSuccess)
  {
    if (message != "")
    {
      std::cerr << " Error from : " << message << std::endl;
    }
    std::cerr << " Error when calling hipFree: " << hipGetErrorString(status) << std::endl;
  }
}

} // namespace arch

namespace device
{
boost::multi::array<std::complex<double>, 1, device::device_allocator<std::complex<double>>>* hipsparse_buffer(nullptr);

arch::device_handles base_device_pointer::handles{&arch::afqmc_hipblas_handle,
                                                  //                                         &afqmc_cublasXt_handle,
                                                  &arch::afqmc_hipsparse_handle, &arch::afqmc_rocsolver_handle,
                                                  &arch::afqmc_rocrand_generator};

} // namespace device
