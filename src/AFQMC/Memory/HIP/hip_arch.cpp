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

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <hip/hip_runtime.h>
#include "AFQMC/Memory/HIP/hip_init.h"
#include "AFQMC/Memory/HIP/hip_arch.h"
#include "AFQMC/Memory/device_pointers.hpp"
#include "AFQMC/Memory/buffer_allocators.h"
#include "mpi3/communicator.hpp"
#include "mpi3/shared_communicator.hpp"
#include "hipblas.h"
#include "hipsparse.h"
#include "rocsolver.h"
#include "rocrand/rocrand.h"

namespace qmcplusplus {
namespace afqmc {
  extern std::shared_ptr<device_allocator_generator_type> device_buffer_generator;
  extern std::shared_ptr<localTG_allocator_generator_type> localTG_buffer_generator;
}
}

namespace arch {

  hipblasHandle_t afqmc_hipblas_handle;
//  cublasXtHandle_t afqmc_cublasXt_handle;
  hipsparseHandle_t afqmc_hipsparse_handle;
  rocsolver_handle afqmc_rocsolver_handle;
  //hiprandGenerator_t afqmc_rand_generator;
  rocrand_generator afqmc_rocrand_generator;

  hipMemcpyKind tohipMemcpyKind(MEMCOPYKIND v) {
    switch(v)
    {
      case memcopyH2H:
      { return hipMemcpyHostToHost; }
      case memcopyH2D:
      { return hipMemcpyHostToDevice; }
      case memcopyD2H:
      { return hipMemcpyDeviceToHost;}
      case memcopyD2D:
      { return hipMemcpyDeviceToDevice; }
      case memcopyDefault:
      { return hipMemcpyDefault; }
    }
    return hipMemcpyDefault;
  }

  void INIT(boost::mpi3::shared_communicator& node, unsigned long long int iseed)
  {
    qmc_hip::HIP_INIT(node,iseed);
    using qmcplusplus::afqmc::device_buffer_generator;
    using qmcplusplus::afqmc::localTG_buffer_generator;
    using qmcplusplus::afqmc::device_allocator_generator_type;
    if(device_buffer_generator == nullptr ) {
      device_buffer_generator = std::make_shared<device_allocator_generator_type> (
                      device::memory_resource{},
                      std::size_t(20*1024*1024),
                      device::constructor<char>{} );
      // same memory space, so use same memory resource
      localTG_buffer_generator = device_buffer_generator;
    } else {
       std::cerr<<" Warning: device_buffer_generator already initialized in arch::INIT."
                <<std::endl;
    }
  }

  void memcopy(void* dst, const void* src, size_t count,
               MEMCOPYKIND kind) {
    if(hipSuccess != hipMemcpy(dst,src,count,tohipMemcpyKind(kind)))
      throw std::runtime_error("Error: hipMemcpy returned error code.");
  }

  void memcopy2D( void* dst, size_t dpitch, const void* src,
                  size_t spitch, size_t width, size_t height,
                  MEMCOPYKIND kind) {
    if(hipSuccess != hipMemcpy2D(dst,dpitch,src,spitch,width,height,tohipMemcpyKind(kind)))
      throw std::runtime_error("Error: hipMemcpy2D returned error code.");
  }

  void malloc(void** devPtr, size_t size) {
    if(hipSuccess != hipMalloc (devPtr,size)) {
      std::cerr<<" Error allocating " <<size*1024.0/1024.0 <<" MBs on GPU." <<std::endl;
      throw std::runtime_error("Error: hipMalloc returned error code.");
    }
  }

  void free(void *p) {
    hipFree(p);
  }

}

namespace device {

  boost::multi::array<std::complex<double>,1,device::device_allocator<std::complex<double>>>
                            *hipsparse_buffer(nullptr);

  arch::device_handles base_device_pointer::handles{&arch::afqmc_hipblas_handle,
//                                         &afqmc_cublasXt_handle,
                                         &arch::afqmc_hipsparse_handle,
                                         &arch::afqmc_rocsolver_handle,
                                         &arch::afqmc_rocrand_generator
                                        };

}

