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

#ifndef AFQMC_HIP_ARCH_HPP
#define AFQMC_HIP_ARCH_HPP

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <hip/hip_runtime.h>
#include "hip_init.h"
#include "mpi3/communicator.hpp"
#include "mpi3/shared_communicator.hpp"
#include <hipblas/hipblas.h>
#include <hipsparse/hipsparse.h>
#include <rocsolver/rocsolver.h>
#include <rocrand/rocrand.h>


namespace arch
{
extern rocrand_generator afqmc_rocrand_generator;

struct device_handles
{
  hipblasHandle_t* hipblas_handle;
  //    cublasXtHandle_t* cublasXt_handle;
  hipsparseHandle_t* hipsparse_handle;
  //hiprandGenerator_t* rand_generator;
  rocsolver_handle* rocsolver_handle_;
  rocrand_generator* rocrand_generator_;
  bool operator==(device_handles const& other) const
  {
    return (hipblas_handle == other.hipblas_handle &&
            //              cublasXt_handle==other.cublasXt_handle &&
            hipsparse_handle == other.hipsparse_handle && rocsolver_handle_ == other.rocsolver_handle_ &&
            rocrand_generator_ == other.rocrand_generator_);
  }
  bool operator!=(device_handles const& other) const { return not(*this == other); }
};


enum MEMCOPYKIND
{
  memcopyH2H     = hipMemcpyHostToHost,
  memcopyH2D     = hipMemcpyHostToDevice,
  memcopyD2H     = hipMemcpyDeviceToHost,
  memcopyD2D     = hipMemcpyDeviceToDevice,
  memcopyDefault = hipMemcpyDefault
};

hipMemcpyKind tohipMemcpyKind(MEMCOPYKIND v);

void INIT(boost::mpi3::shared_communicator& node, unsigned long long int iseed = 911ULL);

void memcopy(void* dst,
             const void* src,
             size_t count,
             MEMCOPYKIND kind           = memcopyDefault,
             const std::string& message = "");

void memcopy2D(void* dst,
               size_t dpitch,
               const void* src,
               size_t spitch,
               size_t width,
               size_t height,
               MEMCOPYKIND kind           = memcopyDefault,
               const std::string& message = "");

void malloc(void** devPtr, size_t size, const std::string& message = "");

void free(void* p, const std::string& message = "");

} // namespace arch

#endif
