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

#ifndef AFQMC_CUDA_ARCH_HPP
#define AFQMC_CUDA_ARCH_HPP

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <cuda_runtime.h>
#include "AFQMC/Memory/CUDA/cuda_init.h"
#include "mpi3/communicator.hpp"
#include "mpi3/shared_communicator.hpp"
#include "cublas_v2.h"
//#include "cublasXt.h"
#include "cusparse.h"
#include "cusolverDn.h"
#include "curand.h"


namespace arch
{
extern curandGenerator_t afqmc_curand_generator;

struct device_handles
{
  cublasHandle_t* cublas_handle;
  //    cublasXtHandle_t* cublasXt_handle;
  cusparseHandle_t* cusparse_handle;
  cusolverDnHandle_t* cusolverDn_handle;
  curandGenerator_t* curand_generator;
  bool operator==(device_handles const& other) const
  {
    return (cublas_handle == other.cublas_handle &&
            //              cublasXt_handle==other.cublasXt_handle &&
            cusparse_handle == other.cusparse_handle && cusolverDn_handle == other.cusolverDn_handle &&
            curand_generator == other.curand_generator);
  }
  bool operator!=(device_handles const& other) const { return not(*this == other); }
};


enum MEMCOPYKIND
{
  memcopyH2H     = cudaMemcpyHostToHost,
  memcopyH2D     = cudaMemcpyHostToDevice,
  memcopyD2H     = cudaMemcpyDeviceToHost,
  memcopyD2D     = cudaMemcpyDeviceToDevice,
  memcopyDefault = cudaMemcpyDefault
};

cudaMemcpyKind tocudaMemcpyKind(MEMCOPYKIND v);

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
