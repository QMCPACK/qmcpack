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

#ifndef AFQMC_HIP_ARCH_HPP
#define AFQMC_HIP_ARCH_HPP

#include<cassert>
#include<cstdlib>
#include<iostream>
#include<stdexcept>
#include <hip/hip_runtime.h>
#include "AFQMC/Memory/HIP/hip_init.h"
#include "mpi3/communicator.hpp"
#include "mpi3/shared_communicator.hpp"
#include "hipblas.h"
//#include "cublasXt.h"
#include "hipsparse.h"
#include "cusolverDn.h"
#ifdef ENABLE_ROCM
#include "rocsolver.h"
#endif
#include "hiprand.hpp"


namespace arch {

  extern hiprandGenerator_t afqmc_curand_generator;

  struct device_handles {
    hipblasHandle_t* blas_handle;
//    cublasXtHandle_t* cublasXt_handle;
    hipsparseHandle_t* sparse_handle;
    hiprandGenerator_t* rand_generator;
#ifdef ENABLE_ROCM
    rocsolver_handle solver_handle;
#endif
    bool operator==(device_handles const& other) const {
      return (blas_handle==other.blas_handle &&
//              cublasXt_handle==other.cublasXt_handle &&
              sparse_handle==other.sparse_handle &&
              solver_handle==other.solver_handle &&
              rand_generator==other.rand_generator);
    }
    bool operator!=(device_handles const& other) const { return not (*this == other); }
  };


  enum MEMCOPYKIND { memcopyH2H = hipMemcpyHostToHost, 
                     memcopyH2D = hipMemcpyHostToDevice, 
                     memcopyD2H = hipMemcpyDeviceToHost, 
                     memcopyD2D = hipMemcpyDeviceToDevice,
                     memcopyDefault = hipMemcpyDefault };

  hipMemcpyKind tohipMemcpyKind(MEMCOPYKIND v);

  void INIT(boost::mpi3::shared_communicator& node, unsigned long long int iseed = 911ULL);

  void memcopy(void* dst, const void* src, size_t count, 
               MEMCOPYKIND kind = memcopyDefault );

  void memcopy2D( void* dst, size_t dpitch, const void* src, 
                  size_t spitch, size_t width, size_t height, 
                  MEMCOPYKIND kind = memcopyDefault);

  void malloc(void** devPtr, size_t size);

  void free(void *p);

}

#endif

