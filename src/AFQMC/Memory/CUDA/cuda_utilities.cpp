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

//#ifndef QMC_CUDA
//#error
//#endif

#include<cassert>
#include<cstdlib>
#include<stdexcept>
#include "Utilities/OutputManager.h"
#include "AFQMC/Memory/CUDA/cuda_utilities.h"
#include "AFQMC/Memory/CUDA/cuda_gpu_pointer.hpp"
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include "cublasXt.h"
#include "cusolverDn.h"
#include "curand.h"
#include "mpi3/communicator.hpp"
#include "mpi3/shared_communicator.hpp"

namespace qmc_cuda {

  cublasHandle_t afqmc_cublas_handle;
  cublasXtHandle_t afqmc_cublasXt_handle;
  cusolverDnHandle_t afqmc_cusolverDn_handle;
  curandGenerator_t afqmc_curand_generator;

  gpu_handles base_cuda_gpu_ptr::handles{&afqmc_cublas_handle,
                                         &afqmc_cublasXt_handle,
                                         &afqmc_cusolverDn_handle,
                                         &afqmc_curand_generator  
                                        };

  void CUDA_INIT(boost::mpi3::shared_communicator& node)
  {

    int num_devices=0;
    cudaGetDeviceCount(&num_devices);
    if(num_devices < node.size()) {
      qmcplusplus::app_error()<<"Error: # GPU < # tasks in node. " <<std::endl;
      qmcplusplus::app_error()<<"# GPU: " <<num_devices <<std::endl;
      qmcplusplus::app_error()<<"# tasks: " <<node.size() <<std::endl;
      APP_ABORT("");
    }

    qmcplusplus::app_log()<<" Running in node with " <<num_devices <<" GPUs. \n";

    cuda_check(cudaSetDevice(node.rank()),"cudaSetDevice()");

    cublas_check(cublasCreate (& afqmc_cublas_handle ), "cublasCreate");
    cublas_check(cublasXtCreate (& afqmc_cublasXt_handle ), "cublasXtCreate");
    int devID[8] {0,1,2,3,4,5,6,7};
    cublas_check(cublasXtDeviceSelect(afqmc_cublasXt_handle, 1, devID), "cublasXtDeviceSelect");
    cublas_check(cublasXtSetPinningMemMode(afqmc_cublasXt_handle, CUBLASXT_PINNING_ENABLED), 
                                            "cublasXtSetPinningMemMode");
    cusolver_check(cusolverDnCreate (& afqmc_cusolverDn_handle ), "cusolverDnCreate");
    curand_check(curandCreateGenerator(&afqmc_curand_generator, CURAND_RNG_PSEUDO_DEFAULT),
                                            "curandCreateGenerator");
    curand_check(curandSetPseudoRandomGeneratorSeed(afqmc_curand_generator,1234ULL),
                                            "curandSetPseudoRandomGeneratorSeed");

  }

  void cuda_check(cudaError_t sucess, std::string message)
  {
    if(cudaSuccess != sucess) {
      std::cerr<<message <<std::endl;
      std::cerr.flush();
      throw std::runtime_error(" Error code returned by cuda. \n");
    }
  }

  void cublas_check(cublasStatus_t sucess, std::string message)
  { 
    if(CUBLAS_STATUS_SUCCESS != sucess) {
      std::cerr<<message <<std::endl;
      std::cerr.flush();
      throw std::runtime_error(" Error code returned by cublas. \n");
    }
  }

  void curand_check(curandStatus_t sucess, std::string message)
  {
    if(CURAND_STATUS_SUCCESS != sucess) {
      std::cerr<<message <<std::endl;
      std::cerr.flush();
      throw std::runtime_error(" Error code returned by curand. \n");
    }
  }

  void cusolver_check(cusolverStatus_t sucess, std::string message)
  {
    if(CUSOLVER_STATUS_SUCCESS != sucess) {
      std::cerr<<message <<std::endl;
      std::cerr.flush();
      throw std::runtime_error(" Error code returned by cusolver. \n");
    }
  }

  cublasOperation_t cublasOperation(char A) {
    if(A=='N' or A=='n')
      return CUBLAS_OP_N;
    else if(A=='T' or A=='t')
      return CUBLAS_OP_T;
    else if(A=='C' or A=='c')
      return CUBLAS_OP_C;
    else {
      throw std::runtime_error("unknown cublasOperation option"); 
      //return CUBLAS_OP_N;
    }
  }

}

