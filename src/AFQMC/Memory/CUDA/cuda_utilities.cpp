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

//#ifndef ENABLE_CUDA
//#error
//#endif

#include<cassert>
#include<complex>
#include<cstdlib>
#include<stdexcept>
#include "Utilities/OutputManager.h"
#include "AFQMC/Memory/CUDA/cuda_utilities.h"
#include "AFQMC/Memory/CUDA/cuda_gpu_pointer.hpp"
#include <cuda_runtime.h>
#include "cublas_v2.h"
//#include "cublasXt.h"
#include "cusparse.h"
#include "cusolverDn.h"
#include "curand.h"
#include "mpi3/communicator.hpp"
#include "mpi3/shared_communicator.hpp"

#include "multi/array.hpp"
#include "multi/array_ref.hpp"

namespace qmc_cuda {

  // work around for problem with csrmm 
  boost::multi::array<std::complex<double>,1,qmc_cuda::cuda_gpu_allocator<std::complex<double>>> 
                            *cusparse_buffer(nullptr);
                                        //(typename boost::multi::layout_t<1u>::extensions_type{1},
                                        //qmc_cuda::cuda_gpu_allocator<std::complex<double>>{});

  cublasHandle_t afqmc_cublas_handle;
//  cublasXtHandle_t afqmc_cublasXt_handle;
  cusparseHandle_t afqmc_cusparse_handle;
  cusolverDnHandle_t afqmc_cusolverDn_handle;
  curandGenerator_t afqmc_curand_generator;
  bool afqmc_cuda_handles_init = false;
  cusparseMatDescr_t afqmc_cusparse_matrix_descr;

  std::vector<cudaStream_t> afqmc_cuda_streams;

  gpu_handles base_cuda_gpu_ptr::handles{&afqmc_cublas_handle,
//                                         &afqmc_cublasXt_handle,
                                         &afqmc_cusparse_handle,
                                         &afqmc_cusolverDn_handle,
                                         &afqmc_curand_generator  
                                        }; 
  // need a cleanup routine
  void CUDA_INIT(boost::mpi3::shared_communicator& node, unsigned long long int iseed)
  {

    if(afqmc_cuda_handles_init) return;
    afqmc_cuda_handles_init=true;

    int num_devices=0;
    cudaGetDeviceCount(&num_devices);
    qmcplusplus::app_log()<<" Running in node with " <<num_devices <<" GPUs. \n";
    if(num_devices < node.size()) {
      qmcplusplus::app_error()<<"Error: # GPU < # tasks in node. " <<std::endl;
      qmcplusplus::app_error()<<"# GPU: " <<num_devices <<std::endl;
      qmcplusplus::app_error()<<"# tasks: " <<node.size() <<std::endl;
      APP_ABORT("");
    } else if(num_devices > node.size()) {
      qmcplusplus::app_log()<<"WARNING: Unused devices !!!!!!!!!!!!!! \n"
                                <<"         # tasks: " <<node.size() <<"\n"
                                <<"         num_devices: " <<num_devices <<std::endl;
    }

    cuda_check(cudaSetDevice(node.rank()),"cudaSetDevice()");

    cublas_check(cublasCreate (& afqmc_cublas_handle ), "cublasCreate");
//    cublas_check(cublasXtCreate (& afqmc_cublasXt_handle ), "cublasXtCreate");
    int devID[8] {0,1,2,3,4,5,6,7};
//    cublas_check(cublasXtDeviceSelect(afqmc_cublasXt_handle, 1, devID), "cublasXtDeviceSelect");
//    cublas_check(cublasXtSetPinningMemMode(afqmc_cublasXt_handle, CUBLASXT_PINNING_ENABLED), 
//                                            "cublasXtSetPinningMemMode");
    cusolver_check(cusolverDnCreate (& afqmc_cusolverDn_handle ), "cusolverDnCreate");
    //curand_check(curandCreateGenerator(&afqmc_curand_generator, CURAND_RNG_PSEUDO_DEFAULT),
    curand_check(curandCreateGenerator(&afqmc_curand_generator, CURAND_RNG_PSEUDO_MT19937),
                                            "curandCreateGenerator");
    curand_check(curandSetPseudoRandomGeneratorSeed(afqmc_curand_generator,iseed),
                                            "curandSetPseudoRandomGeneratorSeed");

    cusparse_check(cusparseCreate (& afqmc_cusparse_handle ), "cusparseCreate");
    cusparse_check(cusparseCreateMatDescr(&afqmc_cusparse_matrix_descr), 
            "cusparseCreateMatDescr: Matrix descriptor initialization failed"); 
    cusparseSetMatType(afqmc_cusparse_matrix_descr,CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(afqmc_cusparse_matrix_descr,CUSPARSE_INDEX_BASE_ZERO); 

    cusparse_buffer = new boost::multi::array<std::complex<double>,1,
                                 qmc_cuda::cuda_gpu_allocator<std::complex<double>>>(
                                 (typename boost::multi::layout_t<1u>::extensions_type{1},
                                 qmc_cuda::cuda_gpu_allocator<std::complex<double>>{}));

  }

  void cuda_check(cudaError_t sucess, std::string message)
  {
    if(cudaSuccess != sucess) {
      std::cerr<<message <<std::endl;
      std::cerr<<" cudaGetErrorName: " <<cudaGetErrorName(sucess) <<std::endl;
      std::cerr<<" cudaGetErrorString: " <<cudaGetErrorString(sucess) <<std::endl;
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

  void cusparse_check(cusparseStatus_t sucess, std::string message)
  {
    if(CUSPARSE_STATUS_SUCCESS != sucess) {
      std::cerr<<message <<std::endl;
      std::cerr.flush();
      throw std::runtime_error(" Error code returned by cusparse. \n");
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
    }
  }

  cusparseOperation_t cusparseOperation(char A) {
    if(A=='N' or A=='n')
      return CUSPARSE_OPERATION_NON_TRANSPOSE;
    else if(A=='T' or A=='t')
      return CUSPARSE_OPERATION_TRANSPOSE;
    else if(A=='C' or A=='c')
      return CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE;
    else {
      throw std::runtime_error("unknown cusparseOperation option"); 
    }
  }

}

