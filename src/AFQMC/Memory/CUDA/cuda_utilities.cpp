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

#include <cassert>
#include <complex>
#include <cstdlib>
#include <stdexcept>
#include "Platforms/Host/OutputManager.h"
#include "cuda_utilities.h"
#include "AFQMC/Memory/device_pointers.hpp"
#include <cuda_runtime.h>
#include "cublas_v2.h"
//#include "cublasXt.h"
#include "cusparse.h"
#include "cusolverDn.h"
#include "curand.h"

#include "multi/array.hpp"
#include "multi/array_ref.hpp"


namespace qmc_cuda
{
bool afqmc_cuda_handles_init = false;
cusparseMatDescr_t afqmc_cusparse_matrix_descr;

std::vector<cudaStream_t> afqmc_cuda_streams;

void cuda_check(cudaError_t success, std::string message)
{
  if (cudaSuccess != success)
  {
    std::cerr << message << std::endl;
    std::cerr << " cudaGetErrorName: " << cudaGetErrorName(success) << std::endl;
    std::cerr << " cudaGetErrorString: " << cudaGetErrorString(success) << std::endl;
    std::cerr.flush();
    throw std::runtime_error(" Error code returned by cuda. \n");
  }
}

void cublas_check(cublasStatus_t success, std::string message)
{
  if (CUBLAS_STATUS_SUCCESS != success)
  {
    std::cerr << message << std::endl;
    std::cerr.flush();
    throw std::runtime_error(" Error code returned by cublas. \n");
  }
}

void cusparse_check(cusparseStatus_t success, std::string message)
{
  if (CUSPARSE_STATUS_SUCCESS != success)
  {
    std::cerr << message << std::endl;
    std::cerr.flush();
    throw std::runtime_error(" Error code returned by cusparse. \n");
  }
}

void curand_check(curandStatus_t success, std::string message)
{
  if (CURAND_STATUS_SUCCESS != success)
  {
    std::cerr << message << std::endl;
    std::cerr.flush();
    throw std::runtime_error(" Error code returned by curand. \n");
  }
}

void cusolver_check(cusolverStatus_t success, std::string message)
{
  if (CUSOLVER_STATUS_SUCCESS != success)
  {
    std::cerr << message << std::endl;
    std::cerr.flush();
    throw std::runtime_error(" Error code returned by cusolver. \n");
  }
}

cublasOperation_t cublasOperation(char A)
{
  if (A == 'N' or A == 'n')
    return CUBLAS_OP_N;
  else if (A == 'T' or A == 't')
    return CUBLAS_OP_T;
  else if (A == 'C' or A == 'c')
    return CUBLAS_OP_C;
  else
  {
    throw std::runtime_error("unknown cublasOperation option");
  }
}

cusparseOperation_t cusparseOperation(char A)
{
  if (A == 'N' or A == 'n')
    return CUSPARSE_OPERATION_NON_TRANSPOSE;
  else if (A == 'T' or A == 't')
    return CUSPARSE_OPERATION_TRANSPOSE;
  else if (A == 'C' or A == 'c')
    return CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE;
  else
  {
    throw std::runtime_error("unknown cusparseOperation option");
  }
}

} // namespace qmc_cuda
