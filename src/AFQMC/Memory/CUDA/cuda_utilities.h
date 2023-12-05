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

#ifndef AFQMC_CUDA_UTILITIES_HPP
#define AFQMC_CUDA_UTILITIES_HPP

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <cuda_runtime.h>
#include "cublas_v2.h"
//#include "cublasXt.h"
#include "cusparse.h"
#include "cusolverDn.h"
#include "curand.h"

#ifdef BUILD_AFQMC_WITH_NCCL
#include "nccl.h"

#define NCCLCHECK(cmd)                                                                      \
  do                                                                                        \
  {                                                                                         \
    ncclResult_t r = cmd;                                                                   \
    if (r != ncclSuccess)                                                                   \
    {                                                                                       \
      printf("Failed, NCCL error %s:%d '%s'\n", __FILE__, __LINE__, ncclGetErrorString(r)); \
      exit(EXIT_FAILURE);                                                                   \
    }                                                                                       \
  } while (0)

#endif

namespace qmc_cuda
{
//  extern curandGenerator_t afqmc_curand_generator;
extern cusparseMatDescr_t afqmc_cusparse_matrix_descr;

extern std::vector<cudaStream_t> afqmc_cuda_streams;

void cuda_check_error();
void cuda_check(cudaError_t success, std::string message = "");
void cublas_check(cublasStatus_t success, std::string message = "");
void cusparse_check(cusparseStatus_t success, std::string message = "");
void curand_check(curandStatus_t success, std::string message = "");
void cusolver_check(cusolverStatus_t success, std::string message = "");
cublasOperation_t cublasOperation(char A);
cusparseOperation_t cusparseOperation(char A);

// since when do these exist? Not sure, so keep guard for now.
#if CUSPARSE_VER_MAJOR > 10
template<typename T>
cudaDataType_t cusparse_data_type()
{
  return cudaDataType_t{};
}
template<>
inline cudaDataType_t cusparse_data_type<float>()
{
  return CUDA_R_32F;
}
template<>
inline cudaDataType_t cusparse_data_type<double>()
{
  return CUDA_R_64F;
}
template<>
inline cudaDataType_t cusparse_data_type<std::complex<float>>()
{
  return CUDA_C_32F;
}
template<>
inline cudaDataType_t cusparse_data_type<std::complex<double>>()
{
  return CUDA_C_64F;
}
#endif


} // namespace qmc_cuda

#endif
