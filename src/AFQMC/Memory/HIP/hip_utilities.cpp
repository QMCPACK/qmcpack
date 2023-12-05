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
#include <complex>
#include <cstdlib>
#include <stdexcept>
#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
#include <hipsparse/hipsparse.h>
#include <rocsolver/rocsolver.h>
#include <rocrand/rocrand.h>
#include "hip_utilities.h"

#include "multi/array.hpp"
#include "multi/array_ref.hpp"


namespace qmc_hip
{
bool afqmc_hip_handles_init = false;
hipsparseMatDescr_t afqmc_hipsparse_matrix_descr;

std::vector<hipStream_t> afqmc_hip_streams;

void hip_check(hipError_t success, std::string message)
{
  if (hipSuccess != success)
  {
    std::cerr << message << std::endl;
    std::cerr << " hipGetErrorName: " << hipGetErrorName(success) << std::endl;
    std::cerr << " hipGetErrorString: " << hipGetErrorString(success) << std::endl;
    std::cerr.flush();
    throw std::runtime_error(" Error code returned by hip. \n");
  }
}

void hipblas_check(hipblasStatus_t success, std::string message)
{
  if (HIPBLAS_STATUS_SUCCESS != success)
  {
    std::cerr << message << std::endl;
    std::cerr.flush();
    throw std::runtime_error(" Error code returned by cublas. \n");
  }
}

void hipsparse_check(hipsparseStatus_t success, std::string message)
{
  if (HIPSPARSE_STATUS_SUCCESS != success)
  {
    std::cerr << message << std::endl;
    std::cerr.flush();
    throw std::runtime_error(" Error code returned by cusparse. \n");
  }
}

void hiprand_check(hiprandStatus_t success, std::string message)
{
  if (ROCRAND_STATUS_SUCCESS != success)
  {
    std::cerr << message << std::endl;
    std::cerr.flush();
    throw std::runtime_error(" Error code returned by hiprand. \n");
  }
}

void hipsolver_check(hipsolverStatus_t success, std::string message)
{
  if (rocblas_status_success != success)
  {
    std::cerr << message << std::endl;
    std::cerr.flush();
    throw std::runtime_error(" Error code returned by hipsolver. \n");
  }
}

hipblasOperation_t hipblasOperation(char A)
{
  if (A == 'N' or A == 'n')
    return HIPBLAS_OP_N;
  else if (A == 'T' or A == 't')
    return HIPBLAS_OP_T;
  else if (A == 'C' or A == 'c')
    return HIPBLAS_OP_C;
  else
  {
    throw std::runtime_error("unknown hipblasOperation option");
  }
}

rocblasOperation_t rocblasOperation(char A)
{
  if (A == 'N' or A == 'n')
    return rocblas_operation_none;
  else if (A == 'T' or A == 't')
    return rocblas_operation_transpose;
  else if (A == 'C' or A == 'c')
    return rocblas_operation_conjugate_transpose;
  else
  {
    throw std::runtime_error("unknown hipblasOperation option");
  }
}

hipsparseOperation_t hipsparseOperation(char A)
{
  if (A == 'N' or A == 'n')
    return HIPSPARSE_OPERATION_NON_TRANSPOSE;
  else if (A == 'T' or A == 't')
    return HIPSPARSE_OPERATION_TRANSPOSE;
  else if (A == 'C' or A == 'c')
    return HIPSPARSE_OPERATION_CONJUGATE_TRANSPOSE;
  else
  {
    throw std::runtime_error("unknown hipsparseOperation option");
  }
}

} // namespace qmc_hip
