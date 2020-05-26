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

#ifndef AFQMC_HIP_UTILITIES_HPP
#define AFQMC_HIP_UTILITIES_HPP

#include<cassert>
#include<cstdlib>
#include<iostream>
#include<stdexcept>
#include<vector>
#include <hip/hip_runtime.h>
#include "hipblas.h"
//#include "cublasXt.h"
#include "hipsparse.h"
#ifdef ENABLE_ROCM
#include "rocsolver.h"
#endif
#include "hiprand.h"

namespace qmc_hip {

//  extern hiprandGenerator_t afqmc_curand_generator;
  extern hipsparseMatDescr_t afqmc_sparse_matrix_descr;

  extern std::vector<hipStream_t> afqmc_hip_streams;

  // FDM: Temprorary hack to allow for easier grepping.
  typedef rocsolver_status = rocsolverStatus_t;
  typedef rocsolver_handle = rocsolverHandle_t;

  void hip_check_error();
  void hip_check(hipError_t sucess, std::string message="");
  void hipblas_check(hipblasStatus_t sucess, std::string message="");
  void hipsparse_check(hipsparseStatus_t sucess, std::string message="");
  void hiprand_check(hiprandStatus_t sucess, std::string message="");
  void hipsolver_check(hipsolverStatus_t sucess, std::string message="");
  hipblasOperation_t hipblasOperation(char A);
  hipsparseOperation_t hipsparseOperation(char A);

}

#endif
