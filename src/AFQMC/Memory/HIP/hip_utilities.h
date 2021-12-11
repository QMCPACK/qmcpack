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

#ifndef AFQMC_HIP_UTILITIES_HPP
#define AFQMC_HIP_UTILITIES_HPP

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <hip/hip_runtime.h>
#include "hipblas.h"
#include "hipsparse.h"
#include "rocsolver.h"
#include "rocrand/rocrand.h"

namespace qmc_hip
{
//  extern hiprandGenerator_t afqmc_curand_generator;
extern hipsparseMatDescr_t afqmc_hipsparse_matrix_descr;

extern std::vector<hipStream_t> afqmc_hip_streams;

// FDM: Temprorary hack to allow for easier grepping.
typedef rocsolver_status rocsolverStatus_t;
typedef rocsolver_status hipsolverStatus_t;
typedef rocsolver_handle hipsolverHandle_t;
typedef rocsolver_handle rocsolverHandle_t;
typedef rocblas_operation_ rocblasOperation_t;
typedef rocrand_status hiprandStatus_t;
typedef rocrand_generator hiprandGenerator_t;

void hip_check_error();
void hip_check(hipError_t success, std::string message = "");
void hipblas_check(hipblasStatus_t success, std::string message = "");
void hipsparse_check(hipsparseStatus_t success, std::string message = "");
void hiprand_check(hiprandStatus_t success, std::string message = "");
void hipsolver_check(hipsolverStatus_t success, std::string message = "");
hipblasOperation_t hipblasOperation(char A);
rocblasOperation_t rocblasOperation(char A);
hipsparseOperation_t hipsparseOperation(char A);

} // namespace qmc_hip

#endif
