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
#include <iostream>
#include <hip/hip_runtime.h>

#include "hip_kernel_utils.h"
#include <rocrand/rocrand.h>

namespace qmc_hip
{
void hip_kernel_check(hipError_t success, std::string message)
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
void rocrand_check(rocrand_status success, std::string message)
{
  if (ROCRAND_STATUS_SUCCESS != success)
  {
    std::cerr << message << std::endl;
    std::cerr.flush();
    throw std::runtime_error(" Error code returned by hiprand. \n");
  }
}

} // namespace qmc_hip
