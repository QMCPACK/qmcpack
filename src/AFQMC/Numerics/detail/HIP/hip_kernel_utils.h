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

#ifndef AFQMC_HIP_KERNEL_UTILITIES_HPP
#define AFQMC_HIP_KERNEL_UTILITIES_HPP

#include <cassert>
#include <hip/hip_runtime.h>

#include "hip_kernel_utils.h"
#include <rocrand/rocrand.h>

namespace qmc_hip
{
void hip_kernel_check(hipError_t success, std::string message = "");
void rocrand_check(rocrand_status success, std::string message = "");
} // namespace qmc_hip

#endif
