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
#include <hip/hip_runtime.h>
#include <rocrand/rocrand.h>
#include "AFQMC/Numerics/detail/HIP/Kernels/hip_settings.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/zero_complex_part.hip.h"
#include "AFQMC/Numerics/detail/HIP/hip_kernel_utils.h"

namespace kernels
{
void sampleGaussianRNG(double* V, int n, rocrand_generator& gen)
{
  qmc_hip::rocrand_check(rocrand_generate_normal_double(gen, V, n, 0.0, 1.0), "rocrand_generate_normal_double");
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

// Convert to double if really necessary
void sampleGaussianRNG(float* V, int n, rocrand_generator& gen)
{
  qmc_hip::rocrand_check(rocrand_generate_normal(gen, V, n, float(0.0), float(1.0)), "rocrand_generate_normal");
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void sampleGaussianRNG(std::complex<double>* V, int n, rocrand_generator& gen)
{
  qmc_hip::rocrand_check(rocrand_generate_normal_double(gen, reinterpret_cast<double*>(V), 2 * n, 0.0, 1.0),
                         "rocrand_generate_normal_double");
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
  // hack hack hack!!!
  kernels::zero_complex_part(n, V);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void sampleGaussianRNG(std::complex<float>* V, int n, rocrand_generator& gen)
{
  qmc_hip::rocrand_check(rocrand_generate_normal(gen, reinterpret_cast<float*>(V), 2 * n, float(0.0), float(1.0)),
                         "rocrand_generate_normal");
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
  // hack hack hack!!!
  kernels::zero_complex_part(n, V);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

} // namespace kernels
