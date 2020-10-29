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

#include <complex>
#include <thrust/complex.h>
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>
#include "AFQMC/Numerics/detail/HIP/Kernels/strided_range.hpp"
#include "AFQMC/Numerics/detail/HIP/Kernels/strided_2Drange.hpp"
#include "AFQMC/Numerics/detail/HIP/hip_kernel_utils.h"

namespace kernels
{
/*
 * There is a bug in HIP9.2 when calling reduce with thrust::complex.
 * Temporary hack until bug is fixed...
 */


// WRITE A KERNEL!!!!

double sum(int n, double const* x, int incx)
{
  thrust::device_ptr<double const> x_(x);
  strided_range<thrust::device_ptr<double const>> strided(x_, x_ + n * incx, incx);
  double res = thrust::reduce(strided.begin(), strided.end());
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
  return res;
}

std::complex<double> sum(int n, std::complex<double> const* x, int incx)
{
  thrust::device_ptr<double const> x_(reinterpret_cast<double const*>(x));
  strided_range<thrust::device_ptr<double const>> Rstrided(x_, x_ + 2 * n * incx, 2 * incx);
  double R = thrust::reduce(Rstrided.begin(), Rstrided.end());
  strided_range<thrust::device_ptr<double const>> Istrided(x_ + 1, x_ + 1 + 2 * n * incx, 2 * incx);
  double I = thrust::reduce(Istrided.begin(), Istrided.end());
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
  return std::complex<double>(R, I);
}

float sum(int n, float const* x, int incx)
{
  thrust::device_ptr<float const> x_(x);
  strided_range<thrust::device_ptr<float const>> strided(x_, x_ + n * incx, incx);
  float res = thrust::reduce(strided.begin(), strided.end());
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
  return res;
}

std::complex<float> sum(int n, std::complex<float> const* x, int incx)
{
  thrust::device_ptr<float const> x_(reinterpret_cast<float const*>(x));
  strided_range<thrust::device_ptr<float const>> Rstrided(x_, x_ + 2 * n * incx, 2 * incx);
  float R = thrust::reduce(Rstrided.begin(), Rstrided.end());
  strided_range<thrust::device_ptr<float const>> Istrided(x_ + 1, x_ + 1 + 2 * n * incx, 2 * incx);
  float I = thrust::reduce(Istrided.begin(), Istrided.end());
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
  return std::complex<float>(R, I);
}

double sum(int m, int n, double const* x, int lda)
{
  thrust::device_ptr<double const> x_(x);
  strided_2Drange<thrust::device_ptr<double const>> strided(x_, x_ + n * lda, lda, m);
  double res = thrust::reduce(strided.begin(), strided.end());
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
  return res;
}

std::complex<double> sum(int m, int n, std::complex<double> const* x, int lda)
{
  std::complex<double> res;
  for (int i = 0; i < m; i++)
    res += sum(n, x + i, lda);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
  return res;
}

float sum(int m, int n, float const* x, int lda)
{
  thrust::device_ptr<float const> x_(x);
  strided_2Drange<thrust::device_ptr<float const>> strided(x_, x_ + n * lda, lda, m);
  return thrust::reduce(strided.begin(), strided.end());
}

std::complex<float> sum(int m, int n, std::complex<float> const* x, int lda)
{
  std::complex<float> res;
  for (int i = 0; i < m; i++)
    res += sum(n, x + i, lda);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
  return res;
}

} // namespace kernels
