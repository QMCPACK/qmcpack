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
#include <vector>
#include <hip/hip_runtime.h>
#include <hip/hip_runtime.h>
#include "Numerics/helpers/determinant.hpp"

namespace kernels
{
double determinant_from_getrf_gpu(int N, double* m, int lda, int* piv)
{
  std::vector<double> m_(lda * N);
  std::vector<int> piv_(N);
  hipMemcpy(m_.data(), m, lda * N * sizeof(double), hipMemcpyDeviceToHost);
  hipMemcpy(piv_.data(), piv, N * sizeof(int), hipMemcpyDeviceToHost);
  return ma::determinant_from_getrf(N, m_.data(), lda, piv_.data());
}

std::complex<double> determinant_from_getrf_gpu(int N, std::complex<double>* m, int lda, int* piv)
{
  std::vector<std::complex<double>> m_(lda * N);
  std::vector<int> piv_(N);
  hipMemcpy(m_.data(), m, lda * N * sizeof(std::complex<double>), hipMemcpyDeviceToHost);
  hipMemcpy(piv_.data(), piv, N * sizeof(int), hipMemcpyDeviceToHost);
  return ma::determinant_from_getrf(N, m_.data(), lda, piv_.data());
}

void determinant_from_getrf_gpu(int N, double* m, int lda, int* piv, double* res)
{
  std::vector<double> m_(lda * N);
  std::vector<int> piv_(N);
  hipMemcpy(m_.data(), m, lda * N * sizeof(double), hipMemcpyDeviceToHost);
  hipMemcpy(piv_.data(), piv, N * sizeof(int), hipMemcpyDeviceToHost);
  double tmp;
  ma::determinant_from_getrf(N, m_.data(), lda, piv_.data(), &tmp);
  hipMemcpy(res, &tmp, sizeof(double), hipMemcpyHostToDevice);
}

void determinant_from_getrf_gpu(int N, std::complex<double>* m, int lda, int* piv, std::complex<double>* res)
{
  std::vector<std::complex<double>> m_(lda * N);
  std::vector<int> piv_(N);
  hipMemcpy(m_.data(), m, lda * N * sizeof(std::complex<double>), hipMemcpyDeviceToHost);
  hipMemcpy(piv_.data(), piv, N * sizeof(int), hipMemcpyDeviceToHost);
  std::complex<double> tmp;
  ma::determinant_from_getrf(N, m_.data(), lda, piv_.data(), &tmp);
  hipMemcpy(res, &tmp, sizeof(std::complex<double>), hipMemcpyHostToDevice);
}

//template __global__ void kernel_determinant_from_getrf<double>;
//template __global__ void kernel_determinant_from_getrf<std::complex<double>>;

} // namespace kernels
