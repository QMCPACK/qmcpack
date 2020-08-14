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

#include <cassert>
#include <complex>
#include <vector>
#include <cuda.h>
#include <cuda_runtime.h>
#include "Numerics/helpers/determinant.hpp"

namespace kernels
{
double determinant_from_getrf_gpu(int N, double* m, int lda, int* piv)
{
  std::vector<double> m_(lda * N);
  std::vector<int> piv_(N);
  cudaMemcpy(m_.data(), m, lda * N * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(piv_.data(), piv, N * sizeof(int), cudaMemcpyDeviceToHost);
  return ma::determinant_from_getrf(N, m_.data(), lda, piv_.data());
}

std::complex<double> determinant_from_getrf_gpu(int N, std::complex<double>* m, int lda, int* piv)
{
  std::vector<std::complex<double>> m_(lda * N);
  std::vector<int> piv_(N);
  cudaMemcpy(m_.data(), m, lda * N * sizeof(std::complex<double>), cudaMemcpyDeviceToHost);
  cudaMemcpy(piv_.data(), piv, N * sizeof(int), cudaMemcpyDeviceToHost);
  return ma::determinant_from_getrf(N, m_.data(), lda, piv_.data());
}

void determinant_from_getrf_gpu(int N, double* m, int lda, int* piv, double* res)
{
  std::vector<double> m_(lda * N);
  std::vector<int> piv_(N);
  cudaMemcpy(m_.data(), m, lda * N * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(piv_.data(), piv, N * sizeof(int), cudaMemcpyDeviceToHost);
  double tmp;
  ma::determinant_from_getrf(N, m_.data(), lda, piv_.data(), &tmp);
  cudaMemcpy(res, &tmp, sizeof(double), cudaMemcpyHostToDevice);
}

void determinant_from_getrf_gpu(int N, std::complex<double>* m, int lda, int* piv, std::complex<double>* res)
{
  std::vector<std::complex<double>> m_(lda * N);
  std::vector<int> piv_(N);
  cudaMemcpy(m_.data(), m, lda * N * sizeof(std::complex<double>), cudaMemcpyDeviceToHost);
  cudaMemcpy(piv_.data(), piv, N * sizeof(int), cudaMemcpyDeviceToHost);
  std::complex<double> tmp;
  ma::determinant_from_getrf(N, m_.data(), lda, piv_.data(), &tmp);
  cudaMemcpy(res, &tmp, sizeof(std::complex<double>), cudaMemcpyHostToDevice);
}

//template __global__ void kernel_determinant_from_getrf<double>;
//template __global__ void kernel_determinant_from_getrf<std::complex<double>>;

} // namespace kernels
