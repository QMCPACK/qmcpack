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
#include <thrust/complex.h>
#include <thrust/device_ptr.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include "AFQMC/Numerics/detail/HIP/hip_kernel_utils.h"

namespace kernels
{
template<typename T>
struct axty_functor
{
  const T a;

  axty_functor(T _a) : a(_a) {}

  __host__ __device__ T operator()(const T& x, const T& y) const { return a * x * y; }
};

template<typename T>
inline static void kernel_axty(int n, T const alpha, T const* x, T* y)
{
  thrust::device_ptr<T const> x_(x);
  thrust::device_ptr<T> y_(y);
  thrust::transform(x_, x_ + n, y_, y_, axty_functor<T>(alpha));
}
template<typename T>
inline static void kernel_axty(int n, std::complex<T> const alpha, std::complex<T> const* x, std::complex<T>* y)
{
  thrust::device_ptr<thrust::complex<T> const> x_(reinterpret_cast<thrust::complex<T> const*>(x));
  thrust::device_ptr<thrust::complex<T>> y_(reinterpret_cast<thrust::complex<T>*>(y));
  thrust::transform(x_, x_ + n, y_, y_, axty_functor<thrust::complex<T>>(static_cast<thrust::complex<T>>(alpha)));
}

void axty(int n, float alpha, float const* x, float* y)
{
  kernel_axty(n, alpha, x, y);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void axty(int n, double alpha, double const* x, double* y)
{
  kernel_axty(n, alpha, x, y);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void axty(int n, std::complex<float> alpha, std::complex<float> const* x, std::complex<float>* y)
{
  kernel_axty(n, alpha, x, y);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}
void axty(int n, std::complex<double> alpha, std::complex<double> const* x, std::complex<double>* y)
{
  kernel_axty(n, alpha, x, y);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}


} // namespace kernels
