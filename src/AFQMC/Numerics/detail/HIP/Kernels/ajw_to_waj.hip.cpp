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
#include <thrust/complex.h>
#include <hip/hip_runtime.h>
#include "AFQMC/Numerics/detail/HIP/Kernels/hip_settings.h"
#include "AFQMC/Numerics/detail/HIP/hip_kernel_utils.h"

namespace kernels
{
// very sloppy, needs improvement!!!!
template<typename T, typename T2>
__global__ void kernel_ajw_to_waj(int na, int nj, int nw, int inca, T const* A, T2* B)
{
  int a = blockIdx.x;
  if (a >= na)
    return;
  T const* A_(A + inca * a);
  int i   = threadIdx.x;
  int njw = nj * nw;
  while (i < njw)
  {
    int j                    = i / nw;
    int w                    = i % nw;
    B[(w * na + a) * nj + j] = static_cast<T2>(A_[i]);
    i += blockDim.x;
  }
}

template<typename T, typename T2>
__global__ void kernel_ajw_to_waj(int na, int nj, int nw, int inca, thrust::complex<T> const* A, thrust::complex<T2>* B)
{
  int a = blockIdx.x;
  if (a >= na)
    return;
  thrust::complex<T> const* A_(A + inca * a);
  int i   = threadIdx.x;
  int njw = nj * nw;
  while (i < njw)
  {
    int j                    = i / nw;
    int w                    = i % nw;
    B[(w * na + a) * nj + j] = static_cast<thrust::complex<T2>>(A_[i]);
    i += blockDim.x;
  }
}

template<typename T, typename T2>
__global__ void kernel_transpose_wabn_to_wban(int nwalk,
                                              int na,
                                              int nb,
                                              int nchol,
                                              thrust::complex<T> const* Tab,
                                              thrust::complex<T2>* Tba)
{
  int w = blockIdx.x;
  int a = blockIdx.y;
  int b = blockIdx.z;

  thrust::complex<T> const* A_(Tab + ((w * na + a) * nb + b) * nchol);
  thrust::complex<T2>* B_(Tba + ((w * nb + b) * na + a) * nchol);
  int i = threadIdx.x;
  while (i < nchol)
  {
    B_[i] = static_cast<thrust::complex<T2>>(A_[i]);
    i += blockDim.x;
  }
}

void ajw_to_waj(int na, int nj, int nw, int inca, double const* A, double* B)
{
  hipLaunchKernelGGL(kernel_ajw_to_waj, dim3(na), dim3(128), 0, 0, na, nj, nw, inca, A, B);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void ajw_to_waj(int na, int nj, int nw, int inca, float const* A, float* B)
{
  hipLaunchKernelGGL(kernel_ajw_to_waj, dim3(na), dim3(128), 0, 0, na, nj, nw, inca, A, B);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void ajw_to_waj(int na, int nj, int nw, int inca, std::complex<double> const* A, std::complex<double>* B)
{
  hipLaunchKernelGGL(kernel_ajw_to_waj, dim3(na), dim3(128), 0, 0, na, nj, nw, inca,
                     reinterpret_cast<thrust::complex<double> const*>(A),
                     reinterpret_cast<thrust::complex<double>*>(B));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void ajw_to_waj(int na, int nj, int nw, int inca, std::complex<float> const* A, std::complex<float>* B)
{
  hipLaunchKernelGGL(kernel_ajw_to_waj, dim3(na), dim3(128), 0, 0, na, nj, nw, inca,
                     reinterpret_cast<thrust::complex<float> const*>(A), reinterpret_cast<thrust::complex<float>*>(B));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void transpose_wabn_to_wban(int nwalk,
                            int na,
                            int nb,
                            int nchol,
                            std::complex<double> const* Tab,
                            std::complex<double>* Tba)
{
  dim3 grid_dim(nwalk, na, nb);
  hipLaunchKernelGGL(kernel_transpose_wabn_to_wban, dim3(grid_dim), dim3(32), 0, 0, nwalk, na, nb, nchol,
                     reinterpret_cast<thrust::complex<double> const*>(Tab),
                     reinterpret_cast<thrust::complex<double>*>(Tba));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void transpose_wabn_to_wban(int nwalk,
                            int na,
                            int nb,
                            int nchol,
                            std::complex<float> const* Tab,
                            std::complex<float>* Tba)
{
  dim3 grid_dim(nwalk, na, nb);
  hipLaunchKernelGGL(kernel_transpose_wabn_to_wban, dim3(grid_dim), dim3(32), 0, 0, nwalk, na, nb, nchol,
                     reinterpret_cast<thrust::complex<float> const*>(Tab),
                     reinterpret_cast<thrust::complex<float>*>(Tba));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void transpose_wabn_to_wban(int nwalk,
                            int na,
                            int nb,
                            int nchol,
                            std::complex<double> const* Tab,
                            std::complex<float>* Tba)
{
  dim3 grid_dim(nwalk, na, nb);
  hipLaunchKernelGGL(kernel_transpose_wabn_to_wban, dim3(grid_dim), dim3(32), 0, 0, nwalk, na, nb, nchol,
                     reinterpret_cast<thrust::complex<double> const*>(Tab),
                     reinterpret_cast<thrust::complex<float>*>(Tba));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void transpose_wabn_to_wban(int nwalk,
                            int na,
                            int nb,
                            int nchol,
                            std::complex<float> const* Tab,
                            std::complex<double>* Tba)
{
  dim3 grid_dim(nwalk, na, nb);
  hipLaunchKernelGGL(kernel_transpose_wabn_to_wban, dim3(grid_dim), dim3(32), 0, 0, nwalk, na, nb, nchol,
                     reinterpret_cast<thrust::complex<float> const*>(Tab),
                     reinterpret_cast<thrust::complex<double>*>(Tba));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

} // namespace kernels
