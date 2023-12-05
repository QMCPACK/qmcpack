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
#include "AFQMC/Numerics/detail/HIP/hip_kernel_utils.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/hip_settings.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/copy_n_cast.hip.h"

namespace kernels
{
template<typename T, typename Q, typename Size>
__global__ void kernel_inplace_cast(Size n, thrust::complex<T>* A, thrust::complex<Q>* B)
{
  Size nb(blockDim.x);
  Size ni(threadIdx.x);
  thrust::complex<Q> Bi;
  if (sizeof(T) >= sizeof(Q))
  {
    for (Size i = 0; i < n; i += nb, ni += nb)
    {
      if (ni < n)
        Bi = static_cast<thrust::complex<Q>>(*(A + ni));
      __syncthreads();
      if (ni < n)
        *(B + ni) = Bi;
    }
  }
  else if (sizeof(T) < sizeof(Q))
  {
    ni = n - 1 - ni;
    for (Size i = 0; i < n; i += nb, ni -= nb)
    {
      if (ni >= 0)
        Bi = static_cast<thrust::complex<Q>>(*(A + ni));
      __syncthreads();
      if (ni >= 0)
        *(B + ni) = Bi;
    }
  }
}

void inplace_cast(unsigned long n, std::complex<float>* A, std::complex<double>* B)
{
  while (n > 4ul * 32ul)
  {
    unsigned long ni = n / 2uL; // number of elements to be copied in this iteration, ni <= n/2
    // copy_n_cast last ni elements
    copy_n_cast(A + (n - ni), int(ni), B + (n - ni));
    n -= ni;
  }
  hipLaunchKernelGGL(kernel_inplace_cast, dim3(1), dim3(32), 0, 0, n, reinterpret_cast<thrust::complex<float>*>(A),
                     reinterpret_cast<thrust::complex<double>*>(B));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void inplace_cast(unsigned long n, std::complex<double>* A, std::complex<float>* B)
{
  hipLaunchKernelGGL(kernel_inplace_cast, dim3(1), dim3(MAX_THREADS_PER_DIM), 0, 0, n,
                     reinterpret_cast<thrust::complex<double>*>(A), reinterpret_cast<thrust::complex<float>*>(B));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void inplace_cast(long n, std::complex<float>* A, std::complex<double>* B)
{
  while (n > long(4 * 32))
  {
    long ni = n / 2uL; // number of elements to be copied in this iteration, ni <= n/2
    // copy_n_cast last ni elements
    copy_n_cast(A + (n - ni), int(ni), B + (n - ni));
    n -= ni;
  }
  hipLaunchKernelGGL(kernel_inplace_cast, dim3(1), dim3(32), 0, 0, n, reinterpret_cast<thrust::complex<float>*>(A),
                     reinterpret_cast<thrust::complex<double>*>(B));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void inplace_cast(long n, std::complex<double>* A, std::complex<float>* B)
{
  hipLaunchKernelGGL(kernel_inplace_cast, dim3(1), dim3(MAX_THREADS_PER_DIM), 0, 0, n,
                     reinterpret_cast<thrust::complex<double>*>(A), reinterpret_cast<thrust::complex<float>*>(B));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

} // namespace kernels
