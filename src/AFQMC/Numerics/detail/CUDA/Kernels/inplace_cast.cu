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
#include <cuda.h>
#include <thrust/complex.h>
#include <cuda_runtime.h>
#include "AFQMC/Memory/CUDA/cuda_utilities.h"
#include "AFQMC/Numerics/detail/CUDA/Kernels/cuda_settings.h"

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
    // this is wrong if the compiler eliminates the intermediate Bi through optimization
    // will it???
    for (Size i = 0; i < n; i += nb, ni += nb)
    {
      if (ni < n)
        Bi = static_cast<thrust::complex<Q>>(*(A + ni));
      __syncthreads();
      if (ni < n)
        *(B + ni) = Bi;
      //      __syncthreads();
    }
  }
  else if (sizeof(T) < sizeof(Q))
  {
    assert(sizeof(T) * 2 <= sizeof(Q));
    ni = n - 1 - ni;
    for (Size i = 0; i < n; i += nb, ni -= nb)
    {
      if (ni >= 0)
        Bi = static_cast<thrust::complex<Q>>(*(A + ni));
      __syncthreads();
      if (ni >= 0)
        *(B + ni) = Bi;
      //      __syncthreads();
    }
  }
}

void inplace_cast(unsigned long n, std::complex<float>* A, std::complex<double>* B)
{
  kernel_inplace_cast<<<1, MAX_THREADS_PER_DIM>>>(n, reinterpret_cast<thrust::complex<float>*>(A),
                                                  reinterpret_cast<thrust::complex<double>*>(B));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void inplace_cast(unsigned long n, std::complex<double>* A, std::complex<float>* B)
{
  kernel_inplace_cast<<<1, MAX_THREADS_PER_DIM>>>(n, reinterpret_cast<thrust::complex<double>*>(A),
                                                  reinterpret_cast<thrust::complex<float>*>(B));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void inplace_cast(long n, std::complex<float>* A, std::complex<double>* B)
{
  kernel_inplace_cast<<<1, MAX_THREADS_PER_DIM>>>(n, reinterpret_cast<thrust::complex<float>*>(A),
                                                  reinterpret_cast<thrust::complex<double>*>(B));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void inplace_cast(long n, std::complex<double>* A, std::complex<float>* B)
{
  kernel_inplace_cast<<<1, MAX_THREADS_PER_DIM>>>(n, reinterpret_cast<thrust::complex<double>*>(A),
                                                  reinterpret_cast<thrust::complex<float>*>(B));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

} // namespace kernels
