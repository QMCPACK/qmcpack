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
#include <type_traits>
#include <cuda.h>
#include <thrust/complex.h>
#include <cuda_runtime.h>
//#include <thrust/system/cuda/detail/core/util.h>
#include "AFQMC/Memory/CUDA/cuda_utilities.h"
#include "AFQMC/Numerics/detail/CUDA/Kernels/cuda_settings.h"
#include "AFQMC/Numerics/detail/CUDA/Kernels/copy_n_cast.cuh"

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
    // MAM: using shared memory should speed this up considerably
    //      not worring about this now, since it is not currently used in the code
    //      and the implementation is tricky!
    //    Size nel(SM_SIZE_KB/sizeof(thrust::complex<Q>));
    //    __shared__ thrust::cuda_cub::core::uninitialized_array<thrust::complex<Q>,
    //                             SM_SIZE_KB*MAX_THREADS_PER_DIM/sizeof(thrust::complex<Q>)> cache;
    // copy and cast into the cache without need to sync, sync when the cache is full,
    // then copy to B without sync.
    // Alternatively, instead of using shared memory here, you can use the device buffers
    // and get buffer space from there if available. Then just call copy_n_cast followed by memcpy
    // to and from the buffer, instead of using hand written kernels like this.
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
      if (std::is_unsigned<Size>::value || ni >= 0)
        Bi = static_cast<thrust::complex<Q>>(*(A + ni));
      __syncthreads();
      if (std::is_unsigned<Size>::value || ni >= 0)
        *(B + ni) = Bi;
      //      __syncthreads();
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
  kernel_inplace_cast<<<1, 32>>>(n, reinterpret_cast<thrust::complex<float>*>(A),
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
  while (n > long(4 * 32))
  {
    long ni = n / long(2); // number of elements to be copied in this iteration, ni <= n/2
    // copy_n_cast last ni elements
    copy_n_cast(A + (n - ni), int(ni), B + (n - ni));
    n -= ni;
  }
  kernel_inplace_cast<<<1, 32>>>(n, reinterpret_cast<thrust::complex<float>*>(A),
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
