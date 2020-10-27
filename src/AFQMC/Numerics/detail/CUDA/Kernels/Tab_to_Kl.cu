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
#include "AFQMC/Numerics/detail/CUDA/Kernels/cuda_settings.h"
#define ENABLE_CUDA 1
#include "AFQMC/Memory/CUDA/cuda_utilities.h"

namespace kernels
{
// Tab [nbatch][nwalk][nocc][nocc][nchol]
// Klr [nwalk][2*nchol_tot]
template<typename T>
__global__ void kernel_Tab_to_Kl(int nwalk, int nocc, int nchol, thrust::complex<T> const* Tab, thrust::complex<T>* Kl)
{
  int w = blockIdx.x;
  if (w < nwalk)
  {
    for (int a = 0; a < nocc; a++)
    {
      thrust::complex<T> const* Tab_(Tab + ((w * nocc + a) * nocc + a) * nchol);
      thrust::complex<T>* Kl_(Kl + w * nchol);
      int c = threadIdx.x;
      while (c < nchol)
      {
        Kl_[c] += Tab_[c];
        c += blockDim.x;
      }
    }
  }
}

template<typename T>
__global__ void kernel_Tanb_to_Kl(int nwalk,
                                  int nocc,
                                  int nchol,
                                  int nchol_tot,
                                  thrust::complex<T> const* Tab,
                                  thrust::complex<T>* Kl)
{
  int w = blockIdx.x;
  if (w < nwalk)
  {
    for (int a = 0; a < nocc; a++)
    {
      thrust::complex<T> const* Tab_(Tab + ((w * nocc + a) * nocc) * nchol + a);
      thrust::complex<T>* Kl_(Kl + w * nchol_tot);
      int c = threadIdx.x;
      while (c < nchol)
      {
        Kl_[c] += Tab_[c * nocc];
        c += blockDim.x;
      }
    }
  }
}

void Tab_to_Kl(int nwalk, int nocc, int nchol, std::complex<double> const* Tab, std::complex<double>* Kl)
{
  dim3 grid_dim(nwalk, 1, 1);
  int nthr = std::min(256, nchol);
  kernel_Tab_to_Kl<<<grid_dim, nthr>>>(nwalk, nocc, nchol, reinterpret_cast<thrust::complex<double> const*>(Tab),
                                       reinterpret_cast<thrust::complex<double>*>(Kl));
  qmc_cuda::cuda_check(cudaGetLastError(), "Tab_to_Kl");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "Tab_to_Kl");
}

void Tab_to_Kl(int nwalk, int nocc, int nchol, std::complex<float> const* Tab, std::complex<float>* Kl)
{
  dim3 grid_dim(nwalk, 1, 1);
  int nthr = std::min(256, nchol);
  kernel_Tab_to_Kl<<<grid_dim, nthr>>>(nwalk, nocc, nchol, reinterpret_cast<thrust::complex<float> const*>(Tab),
                                       reinterpret_cast<thrust::complex<float>*>(Kl));
  qmc_cuda::cuda_check(cudaGetLastError(), "Tab_to_Kl");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "Tab_to_Kl");
}


void Tanb_to_Kl(int nwalk,
                int nocc,
                int nchol,
                int nchol_tot,
                std::complex<double> const* Tab,
                std::complex<double>* Kl)
{
  dim3 grid_dim(nwalk, 1, 1);
  int nthr = std::min(256, nchol);
  kernel_Tanb_to_Kl<<<grid_dim, nthr>>>(nwalk, nocc, nchol, nchol_tot,
                                        reinterpret_cast<thrust::complex<double> const*>(Tab),
                                        reinterpret_cast<thrust::complex<double>*>(Kl));
  qmc_cuda::cuda_check(cudaGetLastError(), "Tab_to_Kl");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "Tab_to_Kl");
}

void Tanb_to_Kl(int nwalk, int nocc, int nchol, int nchol_tot, std::complex<float> const* Tab, std::complex<float>* Kl)
{
  dim3 grid_dim(nwalk, 1, 1);
  int nthr = std::min(256, nchol);
  kernel_Tanb_to_Kl<<<grid_dim, nthr>>>(nwalk, nocc, nchol, nchol_tot,
                                        reinterpret_cast<thrust::complex<float> const*>(Tab),
                                        reinterpret_cast<thrust::complex<float>*>(Kl));
  qmc_cuda::cuda_check(cudaGetLastError(), "Tab_to_Kl");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "Tab_to_Kl");
}

} // namespace kernels
