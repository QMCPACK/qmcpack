//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Youngjun Lee, leey@anl.gov, Argonne National Laboratory
//
// File created by: Youngjun Lee, leey@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include <DualAllocatorAliases.hpp>
#include <AccelBLAS.hpp>
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include <CPU/BLAS.hpp>
#include "batched_blas_perf.hpp"

namespace qmcplusplus
{

template<PlatformKind PL>
void benchmark_gemm_cases()
{
  const int M_b   = 29;
  const int N_b   = 31;
  const int K_b   = 23;
  const int BATCH = 64;

  // Batched Test
  benchmark::perf_gemm_batched<PL, float>('N', 'N', M_b, N_b, K_b, BATCH);
  benchmark::perf_gemm_batched<PL, double>('N', 'N', M_b, N_b, K_b, BATCH);
  benchmark::perf_gemm_batched<PL, float>('T', 'N', M_b, N_b, K_b, BATCH);
  benchmark::perf_gemm_batched<PL, double>('T', 'N', M_b, N_b, K_b, BATCH);
  benchmark::perf_gemm_batched<PL, float>('N', 'T', M_b, N_b, K_b, BATCH);
  benchmark::perf_gemm_batched<PL, double>('N', 'T', M_b, N_b, K_b, BATCH);
  benchmark::perf_gemm_batched<PL, float>('T', 'T', M_b, N_b, K_b, BATCH);
  benchmark::perf_gemm_batched<PL, double>('T', 'T', M_b, N_b, K_b, BATCH);
#if defined(QMC_COMPLEX)
  benchmark::perf_gemm_batched<PL, std::complex<float>>('N', 'N', M_b, N_b, K_b, BATCH);
  benchmark::perf_gemm_batched<PL, std::complex<double>>('N', 'N', M_b, N_b, K_b, BATCH);
  benchmark::perf_gemm_batched<PL, std::complex<float>>('T', 'N', M_b, N_b, K_b, BATCH);
  benchmark::perf_gemm_batched<PL, std::complex<double>>('T', 'N', M_b, N_b, K_b, BATCH);
  benchmark::perf_gemm_batched<PL, std::complex<float>>('N', 'T', M_b, N_b, K_b, BATCH);
  benchmark::perf_gemm_batched<PL, std::complex<double>>('N', 'T', M_b, N_b, K_b, BATCH);
  benchmark::perf_gemm_batched<PL, std::complex<float>>('T', 'T', M_b, N_b, K_b, BATCH);
  benchmark::perf_gemm_batched<PL, std::complex<double>>('T', 'T', M_b, N_b, K_b, BATCH);
#endif
}

template<PlatformKind PL>
void benchmark_gemv_cases()
{
  const int M_b   = 137;
  const int N_b   = 79;
  const int BATCH = 64;

  // Batched Test
  benchmark::perf_gemv_batched<PL, float>('N', M_b, N_b, BATCH);
  benchmark::perf_gemv_batched<PL, double>('N', M_b, N_b, BATCH);
  benchmark::perf_gemv_batched<PL, float>('T', M_b, N_b, BATCH);
  benchmark::perf_gemv_batched<PL, double>('T', M_b, N_b, BATCH);
#if defined(QMC_COMPLEX)
  benchmark::perf_gemv_batched<PL, std::complex<float>>('N', M_b, N_b, BATCH);
  benchmark::perf_gemv_batched<PL, std::complex<double>>('N', M_b, N_b, BATCH);
  benchmark::perf_gemv_batched<PL, std::complex<float>>('T', M_b, N_b, BATCH);
  benchmark::perf_gemv_batched<PL, std::complex<double>>('T', M_b, N_b, BATCH);
#endif
}

template<PlatformKind PL>
void benchmark_ger_cases()
{
  const int M_b   = 137;
  const int N_b   = 79;
  const int BATCH = 64;

  // Batched Test
  benchmark::perf_ger_batched<PL, float>(M_b, N_b, BATCH);
  benchmark::perf_ger_batched<PL, double>(M_b, N_b, BATCH);
#if defined(QMC_COMPLEX)
  benchmark::perf_ger_batched<PL, std::complex<float>>(M_b, N_b, BATCH);
  benchmark::perf_ger_batched<PL, std::complex<double>>(M_b, N_b, BATCH);
#endif
}


TEST_CASE("AccelBLAS batched benchmark", "[BLAS][benchmark]")
{
  SECTION("gemm_batched")
  {
#if defined(ENABLE_CUDA)
    std::cout << "Benchmarking gemm_batched<PlatformKind::CUDA>" << std::endl;
    benchmark_gemm_cases<PlatformKind::CUDA>();
#endif
#if defined(ENABLE_SYCL)
    std::cout << "Benchmarking gemm_batched<PlatformKind::SYCL>" << std::endl;
    benchmark_gemm_cases<PlatformKind::SYCL>();
#endif
#if defined(ENABLE_OFFLOAD)
    std::cout << "Benchmarking gemm_batched<PlatformKind::OMPTARGET>" << std::endl;
    benchmark_gemm_cases<PlatformKind::OMPTARGET>();
#endif
  }

  SECTION("gemv_batched")
  {
#if defined(ENABLE_CUDA)
    std::cout << "Benchmarking gemv_batched<PlatformKind::CUDA>" << std::endl;
    benchmark_gemv_cases<PlatformKind::CUDA>();
#endif
#if defined(ENABLE_SYCL)
    std::cout << "Benchmarking gemv_batched<PlatformKind::SYCL>" << std::endl;
    benchmark_gemv_cases<PlatformKind::SYCL>();
#endif
#if defined(ENABLE_OFFLOAD)
    std::cout << "Benchmarking gemv_batched<PlatformKind::OMPTARGET>" << std::endl;
    benchmark_gemv_cases<PlatformKind::OMPTARGET>();
#endif
  }

  SECTION("ger_batched")
  {
#if defined(ENABLE_CUDA)
    std::cout << "Benchmarking ger_batched<PlatformKind::CUDA>" << std::endl;
    benchmark_ger_cases<PlatformKind::CUDA>();
#endif
#if defined(ENABLE_SYCL)
    std::cout << "Benchmarking ger_batched<PlatformKind::SYCL>" << std::endl;
    benchmark_ger_cases<PlatformKind::SYCL>();
#endif
#if defined(ENABLE_OFFLOAD)
    std::cout << "Benchmarking ger_batched<PlatformKind::OMPTARGET>" << std::endl;
    benchmark_ger_cases<PlatformKind::OMPTARGET>();
#endif
  }
}

} // namespace qmcplusplus
