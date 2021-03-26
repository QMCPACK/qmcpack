//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include <memory>
#include <iostream>
#include "CUDA/cudaError.h"
#include "CUDA/CUDAallocator.hpp"
#include "CUDA/cuBLAS.hpp"
#include "detail/CUDA/cuBLAS_LU.hpp"
namespace qmcplusplus
{
namespace testing
{
struct CUDAHandles
{
  // CUDA specific variables
  cudaStream_t hstream;
  cublasHandle_t h_cublas;

  CUDAHandles()
  {
    cudaErrorCheck(cudaStreamCreate(&hstream), "cudaStreamCreate failed!");
    cublasErrorCheck(cublasCreate(&h_cublas), "cublasCreate failed!");
    cublasErrorCheck(cublasSetStream(h_cublas, hstream), "cublasSetStream failed!");
  }

  CUDAHandles(const CUDAHandles&) : CUDAHandles() {}

  ~CUDAHandles()
  {
    cublasErrorCheck(cublasDestroy(h_cublas), "cublasDestroy failed!");
    cudaErrorCheck(cudaStreamDestroy(hstream), "cudaStreamDestroy failed!");
  }
};
} // namespace testing

// This is in Configuration.h but we don't want all the dependencies that brings
#ifdef TEST_CASE
#ifdef QMC_COMPLEX
typedef ComplexApprox ValueApprox;
#else
typedef Approx ValueApprox;
#endif
#endif


TEST_CASE("cuBLAS_LU::computeLUDiag_batched", "[wavefunction][CUDA]")
{
  auto cuda_handles = std::make_unique<testing::CUDAHandles>();
  int n             = 4;
  int lda           = 4;
  auto& hstream     = cuda_handles->hstream;

  void* vpM = nullptr;
  cudaErrorCheck(cudaMallocHost(&vpM, sizeof(double) * 16), "cudaMallocHost failed");
  double* M = new (vpM) double[16]{7., 0.28571429, 0.71428571, 0.71428571,  5., 3.57142857, 0.12, -0.44,
                                   6., 6.28571429, -1.04,      -0.46153846, 6., 5.28571429, 3.08, 7.46153846};
  // M[0]      = 3;
  // M[1]      = 3;
  // M[2]      = 0;
  // M[3]      = 0;

  double* devM;
  cudaErrorCheck(cudaMalloc((void**)&devM, sizeof(double) * 16), "cudaMalloc failed");
  double** Ms;
  cudaErrorCheck(cudaMallocHost((void**)&Ms, sizeof(double*)), "cudaMallocHost failed");
  Ms = &devM;
  double** devMs;
  cudaErrorCheck(cudaMalloc((void**)&devMs, sizeof(double*)), "cudaMalloc failed");

  void* vpLUD = nullptr;
  cudaErrorCheck(cudaMallocHost(&vpLUD, sizeof(double) * 4), "cudaMallocHost failed");
  double* LU_diags = new (vpLUD) double[4]{1, 1, 1, 1};
  double* dev_LU_diags;
  cudaErrorCheck(cudaMalloc((void**)&dev_LU_diags, sizeof(double) * 4), "cudaMalloc failed");

  int batch_size = 1;

  cudaErrorCheck(cudaMemcpyAsync(devM, M, sizeof(double) * 16, cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying M to device");
  cudaErrorCheck(cudaMemcpyAsync(devMs, Ms, sizeof(double*), cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying Ms to device");
  // cudaErrorCheck(cudaMemcpyAsync(dev_LU_diags, LU_diags, sizeof(double) * 4, cudaMemcpyHostToDevice, hstream),
  //                "cudaMemcpyAsync failed copying LU_diags to device");


  cuBLAS_LU::computeLUDiag_batched(cuda_handles->hstream, n, lda, devMs, dev_LU_diags, batch_size);

  cudaErrorCheck(cudaMemcpyAsync(LU_diags, dev_LU_diags, sizeof(double) * 4, cudaMemcpyDeviceToHost, hstream),
                 "cudaMemcpyAsync failed copying LU_diags to device");
  cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!");

  std::cout << LU_diags[0] << " " << LU_diags[1] << " " << LU_diags[2] << " " << LU_diags[3] << '\n';

  CHECK(LU_diags[0] == M[0]);
  CHECK(LU_diags[1] == M[5]);
  CHECK(LU_diags[2] == M[10]);
  CHECK(LU_diags[3] == M[15]);
}

TEST_CASE("cuBLAS_LU::computeLogDet", "[wavefunction][CUDA]")
{
  auto cuda_handles = std::make_unique<testing::CUDAHandles>();
  int n             = 4;
  int lda           = 4;
  auto& hstream     = cuda_handles->hstream;

  void* vp_lu = nullptr;
  cudaErrorCheck(cudaMallocHost(&vp_lu, sizeof(double) * 4), "cudaMallocHost failed");
  double* lu = new (vp_lu) double[4]{7., 3.57142857, -1.04, 7.46153846};
  double* dev_lu;
  cudaErrorCheck(cudaMalloc((void**)&dev_lu, sizeof(double) * 4), "cudaMalloc failed");

  std::complex<double>* log_values;
  cudaErrorCheck(cudaMallocHost((void**)&log_values, sizeof(std::complex<double>) * 1), "cudaMallocHost failed");
  std::complex<double>* dev_log_values;
  cudaErrorCheck(cudaMalloc((void**)&dev_log_values, sizeof(std::complex<double>) * 1), "cudaMalloc failed");

  void* vp_pivots;
  cudaErrorCheck(cudaMallocHost((void**)&vp_pivots, sizeof(int) * 4), "cudaMallocHost failed");
  int* pivots = new (vp_pivots) int[4]{3, 3, 4, 4};
  int* dev_pivots;
  cudaErrorCheck(cudaMalloc((void**)&dev_pivots, sizeof(int) * 4), "cudaMalloc failed");

  cudaErrorCheck(cudaMemcpyAsync(dev_lu, lu, sizeof(double) * 4, cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying log_values to device");

  cudaErrorCheck(cudaMemcpyAsync(dev_pivots, pivots, sizeof(int) * 4, cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying log_values to device");

  int batch_size = 1;
  cuBLAS_LU::computeLogDet_batched(cuda_handles->hstream, n, dev_lu, dev_pivots, dev_log_values, batch_size);

  cudaErrorCheck(cudaMemcpyAsync(log_values, dev_log_values, sizeof(std::complex<double>) * 1, cudaMemcpyDeviceToHost,
                                 hstream),
                 "cudaMemcpyAsync failed copying log_values from device");
  cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!");

  CHECK(log_values[0] == ComplexApprox(std::complex<double>{5.267858159063328, 6.283185307179586}));
}

TEST_CASE("cuBLAS_LU::computeLogDet_complex", "[wavefunction][CUDA]")
{
  auto cuda_handles = std::make_unique<testing::CUDAHandles>();
  int n             = 4;
  int lda           = 4;
  auto& hstream     = cuda_handles->hstream;

  void* vp_lu = nullptr;
  cudaErrorCheck(cudaMallocHost(&vp_lu, sizeof(double) * 8), "cudaMallocHost failed");
  std::complex<double>* lu =
    new (vp_lu) std::complex<double>[8] { 7, 0.2, 3.699021  ,1.0657421 ,
      -0.52910995,0.97346985,  6.289567,2.0916762};

  std::complex<double>* dev_lu;
  cudaErrorCheck(cudaMalloc((void**)&dev_lu, sizeof(double) * 8), "cudaMalloc failed");

  std::complex<double>* log_values;
  cudaErrorCheck(cudaMallocHost((void**)&log_values, sizeof(std::complex<double>) * 1), "cudaMallocHost failed");
  std::complex<double>* dev_log_values;
  cudaErrorCheck(cudaMalloc((void**)&dev_log_values, sizeof(std::complex<double>) * 1), "cudaMalloc failed");

  void* vp_pivots;
  cudaErrorCheck(cudaMallocHost((void**)&vp_pivots, sizeof(int) * 4), "cudaMallocHost failed");
  int* pivots = new (vp_pivots) int[4]{3, 3, 4, 4};
  int* dev_pivots;
  cudaErrorCheck(cudaMalloc((void**)&dev_pivots, sizeof(int) * 4), "cudaMalloc failed");

  cudaErrorCheck(cudaMemcpyAsync(dev_lu, lu, sizeof(double) * 8, cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying log_values to device");

  cudaErrorCheck(cudaMemcpyAsync(dev_pivots, pivots, sizeof(int) * 4, cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying log_values to device");

  int batch_size = 1;
  cuBLAS_LU::computeLogDet_batched(cuda_handles->hstream, n, dev_lu, dev_pivots, dev_log_values, batch_size);

  cudaErrorCheck(cudaMemcpyAsync(log_values, dev_log_values, sizeof(std::complex<double>) * 1, cudaMemcpyDeviceToHost,
                                 hstream),
                 "cudaMemcpyAsync failed copying log_values from device");
  cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!");

  CHECK(log_values[0] == ComplexApprox(std::complex<double>{5.603777579195571, 0.1245249740607651}));
}


TEST_CASE("cuBLAS_LU::computeLogDet(batch=2)", "[wavefunction][CUDA]")
{
  auto cuda_handles = std::make_unique<testing::CUDAHandles>();
  int n             = 4;
  int lda           = 4;
  auto& hstream     = cuda_handles->hstream;

  void* vp_lu = nullptr;
  cudaErrorCheck(cudaMallocHost(&vp_lu, sizeof(double) * 8), "cudaMallocHost failed");
  double* lu   = new (vp_lu) double[4]{7., 3.57142857, -1.04, 7.46153846};
  void* vp_lu2 = (void*)(lu + 4);
  double* lu2  = new (vp_lu2) double[4]{7., 3.57142857, -1.04, 7.46153846};
  double* dev_lu;
  cudaErrorCheck(cudaMalloc((void**)&dev_lu, sizeof(double) * 8), "cudaMalloc failed");

  void* vp_log_values;
  cudaErrorCheck(cudaMallocHost(&vp_log_values, sizeof(std::complex<double>) * 2), "cudaMallocHost failed");
  // For values we expect zeroed as a side effect of a call we should poison them to test.
  std::complex<double>* log_values =
      reinterpret_cast<std::complex<double>*>(new (vp_log_values) double[4]{1.0, 1.0, 1.0, 1.0});

  std::complex<double>* dev_log_values;
  cudaErrorCheck(cudaMalloc((void**)&dev_log_values, sizeof(std::complex<double>) * 2), "cudaMalloc failed");

  void* vp_pivots;
  cudaErrorCheck(cudaMallocHost((void**)&vp_pivots, sizeof(int) * 8), "cudaMallocHost failed");
  int* pivots      = new (vp_pivots) int[4]{3, 3, 4, 4};
  void* vp_pivots2 = (void*)(pivots + 4);
  int* pivots2     = new (vp_pivots2) int[4]{3, 3, 4, 4};

  int* dev_pivots;
  cudaErrorCheck(cudaMalloc((void**)&dev_pivots, sizeof(int) * 8), "cudaMalloc failed");

  cudaErrorCheck(cudaMemcpyAsync(dev_lu, lu, sizeof(double) * 8, cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying log_values to device");

  cudaErrorCheck(cudaMemcpyAsync(dev_pivots, pivots, sizeof(int) * 8, cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying log_values to device");

  int batch_size = 2;

  // This is done artificially non zero initialize the log_values as they would be from a previous calculation.
  cudaErrorCheck(cudaMemcpyAsync(dev_log_values, log_values, sizeof(std::complex<double>) * 2, cudaMemcpyHostToDevice,
                                 hstream),
                 "cudaMemcpyAsync failed copying log_values to device");

  cuBLAS_LU::computeLogDet_batched(cuda_handles->hstream, n, dev_lu, dev_pivots, dev_log_values, batch_size);

  cudaErrorCheck(cudaMemcpyAsync(log_values, dev_log_values, sizeof(std::complex<double>) * 2, cudaMemcpyDeviceToHost,
                                 hstream),
                 "cudaMemcpyAsync failed copying log_values from device");
  cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!");

  CHECK(log_values[0] == ComplexApprox(std::complex<double>{5.267858159063328, 6.283185307179586}));
  CHECK(log_values[1] == ComplexApprox(std::complex<double>{5.267858159063328, 6.283185307179586}));
}


TEST_CASE("cuBLAS_LU::getrf_batched_complex", "[wavefunction][CUDA]")
{
  auto cuda_handles = std::make_unique<testing::CUDAHandles>();
  int n             = 4;
  int lda           = 4;
  auto& hstream     = cuda_handles->hstream;

  void* vpM = nullptr;
  cudaErrorCheck(cudaMallocHost(&vpM, sizeof(double) * 32), "cudaMallocHost failed");
  double* M = new (vpM) double[32]{2.,  0.1, 5.,  0.1, 7.0, 0.2,  5.0, 0.0,  5.0, 0.6, 2.0, 0.2, 5.0, 1.,   4.0, -0.1,
                                   8.0, 0.5, 2.0, 0.1, 6.0, -0.2, 4.0, -0.6, 7.0, 1.0, 8.0, 0.5, 6.0, -0.2, 8.0, -2.0};

  double* devM;
  cudaErrorCheck(cudaMalloc((void**)&devM, sizeof(double) * 32), "cudaMalloc failed");
  double** Ms;
  cudaErrorCheck(cudaMallocHost((void**)&Ms, sizeof(double*)), "cudaMallocHost failed");
  Ms = &devM;
  std::complex<double>** devMs;
  cudaErrorCheck(cudaMalloc((void**)&devMs, sizeof(double*)), "cudaMalloc failed");

  void* vp_pivots;
  cudaErrorCheck(cudaMallocHost((void**)&vp_pivots, sizeof(int) * 4), "cudaMallocHost failed");
  int* pivots = new (vp_pivots) int[4]{1, 1, 1, 1};
  int* dev_pivots;
  cudaErrorCheck(cudaMalloc((void**)&dev_pivots, sizeof(int) * 4), "cudaMalloc failed");

  void* vp_infos;
  cudaErrorCheck(cudaMallocHost((void**)&vp_infos, sizeof(int) * 4), "cudaMallocHost failed");
  int* infos = new (vp_infos) int[4]{1, 1, 1, 1};
  int* dev_infos;
  cudaErrorCheck(cudaMalloc((void**)&dev_infos, sizeof(int) * 4), "cudaMalloc failed");

  int batch_size = 1;

  cudaErrorCheck(cudaMemcpyAsync(devM, M, sizeof(double) * 32, cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying M to device");
  cudaErrorCheck(cudaMemcpyAsync(devMs, Ms, sizeof(double*), cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying Ms to device");

  cuBLAS_LU::computeGetrf_batched(cuda_handles->h_cublas, n, lda, devMs, dev_pivots, dev_infos, batch_size);

  cudaErrorCheck(cudaMemcpyAsync(M, devM, sizeof(double) * 32, cudaMemcpyDeviceToHost, hstream),
                 "cudaMemcpyAsync failed copying invM from device");
  cudaErrorCheck(cudaMemcpyAsync(infos, dev_infos, sizeof(int) * 4, cudaMemcpyDeviceToHost, hstream),
                 "cudaMemcpyAsync failed copying infos from device");
  cudaErrorCheck(cudaMemcpyAsync(pivots, dev_pivots, sizeof(int) * 4, cudaMemcpyDeviceToHost, hstream),
                 "cudaMemcpyAsync failed copying pivots from device");

  cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!");

  double lu[32]{7, 0.2,  0.285889, 0.00611746, 0.714111, -0.00611746, 0.713703, -0.0203915,
                5, 1,    3.57667,  -0.216476,  -0.43106, -0.161278,   0.126518, -0.191339,
                6, -0.2, 6.28344,  0.520473,   0.341156, 1.51726,     0.337414, 0.84877,
                6, -0.2, 5.28344,  1.02047,    5.82946,  1.97151,     2.56457,  -6.46618};

  int real_pivot[4]{3, 3, 3, 4};

  auto checkArray = [](auto* A, auto* B, int n) {
    for (int i = 0; i < n; ++i)
    {
      CHECK(A[i] == ValueApprox(B[i]));
    }
  };

  auto complexCheckArray = [](auto* A, auto* B, int n) {
    for (int i = 0; i < n; ++i)
    {
      CHECK(A[i] == ComplexApprox(B[i]));
    }
  };
  for (int ie = 0; ie < 16; ++ie)
    std::cout << M[ie * 2] << ", " << M[ie * 2 + 1] << ", ";
  std::cout << "\n";

  for (int ie = 0; ie < 4; ++ie)
    std::cout << M[ie * 8 + ie * 2] << ", " << M[ie * 8 + ie * 2 + 1] << ", ";
  std::cout << "\n";

  for (int ip = 0; ip < 4; ++ip)
    std::cout << pivots[ip] << ", ";
  std::cout << "\n";

  for (int ip = 0; ip < 4; ++ip)
    std::cout << infos[ip] << ", ";
  std::cout << "\n";

  complexCheckArray((std::complex<double>*)lu, (std::complex<double>*)M, 16);
  checkArray(real_pivot, pivots, 4);
}

TEST_CASE("cuBLAS_LU::getrf_batched", "[wavefunction][CUDA]")
{
  auto cuda_handles = std::make_unique<testing::CUDAHandles>();
  int n             = 4;
  int lda           = 4;
  auto& hstream     = cuda_handles->hstream;

  void* vpM = nullptr;
  cudaErrorCheck(cudaMallocHost(&vpM, sizeof(double) * 16), "cudaMallocHost failed");
  double* M = new (vpM) double[16]{2, 5, 7, 5, 5, 2, 5, 4, 8, 2, 6, 4, 7, 8, 6, 8};
  double* devM;
  cudaErrorCheck(cudaMalloc((void**)&devM, sizeof(double) * 16), "cudaMalloc failed");
  double** Ms;
  cudaErrorCheck(cudaMallocHost((void**)&Ms, sizeof(double*)), "cudaMallocHost failed");
  Ms = &devM;
  double** devMs;
  cudaErrorCheck(cudaMalloc((void**)&devMs, sizeof(double*)), "cudaMalloc failed");

  void* vp_pivots;
  cudaErrorCheck(cudaMallocHost((void**)&vp_pivots, sizeof(int) * 4), "cudaMallocHost failed");
  int* pivots = new (vp_pivots) int[4]{1, 1, 1, 1};
  int* dev_pivots;
  cudaErrorCheck(cudaMalloc((void**)&dev_pivots, sizeof(int) * 4), "cudaMalloc failed");

  void* vp_infos;
  cudaErrorCheck(cudaMallocHost((void**)&vp_infos, sizeof(int) * 4), "cudaMallocHost failed");
  int* infos = new (vp_infos) int[4]{1, 1, 1, 1};
  int* dev_infos;
  cudaErrorCheck(cudaMalloc((void**)&dev_infos, sizeof(int) * 4), "cudaMalloc failed");

  int batch_size = 1;

  cudaErrorCheck(cudaMemcpyAsync(devM, M, sizeof(double) * 16, cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying M to device");
  cudaErrorCheck(cudaMemcpyAsync(devMs, Ms, sizeof(double*), cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying Ms to device");

  cuBLAS_LU::computeGetrf_batched(cuda_handles->h_cublas, n, lda, devMs, dev_pivots, dev_infos, batch_size);


  cudaErrorCheck(cudaMemcpyAsync(M, devM, sizeof(double) * 16, cudaMemcpyDeviceToHost, hstream),
                 "cudaMemcpyAsync failed copying invM from device");
  cudaErrorCheck(cudaMemcpyAsync(infos, dev_infos, sizeof(int) * 4, cudaMemcpyDeviceToHost, hstream),
                 "cudaMemcpyAsync failed copying infos from device");
  cudaErrorCheck(cudaMemcpyAsync(pivots, dev_pivots, sizeof(int) * 4, cudaMemcpyDeviceToHost, hstream),
                 "cudaMemcpyAsync failed copying pivots from device");

  cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!");

  double lu[16]{7., 0.28571429, 0.71428571, 0.71428571,  5., 3.57142857, 0.12, -0.44,
                6., 6.28571429, -1.04,      -0.46153846, 6., 5.28571429, 3.08, 7.46153846};

  int real_pivot[4]{3, 3, 4, 4};

  auto checkArray = [](auto* A, auto* B, int n) {
    for (int i = 0; i < n; ++i)
    {
      CHECK(A[i] == ValueApprox(B[i]));
    }
  };

  checkArray(lu, M, 16);
  checkArray(real_pivot, pivots, 4);
}

TEST_CASE("cuBLAS_LU::getri_batched", "[wavefunction][CUDA]")
{
  auto cuda_handles = std::make_unique<testing::CUDAHandles>();
  int n             = 4;
  int lda           = 4;
  auto& hstream     = cuda_handles->hstream;

  void* vpM = nullptr;
  cudaErrorCheck(cudaMallocHost(&vpM, sizeof(double) * 16), "cudaMallocHost failed");
  double* M = new (vpM) double[16]{7., 0.28571429, 0.71428571, 0.71428571,  5., 3.57142857, 0.12, -0.44,
                                   6., 6.28571429, -1.04,      -0.46153846, 6., 5.28571429, 3.08, 7.46153846};
  double* devM;
  cudaErrorCheck(cudaMalloc((void**)&devM, sizeof(double) * 16), "cudaMalloc failed");
  double** Ms;
  cudaErrorCheck(cudaMallocHost((void**)&Ms, sizeof(double*)), "cudaMallocHost failed");
  Ms = &devM;
  double** devMs;
  cudaErrorCheck(cudaMalloc((void**)&devMs, sizeof(double*)), "cudaMalloc failed");

  void* vp_invM = nullptr;
  cudaErrorCheck(cudaMallocHost(&vp_invM, sizeof(double) * 16), "cudaMallocHost failed");
  double* invM = new (vp_invM) double[16]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  double* dev_invM;
  cudaErrorCheck(cudaMalloc((void**)&dev_invM, sizeof(double) * 16), "cudaMalloc failed");

  double** invMs;
  cudaErrorCheck(cudaMallocHost((void**)&invMs, sizeof(double*)), "cudaMallocHost failed");

  invMs = &dev_invM;
  double** dev_invMs;
  cudaErrorCheck(cudaMalloc((void**)&dev_invMs, sizeof(double*)), "cudaMalloc failed");

  void* vp_pivots;
  cudaErrorCheck(cudaMallocHost((void**)&vp_pivots, sizeof(int) * 4), "cudaMallocHost failed");
  int* pivots = new (vp_pivots) int[4]{3, 3, 4, 4};
  int* dev_pivots;
  cudaErrorCheck(cudaMalloc((void**)&dev_pivots, sizeof(int) * 4), "cudaMalloc failed");

  void* vp_infos;
  cudaErrorCheck(cudaMallocHost((void**)&vp_infos, sizeof(int) * 4), "cudaMallocHost failed");
  int* infos = new (vp_infos) int[4]{1, 1, 1, 1};
  int* dev_infos;
  cudaErrorCheck(cudaMalloc((void**)&dev_infos, sizeof(int) * 4), "cudaMalloc failed");

  int batch_size = 1;

  cudaErrorCheck(cudaMemcpyAsync(devM, M, sizeof(double) * 16, cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying M to device");
  cudaErrorCheck(cudaMemcpyAsync(devMs, Ms, sizeof(double*), cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying Ms to device");
  cudaErrorCheck(cudaMemcpyAsync(dev_invMs, invMs, sizeof(double*), cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying invMs to device");
  cudaErrorCheck(cudaMemcpyAsync(dev_pivots, pivots, sizeof(int) * 4, cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying pivots to device");

  cuBLAS_LU::computeGetri_batched(cuda_handles->h_cublas, n, lda, devMs, invMs, dev_pivots, dev_infos, batch_size);

  cudaErrorCheck(cudaMemcpyAsync(invM, dev_invM, sizeof(double) * 16, cudaMemcpyDeviceToHost, hstream),
                 "cudaMemcpyAsync failed copying invM from device");
  cudaErrorCheck(cudaMemcpyAsync(infos, dev_infos, sizeof(int) * 4, cudaMemcpyDeviceToHost, hstream),
                 "cudaMemcpyAsync failed copying infos from device");
  cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!");

  double invA[16]{-0.08247423, -0.26804124, 0.26804124, 0.05154639,  0.18556701,  -0.89690722, 0.39690722,  0.13402062,
                  0.24742268,  -0.19587629, 0.19587629, -0.15463918, -0.29896907, 1.27835052,  -0.77835052, 0.06185567};

  auto checkArray = [](double* A, double* B, int n) {
    for (int i = 0; i < n; ++i)
    {
      CHECK(A[i] == ValueApprox(B[i]));
    }
  };
  checkArray(invA, invM, 16);
}

} // namespace qmcplusplus
