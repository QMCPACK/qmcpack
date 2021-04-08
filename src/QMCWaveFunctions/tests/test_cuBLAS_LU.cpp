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
#include "CUDA/cuBLAS.hpp"
#include "CUDA/CUDAfill.hpp"
#include "Utilities/for_testing/MatrixAccessor.hpp"
#include "Utilities/for_testing/checkMatrix.hpp"
#include "detail/CUDA/cuBLAS_LU.hpp"

/** \file
 *
 *  These are unit tests for the low level LU factorization used by the full inversion and
 *  calculation of log determinant for dirac determinants. Fundamental testing of these kernels
 *  requires full knowledge of the memory layout and data movement, As such OhmmsMatrices and
 *  custom allocators are not used.  They have their own unit tests (Hopefully!) This is also documentation
 *  of how these calls expect the memory handed to them to look.  Please leave this intact.
 *  Someday those container abstractions will change, if inversion breaks and this stil works you
 *  will have a fighting chance to know how to change these routines or fix the bug you introduced in the
 *  higher level abstractions.
 */
namespace qmcplusplus
{
namespace testing
{
/** Doesn't depend on the resource managment scheme thats out of scope for unit tests */
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

/** Single double computeLogDet */
TEST_CASE("cuBLAS_LU::computeLogDet", "[wavefunction][CUDA]")
{
  auto cuda_handles = std::make_unique<testing::CUDAHandles>();
  int n             = 4;
  int lda           = 4;
  int batch_size    = 1;
  auto& hstream     = cuda_handles->hstream;

  // The calls in cuBLAS_LU still care about type
  // being very careful and clear here.
  void* vp_lu = nullptr;
  cudaCheck(cudaMallocHost(&vp_lu, sizeof(double) * 16));
  // This LU was calculated in numpy from the first M in the getrf test
  double* lu = new (vp_lu) double[16]{7., 0.28571429, 0.71428571, 0.71428571,  5., 3.57142857, 0.12, -0.44,
                                      6., 6.28571429, -1.04,      -0.46153846, 6., 5.28571429, 3.08, 7.46153846};
  double* dev_lu;
  cudaCheck(cudaMalloc((void**)&dev_lu, sizeof(double) * 16));
  void* vp_lus;
  cudaCheck(cudaMallocHost(&vp_lus, sizeof(double**)));
  double** lus = new (vp_lus) double* [1] { nullptr };
  lus[0]       = dev_lu;
  double** dev_lus;
  cudaCheck(cudaMalloc((void**)&dev_lus, sizeof(double**)));

  std::complex<double>* log_values;
  cudaCheck(cudaMallocHost((void**)&log_values, sizeof(std::complex<double>) * 1));
  std::complex<double>* dev_log_values;
  cudaCheck(cudaMalloc((void**)&dev_log_values, sizeof(std::complex<double>) * 1));
  CUDAfill_n(dev_log_values, batch_size, {0,0});

  void* vp_pivots;
  cudaCheck(cudaMallocHost(&vp_pivots, sizeof(int) * 4));
  int* pivots = new (vp_pivots) int[4]{3, 3, 4, 4};
  int* dev_pivots;
  cudaCheck(cudaMalloc((void**)&dev_pivots, sizeof(int) * 4));

  // Transfer and run kernel.
  cudaCheck(cudaMemcpyAsync(dev_lu, lu, sizeof(double) * 16, cudaMemcpyHostToDevice, hstream));
  cudaCheck(cudaMemcpyAsync(dev_lus, lus, sizeof(double**), cudaMemcpyHostToDevice, hstream));
  cudaCheck(cudaMemcpyAsync(dev_pivots, pivots, sizeof(int) * 4, cudaMemcpyHostToDevice, hstream));

  // The types of the pointers passed here matter
  cuBLAS_LU::computeLogDet_batched(cuda_handles->hstream, n, lda, dev_lus, dev_pivots, dev_log_values, batch_size);

  cudaCheck(
      cudaMemcpyAsync(log_values, dev_log_values, sizeof(std::complex<double>) * 1, cudaMemcpyDeviceToHost, hstream));
  cudaCheck(cudaStreamSynchronize(hstream));


  CHECK(*log_values == ComplexApprox(std::complex<double>{5.267858159063328, 6.283185307179586}));

  // Free memory
  cudaCheck(cudaFree(dev_pivots));
  cudaCheck(cudaFreeHost(vp_pivots));
  cudaCheck(cudaFree(dev_log_values));
  cudaCheck(cudaFreeHost(log_values));
  cudaCheck(cudaFree(dev_lus));
  cudaCheck(cudaFreeHost(vp_lus));
  cudaCheck(cudaFree(dev_lu));
  cudaCheck(cudaFreeHost(lu));
}

TEST_CASE("cuBLAS_LU::computeLogDet_complex", "[wavefunction][CUDA]")
{
  auto cuda_handles = std::make_unique<testing::CUDAHandles>();
  int n             = 4;
  int lda           = 4;
  int batch_size = 1;
  auto& hstream     = cuda_handles->hstream;

  void* vp_lu = nullptr;
  cudaErrorCheck(cudaMallocHost(&vp_lu, sizeof(double) * 32), "cudaMallocHost failed");
  // clang-format off
  double* lu   = new (vp_lu) double[32]{8.0,                   0.5,
                                        0.8793774319066148,    0.07003891050583658,
                                        0.24980544747081712,   -0.0031128404669260694,
                                        0.6233463035019455,    -0.026459143968871595,
                                        2.0,                   0.1,
                                        6.248249027237354,     0.2719844357976654,
                                        0.7194170575332381,    -0.01831314754114669,
                                        0.1212375092639108,    0.02522449751055713,
                                        6.0,                   -0.2,
                                        0.7097276264591441,    -0.4443579766536965,
                                        4.999337315778741,     0.6013141870887196,
                                        0.26158183940834034,   0.23245112532996867,
                                        4.0,                   -0.6,
                                        4.440466926070039,     -1.7525291828793774,
                                        0.840192589866152,     1.5044529443071093,
                                        1.0698651110730424,    -0.10853319738453365};
  // clang-format on
  std::complex<double>* dev_lu;
  cudaErrorCheck(cudaMalloc((void**)&dev_lu, sizeof(double) * 32), "cudaMalloc failed");
  void* vp_lus;
  cudaErrorCheck(cudaMallocHost(&vp_lus, sizeof(double**)), "cudaMallocHost failed");
  std::complex<double>** lus = new (vp_lus) std::complex<double>* [1] { nullptr };
  lus[0]                     = dev_lu;
  std::complex<double>** dev_lus;
  cudaErrorCheck(cudaMalloc((void**)&dev_lus, sizeof(double**)), "cudaMallocHost failed");

  std::complex<double>* log_values;
  cudaErrorCheck(cudaMallocHost((void**)&log_values, sizeof(std::complex<double>) * 1), "cudaMallocHost failed");
  std::complex<double>* dev_log_values;

  cudaErrorCheck(cudaMalloc((void**)&dev_log_values, sizeof(std::complex<double>) * 1), "cudaMalloc failed");

  void* vp_pivots;
  cudaErrorCheck(cudaMallocHost((void**)&vp_pivots, sizeof(int) * 4), "cudaMallocHost failed");
  int* pivots = new (vp_pivots) int[4]{3, 4, 3, 4};
  int* dev_pivots;
  cudaErrorCheck(cudaMalloc((void**)&dev_pivots, sizeof(int) * 4), "cudaMalloc failed");

  cudaErrorCheck(cudaMemcpyAsync(dev_lu, lu, sizeof(double) * 32, cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying log_values to device");

  cudaErrorCheck(cudaMemcpyAsync(dev_lus, lus, sizeof(double*), cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying lus to device");
  cudaErrorCheck(cudaMemcpyAsync(dev_pivots, pivots, sizeof(int) * 4, cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying log_values to device");

  cuBLAS_LU::computeLogDet_batched(cuda_handles->hstream, n, lda, dev_lus, dev_pivots, dev_log_values, batch_size);

  cudaErrorCheck(cudaMemcpyAsync(log_values, dev_log_values, sizeof(std::complex<double>) * 1, cudaMemcpyDeviceToHost,
                                 hstream),
                 "cudaMemcpyAsync failed copying log_values from device");
  cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!");

  CHECK(log_values[0] == ComplexApprox(std::complex<double>{5.603777579195571, -6.1586603331188225}));

  // There are no destructors to call since all the inplace news were on POD arrays.
  cudaCheck(cudaFree(dev_pivots));
  cudaCheck(cudaFreeHost(vp_pivots));
  cudaCheck(cudaFree(dev_log_values));
  cudaCheck(cudaFreeHost(log_values));
  cudaCheck(cudaFree(dev_lus));
  cudaCheck(cudaFreeHost(vp_lus));
  cudaCheck(cudaFree(dev_lu));
  cudaCheck(cudaFreeHost(lu));
}

/** while this working is a good test, in production code its likely we want to
 *  widen the matrix M to double and thereby the LU matrix as well.
 */
TEST_CASE("cuBLAS_LU::computeLogDet_float", "[wavefunction][CUDA]")
{
  auto cuda_handles = std::make_unique<testing::CUDAHandles>();
  int n             = 4;
  int lda           = 4;
  int batch_size    = 1;
  auto& hstream     = cuda_handles->hstream;

  // The calls in cuBLAS_LU still care about type
  // being very careful and clear here.
  void* vp_lu = nullptr;
  cudaCheck(cudaMallocHost(&vp_lu, sizeof(float) * 16));
  // This is from numpy which did the LU factorization using float128
  float* lu = new (vp_lu) float[16]{7., 0.28571429, 0.71428571, 0.71428571,  5., 3.57142857, 0.12, -0.44,
                                      6., 6.28571429, -1.04,      -0.46153846, 6., 5.28571429, 3.08, 7.46153846};  
  float* dev_lu;
  cudaCheck(cudaMalloc((void**)&dev_lu, sizeof(float) * 16));
  void* vp_lus;
  cudaCheck(cudaMallocHost(&vp_lus, sizeof(float**)));
  float** lus = new (vp_lus) float* [batch_size] { nullptr };
  lus[0]       = dev_lu;
  float** dev_lus;
  cudaCheck(cudaMalloc((void**)&dev_lus, sizeof(float**)));

  std::complex<double>* log_values;
  cudaCheck(cudaMallocHost((void**)&log_values, sizeof(std::complex<double>) * 1));
  std::complex<double>* dev_log_values;
  cudaCheck(cudaMalloc((void**)&dev_log_values, sizeof(std::complex<double>) * 1));
  CUDAfill_n(dev_log_values, batch_size, {0,0});

  void* vp_pivots;
  cudaCheck(cudaMallocHost(&vp_pivots, sizeof(int) * 4));
  int* pivots = new (vp_pivots) int[4]{3, 3, 4, 4};
  int* dev_pivots;
  cudaCheck(cudaMalloc((void**)&dev_pivots, sizeof(int) * 4));

  // Transfer and run kernel.
  cudaCheck(cudaMemcpyAsync(dev_lu, lu, sizeof(float) * 16, cudaMemcpyHostToDevice, hstream));
  cudaCheck(cudaMemcpyAsync(dev_lus, lus, sizeof(float**), cudaMemcpyHostToDevice, hstream));
  cudaCheck(cudaMemcpyAsync(dev_pivots, pivots, sizeof(int) * 4, cudaMemcpyHostToDevice, hstream));

  // The types of the pointers passed here matter
  cuBLAS_LU::computeLogDet_batched(hstream, n, lda, dev_lus, dev_pivots, dev_log_values, batch_size);

  cudaCheck(
      cudaMemcpyAsync(log_values, dev_log_values, sizeof(std::complex<double>) * 1, cudaMemcpyDeviceToHost, hstream));
  cudaCheck(cudaStreamSynchronize(hstream));


  CHECK(*log_values == ComplexApprox(std::complex<double>{5.267858159063328, 6.283185307179586}));

  // Free memory
  cudaCheck(cudaFree(dev_pivots));
  cudaCheck(cudaFreeHost(vp_pivots));
  cudaCheck(cudaFree(dev_log_values));
  cudaCheck(cudaFreeHost(log_values));
  cudaCheck(cudaFree(dev_lus));
  cudaCheck(cudaFreeHost(vp_lus));
  cudaCheck(cudaFree(dev_lu));
  cudaCheck(cudaFreeHost(lu));
}

TEST_CASE("cuBLAS_LU::computeLogDet(batch=2)", "[wavefunction][CUDA]")
{
  auto cuda_handles = std::make_unique<testing::CUDAHandles>();
  int n             = 4;
  int lda           = 4;
  auto& hstream     = cuda_handles->hstream;

  void* vp_lu = nullptr;
  // clang-format off
  cudaErrorCheck(cudaMallocHost(&vp_lu, sizeof(double) * 64), "cudaMallocHost failed");
  double* lu   = new (vp_lu) double[32]{8.0,                   0.5,
                                        0.8793774319066148,    0.07003891050583658,
                                        0.24980544747081712,   -0.0031128404669260694,
                                        0.6233463035019455,    -0.026459143968871595,
                                        2.0,                   0.1,
                                        6.248249027237354,     0.2719844357976654,
                                        0.7194170575332381,    -0.01831314754114669,
                                        0.1212375092639108,    0.02522449751055713,
                                        6.0,                   -0.2,
                                        0.7097276264591441,    -0.4443579766536965,
                                        4.999337315778741,     0.6013141870887196,
                                        0.26158183940834034,   0.23245112532996867,
                                        4.0,                   -0.6,
                                        4.440466926070039,     -1.7525291828793774,
                                        0.840192589866152,     1.5044529443071093,
                                        1.0698651110730424,    -0.10853319738453365};
  void* vp_lu2 = reinterpret_cast<void*>(lu + 32);
  double* lu2  = new (vp_lu2) double[32]{8.0, 0.5,
					 0.8793774319066148, 0.07003891050583658,
					 0.49883268482490273, -0.01867704280155642,
					 0.24980544747081712, -0.0031128404669260694,
					 2.0, 0.1,
					 6.248249027237354, 0.2719844357976654,
					 0.800088933543564, -0.004823898651572499,
					 0.2401906003014191, 0.0025474386841018853,
					 3.0, -0.2,
					 3.3478599221789884, -0.23424124513618677,
					 0.8297816353227319, 1.3593612303468308,
					 0.6377685195602139, -0.6747848919351336,
					 4.0, -0.6,
					 4.440466926070039, -1.7525291828793774,
					 -1.5284389377713894, 1.6976073494521235,
					 2.7608934839023482, -1.542084179899335};
  // clang-format off
  std::complex<double>* dev_lu;
  cudaErrorCheck(cudaMalloc((void**)&dev_lu, sizeof(double) * 64), "cudaMalloc failed");

  void* vp_lus;
  cudaErrorCheck(cudaMallocHost(&vp_lus, sizeof(double**) * 2), "cudaMallocHost failed");
  std::complex<double>** lus = new (vp_lus) std::complex<double>* [2] { nullptr, nullptr };
  lus[0]                     = dev_lu;
  lus[1]                     = dev_lu + 16;
  std::complex<double>** dev_lus;
  cudaErrorCheck(cudaMalloc((void**)&dev_lus, sizeof(double*) * 2), "cudaMallocHost failed");

  void* vp_log_values;
  cudaErrorCheck(cudaMallocHost(&vp_log_values, sizeof(std::complex<double>) * 2), "cudaMallocHost failed");
  // For values we expect zeroed as a side effect of a call we should poison them to test.
  std::complex<double>* log_values =
      reinterpret_cast<std::complex<double>*>(new (vp_log_values) double[4]{1.0, 1.0, 1.0, 1.0});

  std::complex<double>* dev_log_values;
  cudaErrorCheck(cudaMalloc((void**)&dev_log_values, sizeof(std::complex<double>) * 2), "cudaMalloc failed");

  void* vp_pivots;
  cudaErrorCheck(cudaMallocHost((void**)&vp_pivots, sizeof(int) * 8), "cudaMallocHost failed");
  int* pivots      = new (vp_pivots) int[4]{3, 4, 3, 4};
  void* vp_pivots2 = (void*)(pivots + 4);
  int* pivots2     = new (vp_pivots2) int[4]{3, 4, 4, 4};

  int* dev_pivots;
  cudaErrorCheck(cudaMalloc((void**)&dev_pivots, sizeof(int) * 8), "cudaMalloc failed");

  cudaErrorCheck(cudaMemcpyAsync(dev_lu, vp_lu, sizeof(double) * 64, cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying log_values to device");

  cudaErrorCheck(cudaMemcpyAsync(dev_lus, lus, sizeof(std::complex<double>*) * 2, cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying log_values to device");

  cudaErrorCheck(cudaMemcpyAsync(dev_pivots, vp_pivots, sizeof(int) * 8, cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying log_values to device");

  int batch_size = 2;

  cuBLAS_LU::computeLogDet_batched(cuda_handles->hstream, n, lda, dev_lus, dev_pivots, dev_log_values, batch_size);
  cudaErrorCheck(cudaMemcpyAsync(log_values, dev_log_values, sizeof(std::complex<double>) * 2, cudaMemcpyDeviceToHost,
                                 hstream),
                 "cudaMemcpyAsync failed copying log_values from device");
  cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!");

  CHECK(log_values[0] == ComplexApprox(std::complex<double>{ 5.603777579195571, -6.1586603331188225 }));
  CHECK(log_values[1] == ComplexApprox(std::complex<double>{ 5.531331998282581, -8.805487075984523  }));
  cudaCheck(cudaFree(dev_pivots));
  cudaCheck(cudaFreeHost(vp_pivots));
  cudaCheck(cudaFree(dev_log_values));
  cudaCheck(cudaFreeHost(log_values));
  cudaCheck(cudaFree(dev_lus));
  cudaCheck(cudaFreeHost(vp_lus));
  cudaCheck(cudaFree(dev_lu));
  cudaCheck(cudaFreeHost(lu));

}


TEST_CASE("cuBLAS_LU::getrf_batched_complex", "[wavefunction][CUDA]")
{
  auto cuda_handles = std::make_unique<testing::CUDAHandles>();
  int n             = 4;
  int lda           = 4;
  auto& hstream     = cuda_handles->hstream;

  void* vpM = nullptr;
  cudaErrorCheck(cudaMallocHost(&vpM, sizeof(double) * 32), "cudaMallocHost failed");

  // this should be read as a column major matrix by cuBLAS
  double* M = new (vpM) double[32]{
      2.0, 0.1, 5.0, 0.1, 8.0, 0.5,  7.0, 1.0,  5.0, 0.1, 2.0, 0.2,  2.0, 0.1,  8.0, 0.5,
      7.0, 0.2, 5.0, 1.0, 6.0, -0.2, 6.0, -0.2, 5.0, 0.0, 4.0, -0.1, 4.0, -0.6, 8.0, -2.0,
  };

  double* devM;
  cudaErrorCheck(cudaMalloc((void**)&devM, sizeof(double) * 32), "cudaMalloc failed");
  std::complex<double>** Ms;
  cudaErrorCheck(cudaMallocHost((void**)&Ms, sizeof(double*)), "cudaMallocHost failed");
  Ms = reinterpret_cast<std::complex<double>**>(&devM);
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

  cuBLAS_LU::computeGetrf_batched(cuda_handles->h_cublas, cuda_handles->hstream,  n, lda, devMs, dev_pivots, infos, dev_infos, batch_size);

  cudaErrorCheck(cudaMemcpyAsync(M, devM, sizeof(double) * 32, cudaMemcpyDeviceToHost, hstream),
                 "cudaMemcpyAsync failed copying invM from device");
  cudaErrorCheck(cudaMemcpyAsync(pivots, dev_pivots, sizeof(int) * 4, cudaMemcpyDeviceToHost, hstream),
                 "cudaMemcpyAsync failed copying pivots from device");

  cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!");

  int real_pivot[4]{3, 4, 3, 4};

  auto checkArray = [](auto* A, auto* B, int n) {
    for (int i = 0; i < n; ++i)
    {
      CHECK(A[i] == B[i]);
    }
  };
  checkArray(real_pivot, pivots, 4);

  double lu[32]{8.0,
                0.5,
                0.8793774319066148,
                0.07003891050583658,
                0.24980544747081712,
                -0.0031128404669260694,
                0.6233463035019455,
                -0.026459143968871595,
                2.0,
                0.1,
                6.248249027237354,
                0.2719844357976654,
                0.7194170575332381,
                -0.01831314754114669,
                0.1212375092639108,
                0.02522449751055713,
                6.0,
                -0.2,
                0.7097276264591441,
                -0.4443579766536965,
                4.999337315778741,
                0.6013141870887196,
                0.26158183940834034,
                0.23245112532996867,
                4.0,
                -0.6,
                4.440466926070039,
                -1.7525291828793774,
                0.840192589866152,
                1.5044529443071093,
                1.0698651110730424,
                -0.10853319738453365};

  // This could actually be any container that supported the concept of
  // access via operator()(i, j) and had <T, ALLOCT> template signature
  testing::MatrixAccessor<std::complex<double>> lu_mat(reinterpret_cast<std::complex<double>*>(lu), 4, 4);
  testing::MatrixAccessor<std::complex<double>> M_mat(reinterpret_cast<std::complex<double>*>(M), 4, 4);
  auto check_matrix_result = checkMatrix(lu_mat, M_mat);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
}

TEST_CASE("cuBLAS_LU::getrf_batched(batch=2)", "[wavefunction][CUDA]")
{
  auto cuda_handles = std::make_unique<testing::CUDAHandles>();
  int n             = 4;
  int lda           = 4;
  auto& hstream     = cuda_handles->hstream;

  int batch_size = 2;

  // This is the array of double** with the dev pointer address ie (double*[])
  double** Ms;
  cudaErrorCheck(cudaMallocHost((void**)&Ms, sizeof(double*) * 2), "cudaMallocHost failed");
  double** devMs;
  cudaErrorCheck(cudaMalloc((void**)&devMs, sizeof(double*) * 2), "cudaMalloc failed");

  // These are the M matrices
  void* vpM = nullptr;
  cudaErrorCheck(cudaMallocHost(&vpM, sizeof(double) * 16), "cudaMallocHost failed");
  double* M = new (vpM) double[16]{2, 5, 7, 5, 5, 2, 5, 4, 8, 2, 6, 4, 7, 8, 6, 8};
  double* devM;
  cudaErrorCheck(cudaMalloc((void**)&devM, sizeof(double) * 16), "cudaMalloc failed");
  void* vpM2 = nullptr;
  cudaErrorCheck(cudaMallocHost(&vpM2, sizeof(double) * 16), "cudaMallocHost failed");
  double* M2 = new (vpM2) double[16]{6, 5, 7, 5, 2, 2, 5, 4, 8, 2, 6, 4, 3, 8, 6, 8};
  double* devM2;
  cudaErrorCheck(cudaMalloc((void**)&devM2, sizeof(double) * 16), "cudaMalloc failed");

  // Now we get their pointers on the device and store them in the host double*[]
  *Ms       = devM;
  *(Ms + 1) = devM2;

  // Now we allocate the pivot and info arrays for both matrices
  void* vp_pivots;
  cudaErrorCheck(cudaMallocHost((void**)&vp_pivots, sizeof(int) * 8), "cudaMallocHost failed");
  int* pivots = new (vp_pivots) int[8]{-1, -1, -1, -1, -1, -1, -1, -1};
  int* dev_pivots;
  cudaErrorCheck(cudaMalloc((void**)&dev_pivots, sizeof(int) * 8), "cudaMalloc failed");

  void* vp_infos;
  cudaErrorCheck(cudaMallocHost((void**)&vp_infos, sizeof(int) * 8), "cudaMallocHost failed");
  int* infos = new (vp_infos) int[8]{1, 1, 1, 1, 1, 1, 1, 1};
  int* dev_infos;
  cudaErrorCheck(cudaMalloc((void**)&dev_infos, sizeof(int) * 8), "cudaMalloc failed");

  //Now copy the Ms
  cudaErrorCheck(cudaMemcpyAsync(devM, M, sizeof(double) * 16, cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying M to device");
  cudaErrorCheck(cudaMemcpyAsync(devM2, M2, sizeof(double) * 16, cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying M to device");
  // copy the pointer array
  cudaErrorCheck(cudaMemcpyAsync(devMs, Ms, sizeof(double*) * 2, cudaMemcpyHostToDevice, hstream),
                 "cudaMemcpyAsync failed copying Ms to device");

  cuBLAS_LU::computeGetrf_batched(cuda_handles->h_cublas, cuda_handles->hstream, n, lda, devMs, dev_pivots, infos, dev_infos, batch_size);

  // copy back the Ms, infos, pivots
  cudaErrorCheck(cudaMemcpyAsync(M, devM, sizeof(double) * 16, cudaMemcpyDeviceToHost, hstream),
                 "cudaMemcpyAsync failed copying invM from device");
  cudaErrorCheck(cudaMemcpyAsync(M2, devM2, sizeof(double) * 16, cudaMemcpyDeviceToHost, hstream),
                 "cudaMemcpyAsync failed copying invM from device");
  cudaErrorCheck(cudaMemcpyAsync(pivots, dev_pivots, sizeof(int) * 8, cudaMemcpyDeviceToHost, hstream),
                 "cudaMemcpyAsync failed copying pivots from device");

  cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!");

  // clang-format off
  std::vector<double> lu{7.,                    0.28571429,
                         0.71428571,            0.71428571,
                         5.,                    3.57142857,
                         0.12,                  -0.44,
                         6.,                    6.28571429,
                         -1.04,                 -0.46153846,
                         6.,                    5.28571429,
                         3.08,                  7.46153846};

  std::vector<double> lu2{7.0,                  0.8571428571428571,
                          0.7142857142857142,   0.7142857142857142,
                          5.0,                  -2.2857142857142856,
                          0.6874999999999998,   -0.18750000000000022,
                          6.0,                  2.8571428571428577,
                          -4.249999999999999,   -0.05882352941176502,
                          6.0,                  -2.1428571428571423,
                          5.1875,               3.617647058823531};
  // clang-format on
  int real_pivot[8]{3, 3, 4, 4, 3, 3, 3, 4};

  auto checkArray = [](auto* A, auto* B, int n) {
    for (int i = 0; i < n; ++i)
    {
      std::cout << A[i] << ":" << B[i] << ", ";
      CHECK(A[i] == Approx(B[i]));
    }
  };
  std::cout << '\n';

  testing::MatrixAccessor<double> M_mat(M, 4, 4);
  testing::MatrixAccessor<double> lu_mat(lu.data(), 4, 4);
  testing::MatrixAccessor<double> M2_mat(M2, 4, 4);
  testing::MatrixAccessor<double> lu2_mat(lu2.data(), 4, 4);

  checkArray(real_pivot, pivots, 8);
  auto check_matrix_result = checkMatrix(lu_mat, M_mat);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
  check_matrix_result = checkMatrix(lu2_mat, M2_mat);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
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
  void* vp_Ms;
  cudaErrorCheck(cudaMallocHost(&vp_Ms, sizeof(double*)), "cudaMallocHost failed");
  double** Ms = new (vp_Ms) double*[1] {nullptr};
  Ms[0] = devM;
  double** devMs;
  cudaErrorCheck(cudaMalloc((void**)&devMs, sizeof(double*)), "cudaMalloc failed");

  void* vp_invM = nullptr;
  cudaErrorCheck(cudaMallocHost(&vp_invM, sizeof(double) * 16), "cudaMallocHost failed");
  double* invM = new (vp_invM) double[16]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  double* dev_invM;
  cudaErrorCheck(cudaMalloc((void**)&dev_invM, sizeof(double) * 16), "cudaMalloc failed");

  void* vp_invMs;
  cudaErrorCheck(cudaMallocHost(&vp_invMs, sizeof(double*)), "cudaMallocHost failed");
  double** invMs = new (vp_invMs) double*[1] {nullptr};
  invMs[0] = dev_invM;
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
      CHECK(A[i] == Approx(B[i]));
    }
  };

  cudaCheck(cudaFree(devM));
  cudaCheck(cudaFreeHost(vpM));
  cudaCheck(cudaFree(devMs));
  cudaCheck(cudaFreeHost(vp_Ms));
  cudaCheck(cudaFree(dev_invM));
  cudaCheck(cudaFreeHost(vp_invM));
  cudaCheck(cudaFree(dev_invMs));
  cudaCheck(cudaFreeHost(vp_invMs));
  cudaCheck(cudaFree(dev_pivots));
  cudaCheck(cudaFreeHost(vp_pivots));
  cudaCheck(cudaFree(dev_infos));
  cudaCheck(cudaFreeHost(vp_infos));
  
  testing::MatrixAccessor<double> invA_mat(invA, 4, 4);
  testing::MatrixAccessor<double> invM_mat(invM, 4, 4);

  auto check_matrix_result = checkMatrix(invA_mat, invM_mat);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
}

} // namespace qmcplusplus
