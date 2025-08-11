//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
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
template<PlatformKind PL, typename T>
void test_one_gemm(const int M, const int N, const int K, const char transa, const char transb)
{
  const int a0 = transa == 'T' ? M : K;
  const int a1 = transa == 'T' ? K : M;

  const int b0 = transb == 'T' ? K : N;
  const int b1 = transb == 'T' ? N : K;

  using vec_t = Vector<T, PinnedDualAllocator<T>>;
  using mat_t = Matrix<T, PinnedDualAllocator<T>>;

  mat_t A(a0, a1); // Input matrix
  mat_t B(b0, b1); // Input matrix
  mat_t C(N, M);   // Result matrix AccelBLAS
  mat_t D(N, M);   // Result matrix BLAS

  // Fill data
  for (int j = 0; j < a0; j++)
    for (int i = 0; i < a1; i++)
      A[j][i] = i * 3 + j * 4;

  for (int j = 0; j < b0; j++)
    for (int i = 0; i < b1; i++)
      B[j][i] = i * 4 + j * 5;

  A.updateTo();
  B.updateTo();
  C.updateTo();

  T alpha(1);
  T beta(0);

  // U[X,Y] denotes a row-major matrix U with X rows and Y cols
  // element U(i,j) is located at: U.data() + sizeof(U_type) * (i*ldU + j)
  //
  // A,B,C,D are treated as row-major matrices, but the arguments to gemm are treated as col-major
  // so the call below to AccelBLAS::gemm is equivalent to one of the following (with row-major matrices)
  // transa/transb == 'N'/'N':   C[N,M] = A[K,M] * B[N,K]; C = B * A
  // transa/transb == 'N'/'T':   C[N,M] = A[K,M] * B[K,N]; C = B^t * A
  // transa/transb == 'T'/'N':   C[N,M] = A[M,K] * B[N,K]; C = B * A^t
  // transa/transb == 'T'/'T':   C[N,M] = A[M,K] * B[K,N]; C = B^t * A^t

  compute::Queue<PL> queue;
  compute::BLASHandle<PL> h_blas(queue);
  compute::BLAS::gemm(h_blas, transa, transb, M, N, K, alpha, A.device_data(), a1, B.device_data(), b1, beta,
                      C.device_data(), M);
  queue.sync();

  C.updateFrom();

  BLAS::gemm(transa, transb, M, N, K, alpha, A.data(), a1, B.data(), b1, beta, D.data(), M);

  for (int j = 0; j < N; j++)
    for (int i = 0; i < M; i++)
    {
      CHECK(std::real(C[j][i]) == Approx(std::real(D[j][i])));
      CHECK(std::imag(C[j][i]) == Approx(std::imag(D[j][i])));
    }

  mat_t A2(a0, a1); // Input matrix
  mat_t B2(b0, b1); // Input matrix
  mat_t C2(N, M);   // Result matrix AccelBLAS
  mat_t D2(N, M);   // Result matrix BLAS

  // Fill data
  for (int j = 0; j < a0; j++)
    for (int i = 0; i < a1; i++)
      A2[j][i] = j * 3 + i * 4;

  for (int j = 0; j < b0; j++)
    for (int i = 0; i < b1; i++)
      B2[j][i] = j * 4 + i * 5;

  A2.updateTo();
  B2.updateTo();

  Vector<const T*, PinnedDualAllocator<const T*>> Aarr(2), Barr(2);
  Vector<T*, PinnedDualAllocator<T*>> Carr(2);

  Aarr[0] = A2.device_data();
  Aarr[1] = A.device_data();
  Barr[0] = B2.device_data();
  Barr[1] = B.device_data();

  Carr[0] = C.device_data();
  Carr[1] = C2.device_data();

  Aarr.updateTo();
  Barr.updateTo();
  Carr.updateTo();

  compute::BLAS::gemm_batched(h_blas, transa, transb, M, N, K, alpha, Aarr.device_data(), a1, Barr.device_data(), b1,
                              beta, Carr.device_data(), M, 2);
  queue.sync();

  C.updateFrom();
  C2.updateFrom();

  BLAS::gemm(transa, transb, M, N, K, alpha, A2.data(), a1, B2.data(), b1, beta, D2.data(), M);

  for (int j = 0; j < N; j++)
    for (int i = 0; i < M; i++)
    {
      CHECK(std::real(C2[j][i]) == Approx(std::real(D[j][i])));
      CHECK(std::imag(C2[j][i]) == Approx(std::imag(D[j][i])));
    }

  for (int j = 0; j < N; j++)
    for (int i = 0; i < M; i++)
    {
      CHECK(std::real(C[j][i]) == Approx(std::real(D2[j][i])));
      CHECK(std::imag(C[j][i]) == Approx(std::imag(D2[j][i])));
    }
}

template<PlatformKind PL>
void test_gemm_cases()
{
  const int M = 29;
  const int N = 31;
  const int K = 23;

  // Non-batched test
  std::cout << "Testing NN gemm" << std::endl;
  test_one_gemm<PL, float>(M, N, K, 'N', 'N');
  test_one_gemm<PL, double>(M, N, K, 'N', 'N');
#if defined(QMC_COMPLEX)
  test_one_gemm<PL, std::complex<float>>(N, M, K, 'N', 'N');
  test_one_gemm<PL, std::complex<double>>(N, M, K, 'N', 'N');
#endif
  std::cout << "Testing NT gemm" << std::endl;
  test_one_gemm<PL, float>(M, N, K, 'N', 'T');
  test_one_gemm<PL, double>(M, N, K, 'N', 'T');
#if defined(QMC_COMPLEX)
  test_one_gemm<PL, std::complex<float>>(N, M, K, 'N', 'T');
  test_one_gemm<PL, std::complex<double>>(N, M, K, 'N', 'T');
#endif
  std::cout << "Testing TN gemm" << std::endl;
  test_one_gemm<PL, float>(M, N, K, 'T', 'N');
  test_one_gemm<PL, double>(M, N, K, 'T', 'N');
#if defined(QMC_COMPLEX)
  test_one_gemm<PL, std::complex<float>>(N, M, K, 'T', 'N');
  test_one_gemm<PL, std::complex<double>>(N, M, K, 'T', 'N');
#endif
  std::cout << "Testing TT gemm" << std::endl;
  test_one_gemm<PL, float>(M, N, K, 'T', 'T');
  test_one_gemm<PL, double>(M, N, K, 'T', 'T');
#if defined(QMC_COMPLEX)
  test_one_gemm<PL, std::complex<float>>(N, M, K, 'T', 'T');
  test_one_gemm<PL, std::complex<double>>(N, M, K, 'T', 'T');
#endif
}

template<PlatformKind PL, typename T>
void test_one_gemv(const int M_b, const int N_b, const char trans)
{
  const int M = trans == 'T' ? M_b : N_b;
  const int N = trans == 'T' ? N_b : M_b;

  using vec_t = Vector<T, PinnedDualAllocator<T>>;
  using mat_t = Matrix<T, PinnedDualAllocator<T>>;

  vec_t A(N);        // Input vector
  mat_t B(M_b, N_b); // Input matrix
  vec_t C(M);        // Result vector AccelBLAS
  vec_t D(M);        // Result vector BLAS

  // Fill data
  for (int i = 0; i < N; i++)
    A[i] = i;

  for (int j = 0; j < M_b; j++)
    for (int i = 0; i < N_b; i++)
      B[j][i] = i + j * 2;

  // Fill C and D with 0
  for (int i = 0; i < M; i++)
    C[i] = D[i] = T(0);

  A.updateTo();
  B.updateTo();
  C.updateTo();

  T alpha(0.5);
  T beta(0);

  compute::Queue<PL> queue;
  compute::BLASHandle<PL> h_blas(queue);
  // in Fortran, B[M][N] is viewed as B^T
  // when trans == 'T', the actual calculation is B * A[N] = C[M]
  // when trans == 'N', the actual calculation is B^T * A[M] = C[N]
  compute::BLAS::gemv(h_blas, trans, N_b, M_b, alpha, B.device_data(), N_b, A.device_data(), 1, beta, C.device_data(),
                      1);
  queue.sync();

  C.updateFrom();

  if (trans == 'T')
    BLAS::gemv_trans(M_b, N_b, B.data(), A.data(), D.data());
  else
    BLAS::gemv(M_b, N_b, B.data(), A.data(), D.data());

  for (int index = 0; index < M; index++)
    CHECK(C[index] == D[index] * alpha);

  vec_t A2(N);        // Input vector
  mat_t B2(M_b, N_b); // Input matrix
  vec_t C2(M);        // Result vector AccelBLAS
  vec_t D2(M);        // Result vector BLAS

  // Fill data
  for (int i = 0; i < N; i++)
    A2[i] = i + 1;

  for (int j = 0; j < M_b; j++)
    for (int i = 0; i < N_b; i++)
      B2[j][i] = i * 2 + j;

  // Fill C and D with 0
  for (int i = 0; i < M; i++)
    C2[i] = D2[i] = T(0);

  A2.updateTo();
  B2.updateTo();
  C2.updateTo();

  Vector<const T*, PinnedDualAllocator<const T*>> Aarr(2), Barr(2);
  Vector<T*, PinnedDualAllocator<T*>> Carr(2);
  Vector<T, PinnedDualAllocator<T>> alpha_arr(2), beta_arr(2);

  Aarr[0] = A2.device_data();
  Aarr[1] = A.device_data();
  Barr[0] = B2.device_data();
  Barr[1] = B.device_data();

  Carr[0] = C2.device_data();
  Carr[1] = C.device_data();

  alpha_arr[0] = 2;
  alpha_arr[1] = 0.5;
  beta_arr[0]  = 0;
  beta_arr[1]  = 1;

  Aarr.updateTo();
  Barr.updateTo();
  Carr.updateTo();
  alpha_arr.updateTo();
  beta_arr.updateTo();

  // in Fortran, B[M][N] is viewed as B^T
  // when trans == 'T', the actual calculation is B * A[N] = C[M]
  // when trans == 'N', the actual calculation is B^T * A[M] = C[N]
  compute::BLAS::gemv_batched(h_blas, trans, N_b, M_b, alpha_arr.device_data(), Barr.device_data(), N_b,
                              Aarr.device_data(), 1, beta_arr.device_data(), Carr.device_data(), 1, 2);
  queue.sync();

  C.updateFrom();
  C2.updateFrom();

  if (trans == 'T')
    BLAS::gemv_trans(M_b, N_b, B2.data(), A2.data(), D2.data());
  else
    BLAS::gemv(M_b, N_b, B2.data(), A2.data(), D2.data());

  for (int index = 0; index < M; index++)
  {
    CHECK(C[index] == D[index]);
    CHECK(C2[index] == D2[index] * alpha_arr[0]);
  }
}

template<PlatformKind PL>
void test_gemv_cases()
{
  const int M = 137;
  const int N = 79;

  std::cout << "Testing NOTRANS gemv" << std::endl;
  test_one_gemv<PL, float>(M, N, 'N');
  test_one_gemv<PL, double>(M, N, 'N');
#if defined(QMC_COMPLEX)
  test_one_gemv<PL, std::complex<float>>(N, M, 'N');
  test_one_gemv<PL, std::complex<double>>(N, M, 'N');
#endif
  std::cout << "Testing TRANS gemv" << std::endl;
  test_one_gemv<PL, float>(M, N, 'T');
  test_one_gemv<PL, double>(M, N, 'T');
#if defined(QMC_COMPLEX)
  test_one_gemv<PL, std::complex<float>>(N, M, 'T');
  test_one_gemv<PL, std::complex<double>>(N, M, 'T');
#endif
}

template<PlatformKind PL, typename T>
void test_one_ger(const int M, const int N)
{
  using vec_t = Vector<T, PinnedDualAllocator<T>>;
  using mat_t = Matrix<T, PinnedDualAllocator<T>>;

  mat_t Ah(M, N); // Input matrix
  mat_t Ad(M, N); // Input matrix
  vec_t x(M);     // Input vector
  vec_t y(N);     // Input vector

  // Fill data
  for (int i = 0; i < M; i++)
    x[i] = i;
  for (int i = 0; i < N; i++)
    y[i] = N - i;

  for (int j = 0; j < M; j++)
    for (int i = 0; i < N; i++)
    {
      Ah[j][i] = i + j * 2;
      Ad[j][i] = i + j * 2;
    }

  Ad.updateTo();
  x.updateTo();
  y.updateTo();

  T alpha(0.5);

  compute::Queue<PL> queue;
  compute::BLASHandle<PL> h_blas(queue);
  // in Fortran, B[M][N] is viewed as B^T
  compute::BLAS::ger(h_blas, M, N, alpha, x.device_data(), 1, y.device_data(), 1, Ad.device_data(), M);
  queue.sync();
  Ad.updateFrom();

  BLAS::ger(M, N, alpha, x.data(), 1, y.data(), 1, Ah.data(), M);

  for (int j = 0; j < M; j++)
    for (int i = 0; i < N; i++)
      CHECK(Ah[j][i] == Ad[j][i]);

  mat_t Ah2(M, N); // Input matrix
  mat_t Ad2(M, N); // Input matrix
  vec_t x2(M);     // Input vector
  vec_t y2(N);     // Input vector

  // Fill data
  for (int i = 0; i < M; i++)
    x2[i] = i - 1;
  for (int i = 0; i < N; i++)
    y2[i] = N + i;

  for (int j = 0; j < M; j++)
    for (int i = 0; i < N; i++)
    {
      Ah2[j][i] = j + i * 2;
      Ad2[j][i] = j + i * 2;
    }

  Ad2.updateTo();
  x2.updateTo();
  y2.updateTo();

  Vector<T*, PinnedDualAllocator<T*>> Aarr(2);
  Vector<const T*, PinnedDualAllocator<const T*>> Xarr(2), Yarr(2);
  Vector<T, PinnedDualAllocator<T>> alpha_arr(2);

  Aarr[0] = Ad2.device_data();
  Aarr[1] = Ad.device_data();
  Xarr[0] = x2.device_data();
  Xarr[1] = x.device_data();
  Yarr[0] = y2.device_data();
  Yarr[1] = y.device_data();

  alpha_arr[0] = 2;
  alpha_arr[1] = 0.5;

  Aarr.updateTo();
  Xarr.updateTo();
  Yarr.updateTo();
  alpha_arr.updateTo();

  // in Fortran, B[M][N] is viewed as B^T
  compute::BLAS::ger_batched(h_blas, M, N, alpha_arr.device_data(), Xarr.device_data(), 1, Yarr.device_data(), 1,
                             Aarr.device_data(), M, 2);
  queue.sync();
  Ad.updateFrom();
  Ad2.updateFrom();

  BLAS::ger(M, N, alpha_arr[1], x.data(), 1, y.data(), 1, Ah.data(), M);
  BLAS::ger(M, N, alpha_arr[0], x2.data(), 1, y2.data(), 1, Ah2.data(), M);

  for (int j = 0; j < M; j++)
    for (int i = 0; i < N; i++)
    {
      CHECK(Ah[j][i] == Ad[j][i]);
      CHECK(Ah2[j][i] == Ad2[j][i]);
    }
}

template<PlatformKind PL>
void test_ger_cases()
{
  const int M = 137;
  const int N = 79;

  // Batched Test
  std::cout << "Testing ger_batched" << std::endl;
  test_one_ger<PL, float>(M, N);
  test_one_ger<PL, double>(M, N);
#if defined(QMC_COMPLEX)
  test_one_ger<PL, std::complex<float>>(N, M);
  test_one_ger<PL, std::complex<double>>(N, M);
#endif
}

TEST_CASE("AccelBLAS", "[BLAS]")
{
  SECTION("gemm")
  {
#if defined(ENABLE_CUDA)
    std::cout << "Testing gemm<PlatformKind::CUDA>" << std::endl;
    test_gemm_cases<PlatformKind::CUDA>();
#endif
#if defined(ENABLE_SYCL)
    std::cout << "Testing gemm<PlatformKind::SYCL>" << std::endl;
    test_gemm_cases<PlatformKind::SYCL>();
#endif
#if defined(ENABLE_OFFLOAD)
    std::cout << "Testing gemm<PlatformKind::OMPTARGET>" << std::endl;
    test_gemm_cases<PlatformKind::OMPTARGET>();
#endif
  }

  SECTION("gemv")
  {
#if defined(ENABLE_CUDA)
    std::cout << "Testing gemm<PlatformKind::CUDA>" << std::endl;
    test_gemv_cases<PlatformKind::CUDA>();
#endif
#if defined(ENABLE_SYCL)
    std::cout << "Testing gemm<PlatformKind::SYCL>" << std::endl;
    test_gemv_cases<PlatformKind::SYCL>();
#endif
#if defined(ENABLE_OFFLOAD)
    std::cout << "Testing gemm<PlatformKind::OMPTARGET>" << std::endl;
    test_gemv_cases<PlatformKind::OMPTARGET>();
#endif
  }

  SECTION("ger")
  {
#if defined(ENABLE_CUDA)
    std::cout << "Testing ger<PlatformKind::CUDA>" << std::endl;
    test_ger_cases<PlatformKind::CUDA>();
#endif
#if defined(ENABLE_SYCL)
    std::cout << "Testing ger<PlatformKind::SYCL>" << std::endl;
    test_ger_cases<PlatformKind::SYCL>();
#endif
#if defined(ENABLE_OFFLOAD)
    std::cout << "Testing ger<PlatformKind::OMPTARGET>" << std::endl;
    test_ger_cases<PlatformKind::OMPTARGET>();
#endif
  }
}

template<PlatformKind PL>
void benchmark_gemm_cases()
{
  const int M_b = 29;
  const int N_b = 31;
  const int K_b = 23;
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
  const int M_b = 137;
  const int N_b = 79;
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
  const int M_b = 137;
  const int N_b = 79;
  const int BATCH = 64;

  // Batched Test
  benchmark::perf_ger_batched<PL, float>(M_b, N_b, BATCH);
  benchmark::perf_ger_batched<PL, double>(M_b, N_b, BATCH);
#if defined(QMC_COMPLEX)
  benchmark::perf_ger_batched<PL, std::complex<float>>(M_b, N_b, BATCH);
  benchmark::perf_ger_batched<PL, std::complex<double>>(M_b, N_b, BATCH);
#endif
}


TEST_CASE("AccelBLAS batched benchmark", "[BLAS][benchmark]") {

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
