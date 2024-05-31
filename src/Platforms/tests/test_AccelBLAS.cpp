//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include <DualAllocatorAliases.hpp>
#include <AccelBLASAliases.hpp>
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include <CPU/BLAS.hpp>

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
  mat_t C(N, M);   // Result matrix ompBLAS
  mat_t D(N, M);   // Result matrix BLAS

  // Fill data
  for (int j = 0; j < a0; j++)
    for (int i = 0; i < a1; i++)
      A[j][i] = i * 3 + j * 4;

  for (int j = 0; j < b0; j++)
    for (int i = 0; i < b1; i++)
      B[j][i] = i * 4 + j * 5;

  // Fill C and D with 0
  for (int j = 0; j < N; j++)
    for (int i = 0; i < M; i++)
      C[j][i] = D[j][i] = T(0);

  A.updateTo();
  B.updateTo();
  C.updateTo();

  T alpha(1);
  T beta(0);

  // U[X,Y] denotes a row-major matrix U with X rows and Y cols
  // element U(i,j) is located at: U.data() + sizeof(U_type) * (i*ldU + j)
  //
  // A,B,C,D are treated as row-major matrices, but the arguments to gemm are treated as col-major
  // so the call below to ompBLAS::gemm is equivalent to one of the following (with row-major matrices)
  // transa/transb == 'N'/'N':   C[N,M] = A[K,M] * B[N,K]; C = B * A
  // transa/transb == 'N'/'T':   C[N,M] = A[K,M] * B[K,N]; C = B^t * A
  // transa/transb == 'T'/'N':   C[N,M] = A[M,K] * B[N,K]; C = B * A^t
  // transa/transb == 'T'/'T':   C[N,M] = A[M,K] * B[K,N]; C = B^t * A^t

  compute::Queue<PL> queue;
  compute::BLASHandle<PlatformKind::CUDA> h_blas(queue);
  compute::BLAS::gemm(h_blas, transa, transb, M, N, K, alpha, A.device_data(), a1, B.device_data(), b1, beta, C.device_data(),
                M);
  queue.sync();

  C.updateFrom();

  BLAS::gemm(transa, transb, M, N, K, alpha, A.data(), a1, B.data(), b1, beta, D.data(), M);

  for (int j = 0; j < N; j++)
    for (int i = 0; i < M; i++)
    {
      CHECK(std::real(C[j][i]) == Approx(std::real(D[j][i])));
      CHECK(std::imag(C[j][i]) == Approx(std::imag(D[j][i])));
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

TEST_CASE("AccelBLAS", "[BLAS]")
{
  SECTION("gemm")
  {
#if defined(ENABLE_CUDA)
    test_gemm_cases<PlatformKind::CUDA>();
#endif
  }
}
} // namespace qmcplusplus
