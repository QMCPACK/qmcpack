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
  compute::BLASHandle<PL> h_blas(queue);
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

  mat_t A2(a0, a1); // Input matrix
  mat_t B2(b0, b1); // Input matrix
  mat_t C2(N, M);   // Result matrix ompBLAS
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

  compute::BLAS::gemm_batched(h_blas, transa, transb, M, N, K, alpha, Aarr.device_data(), a1, Barr.device_data(), b1, beta,
                        Carr.device_data(), M, 2);
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

/*
template<typename T>
void test_gemv(const int M_b, const int N_b, const char trans)
{
  const int M = trans == 'T' ? M_b : N_b;
  const int N = trans == 'T' ? N_b : M_b;

  using vec_t = Vector<T, PinnedDualAllocator<T>>;
  using mat_t = Matrix<T, PinnedDualAllocator<T>>;

  ompBLAS::ompBLAS_handle handle;

  vec_t A(N);        // Input vector
  mat_t B(M_b, N_b); // Input matrix
  vec_t C(M);        // Result vector ompBLAS
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

  T alpha(1);
  T beta(0);

  // in Fortran, B[M][N] is viewed as B^T
  // when trans == 'T', the actual calculation is B * A[N] = C[M]
  // when trans == 'N', the actual calculation is B^T * A[M] = C[N]
  ompBLAS::gemv(handle, trans, N_b, M_b, alpha, B.device_data(), N_b, A.device_data(), 1, beta, C.device_data(), 1);
  C.updateFrom();

  if (trans == 'T')
    BLAS::gemv_trans(M_b, N_b, B.data(), A.data(), D.data());
  else
    BLAS::gemv(M_b, N_b, B.data(), A.data(), D.data());

  for (int index = 0; index < M; index++)
    CHECK(C[index] == D[index]);
}
*/

template<PlatformKind PL, typename T>
void test_one_gemv(const int M_b, const int N_b, const char trans, const int batch_count)
{
  const int M = trans == 'T' ? M_b : N_b;
  const int N = trans == 'T' ? N_b : M_b;

  using vec_t = Vector<T, PinnedDualAllocator<T>>;
  using mat_t = Matrix<T, PinnedDualAllocator<T>>;

  // Create input vector
  std::vector<vec_t> As;
  Vector<const T*, PinnedDualAllocator<const T*>> Aptrs;

  // Create input matrix
  std::vector<mat_t> Bs;
  Vector<const T*, PinnedDualAllocator<const T*>> Bptrs;

  // Create output vector (ompBLAS)
  std::vector<vec_t> Cs;
  Vector<T*, PinnedDualAllocator<T*>> Cptrs;

  // Create output vector (BLAS)
  std::vector<vec_t> Ds;
  Vector<T*, PinnedDualAllocator<T*>> Dptrs;

  // Resize pointer vectors
  Aptrs.resize(batch_count);
  Bptrs.resize(batch_count);
  Cptrs.resize(batch_count);
  Dptrs.resize(batch_count);

  // Resize data vectors
  As.resize(batch_count);
  Bs.resize(batch_count);
  Cs.resize(batch_count);
  Ds.resize(batch_count);

  // Fill data
  for (int batch = 0; batch < batch_count; batch++)
  {
    As[batch].resize(N);
    Aptrs[batch] = As[batch].device_data();

    Bs[batch].resize(M_b, N_b);
    Bptrs[batch] = Bs[batch].device_data();

    Cs[batch].resize(M);
    Cptrs[batch] = Cs[batch].device_data();

    Ds[batch].resize(M);
    Dptrs[batch] = Ds[batch].data();

    for (int i = 0; i < N; i++)
      As[batch][i] = i;

    for (int j = 0; j < M_b; j++)
      for (int i = 0; i < N_b; i++)
        Bs[batch][j][i] = i + j * 2;

    for (int i = 0; i < M; i++)
      Cs[batch][i] = Ds[batch][i] = T(0);

    As[batch].updateTo();
    Bs[batch].updateTo();
  }

  Aptrs.updateTo();
  Bptrs.updateTo();
  Cptrs.updateTo();

  // Run tests
  Vector<T, PinnedDualAllocator<T>> alpha(batch_count);
  Vector<T, PinnedDualAllocator<T>> beta(batch_count);
  Vector<T, PinnedDualAllocator<T>> beta1(batch_count);

  for (int batch = 0; batch < batch_count; batch++)
  {
    alpha[batch] = T(0.5);
    beta[batch]  = T(0);
    beta1[batch] = T(1);
  }

  alpha.updateTo();
  beta.updateTo();
  beta1.updateTo();

  compute::Queue<PL> queue;
  compute::BLASHandle<PL> h_blas(queue);
  // alpha 0.5, beta 0
  compute::BLAS::gemv_batched(h_blas, trans, N_b, M_b, alpha.device_data(), Bptrs.device_data(), N_b, Aptrs.device_data(), 1,
                        beta.device_data(), Cptrs.device_data(), 1, batch_count);
  // alpha 0.5, beta 1
  compute::BLAS::gemv_batched(h_blas, trans, N_b, M_b, alpha.device_data(), Bptrs.device_data(), N_b, Aptrs.device_data(), 1,
                        beta1.device_data(), Cptrs.device_data(), 1, batch_count);

  queue.sync();

  for (int batch = 0; batch < batch_count; batch++)
  {
    Cs[batch].updateFrom();
    if (trans == 'T')
      BLAS::gemv_trans(M_b, N_b, Bs[batch].data(), As[batch].data(), Ds[batch].data());
    else
      BLAS::gemv(M_b, N_b, Bs[batch].data(), As[batch].data(), Ds[batch].data());

    // Check results
    for (int index = 0; index < M; index++)
      CHECK(Cs[batch][index] == Ds[batch][index]);
  }
}

template<PlatformKind PL>
void test_gemv_cases()
{
  const int M           = 137;
  const int N           = 79;
  const int batch_count = 23;

  std::cout << "Testing NOTRANS gemv_batched" << std::endl;
  test_one_gemv<PL, float>(M, N, 'N', batch_count);
  test_one_gemv<PL, double>(M, N, 'N', batch_count);
#if defined(QMC_COMPLEX)
  test_one_gemv<PL, std::complex<float>>(N, M, 'N', batch_count);
  test_one_gemv<PL, std::complex<double>>(N, M, 'N', batch_count);
#endif
  std::cout << "Testing TRANS gemv_batched" << std::endl;
  test_one_gemv<PL, float>(M, N, 'T', batch_count);
  test_one_gemv<PL, double>(M, N, 'T', batch_count);
#if defined(QMC_COMPLEX)
  test_one_gemv<PL, std::complex<float>>(N, M, 'T', batch_count);
  test_one_gemv<PL, std::complex<double>>(N, M, 'T', batch_count);
#endif
}

/*
template<typename T>
void test_ger(const int M, const int N)
{
  using vec_t = Vector<T, PinnedDualAllocator<T>>;
  using mat_t = Matrix<T, PinnedDualAllocator<T>>;

  ompBLAS::ompBLAS_handle handle;

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

  T alpha(1);

  // in Fortran, B[M][N] is viewed as B^T
  ompBLAS::ger(handle, M, N, alpha, x.device_data(), 1, y.device_data(), 1, Ad.device_data(), M);
  Ad.updateFrom();

  BLAS::ger(M, N, alpha, x.data(), 1, y.data(), 1, Ah.data(), M);

  for (int j = 0; j < M; j++)
    for (int i = 0; i < N; i++)
      CHECK(Ah[j][i] == Ad[j][i]);
}
*/

template<PlatformKind PL, typename T>
void test_one_ger(const int M, const int N, const int batch_count)
{
  using vec_t = Vector<T, PinnedDualAllocator<T>>;
  using mat_t = Matrix<T, PinnedDualAllocator<T>>;

  // Create input vector
  std::vector<vec_t> Xs;
  Vector<const T*, PinnedDualAllocator<const T*>> Xptrs;
  std::vector<vec_t> Ys;
  Vector<const T*, PinnedDualAllocator<const T*>> Yptrs;

  // Create input matrix
  std::vector<mat_t> Ahs;
  Vector<T*, PinnedDualAllocator<T*>> Ahptrs;
  std::vector<mat_t> Ads;
  Vector<T*, PinnedDualAllocator<T*>> Adptrs;

  // Resize pointer vectors
  Xptrs.resize(batch_count);
  Yptrs.resize(batch_count);
  Ahptrs.resize(batch_count);
  Adptrs.resize(batch_count);

  // Resize data vectors
  Xs.resize(batch_count);
  Ys.resize(batch_count);
  Ahs.resize(batch_count);
  Ads.resize(batch_count);

  // Fill data
  for (int batch = 0; batch < batch_count; batch++)
  {
    Xs[batch].resize(M);
    Xptrs[batch] = Xs[batch].device_data();

    Ys[batch].resize(N);
    Yptrs[batch] = Ys[batch].device_data();

    Ads[batch].resize(M, N);
    Adptrs[batch] = Ads[batch].device_data();

    Ahs[batch].resize(M, N);
    Ahptrs[batch] = Ahs[batch].data();

    // Fill data
    for (int i = 0; i < M; i++)
      Xs[batch][i] = i;
    for (int i = 0; i < N; i++)
      Ys[batch][i] = N - i;

    for (int j = 0; j < M; j++)
      for (int i = 0; i < N; i++)
      {
        Ads[batch][j][i] = i + j * 2;
        Ahs[batch][j][i] = i + j * 2;
      }

    Xs[batch].updateTo();
    Ys[batch].updateTo();
    Ads[batch].updateTo();
  }

  Adptrs.updateTo();
  Xptrs.updateTo();
  Yptrs.updateTo();

  // Run tests
  Vector<T, PinnedDualAllocator<T>> alpha(batch_count);

  for (int batch = 0; batch < batch_count; batch++)
    alpha[batch] = T(0.5);

  alpha.updateTo();

  compute::Queue<PL> queue;
  compute::BLASHandle<PL> h_blas(queue);
  compute::BLAS::ger_batched(h_blas, M, N, alpha.device_data(), Xptrs.device_data(), 1, Yptrs.device_data(), 1,
                       Adptrs.device_data(), M, batch_count);
  queue.sync();

  for (int batch = 0; batch < batch_count; batch++)
  {
    Ads[batch].updateFrom();
    BLAS::ger(M, N, alpha[batch], Xs[batch].data(), 1, Ys[batch].data(), 1, Ahs[batch].data(), M);

    // Check results
    for (int j = 0; j < M; j++)
      for (int i = 0; i < N; i++)
        CHECK(Ads[batch][j][i] == Ahs[batch][j][i]);
  }
}

template<PlatformKind PL>
void test_ger_cases()
{
  const int M           = 137;
  const int N           = 79;
  const int batch_count = 23;

  // Batched Test
  std::cout << "Testing ger_batched" << std::endl;
  test_one_ger<PL, float>(M, N, batch_count);
  test_one_ger<PL, double>(M, N, batch_count);
#if defined(QMC_COMPLEX)
  test_one_ger<PL, std::complex<float>>(N, M, batch_count);
  test_one_ger<PL, std::complex<double>>(N, M, batch_count);
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
} // namespace qmcplusplus
