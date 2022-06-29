//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include <memory>
#include <vector>
#include <iostream>
#include "OMPTarget/OMPallocator.hpp"
#include "OMPTarget/ompBLAS.hpp"
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include <CPU/BLAS.hpp>

namespace qmcplusplus
{

template<typename T>
void test_gemv(const int M_b, const int N_b, const char trans)
{
  const int M = trans == 'T' ? M_b : N_b;
  const int N = trans == 'T' ? N_b : M_b;

  using vec_t = Vector<T, OMPallocator<T>>;
  using mat_t = Matrix<T, OMPallocator<T>>;

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

template<typename T>
void test_gemv_batched(const int M_b, const int N_b, const char trans, const int batch_count)
{
  const int M = trans == 'T' ? M_b : N_b;
  const int N = trans == 'T' ? N_b : M_b;

  using vec_t = Vector<T, OMPallocator<T>>;
  using mat_t = Matrix<T, OMPallocator<T>>;

  ompBLAS::ompBLAS_handle handle;

  // Create input vector
  std::vector<vec_t> As;
  Vector<const T*, OMPallocator<const T*>> Aptrs;

  // Create input matrix
  std::vector<mat_t> Bs;
  Vector<const T*, OMPallocator<const T*>> Bptrs;

  // Create output vector (ompBLAS)
  std::vector<vec_t> Cs;
  Vector<T*, OMPallocator<T*>> Cptrs;

  // Create output vector (BLAS)
  std::vector<vec_t> Ds;
  Vector<T*, OMPallocator<T*>> Dptrs;

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
    handle = batch;

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
  Vector<T, OMPallocator<T>> alpha(batch_count);
  Vector<T, OMPallocator<T>> beta(batch_count);

  for (int batch = 0; batch < batch_count; batch++)
  {
    alpha[batch] = T(1);
    beta[batch]  = T(0);
  }

  alpha.updateTo();
  beta.updateTo();

  ompBLAS::gemv_batched(handle, trans, N_b, M_b, alpha.device_data(), Bptrs.device_data(), N_b, Aptrs.device_data(), 1,
                        beta.device_data(), Cptrs.device_data(), 1, batch_count);

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

TEST_CASE("OmpBLAS gemv", "[OMP]")
{
  const int M           = 137;
  const int N           = 79;
  const int batch_count = 23;

  // Non-batched test
  std::cout << "Testing TRANS gemv" << std::endl;
  test_gemv<float>(M, N, 'T');
  test_gemv<double>(M, N, 'T');
#if defined(QMC_COMPLEX)
  test_gemv<std::complex<float>>(N, M, 'T');
  test_gemv<std::complex<double>>(N, M, 'T');
#endif
  // Batched Test
  std::cout << "Testing TRANS gemv_batched" << std::endl;
  test_gemv_batched<float>(M, N, 'T', batch_count);
  test_gemv_batched<double>(M, N, 'T', batch_count);
#if defined(QMC_COMPLEX)
  test_gemv_batched<std::complex<float>>(N, M, 'T', batch_count);
  test_gemv_batched<std::complex<double>>(N, M, 'T', batch_count);
#endif
}

TEST_CASE("OmpBLAS gemv notrans", "[OMP]")
{
  const int M           = 137;
  const int N           = 79;
  const int batch_count = 23;

  // Non-batched test
  std::cout << "Testing NOTRANS gemv" << std::endl;
  test_gemv<float>(M, N, 'N');
  test_gemv<double>(M, N, 'N');
#if defined(QMC_COMPLEX)
  test_gemv<std::complex<float>>(N, M, 'N');
  test_gemv<std::complex<double>>(N, M, 'N');
#endif
  // Batched Test
  std::cout << "Testing NOTRANS gemv_batched" << std::endl;
  test_gemv_batched<float>(M, N, 'N', batch_count);
  test_gemv_batched<double>(M, N, 'N', batch_count);
#if defined(QMC_COMPLEX)
  test_gemv_batched<std::complex<float>>(N, M, 'N', batch_count);
  test_gemv_batched<std::complex<double>>(N, M, 'N', batch_count);
#endif
}

template<typename T>
void test_ger(const int M, const int N)
{
  using vec_t = Vector<T, OMPallocator<T>>;
  using mat_t = Matrix<T, OMPallocator<T>>;

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

template<typename T>
void test_ger_batched(const int M, const int N, const int batch_count)
{
  using vec_t = Vector<T, OMPallocator<T>>;
  using mat_t = Matrix<T, OMPallocator<T>>;

  ompBLAS::ompBLAS_handle handle;

  // Create input vector
  std::vector<vec_t> Xs;
  Vector<const T*, OMPallocator<const T*>> Xptrs;
  std::vector<vec_t> Ys;
  Vector<const T*, OMPallocator<const T*>> Yptrs;

  // Create input matrix
  std::vector<mat_t> Ahs;
  Vector<T*, OMPallocator<T*>> Ahptrs;
  std::vector<mat_t> Ads;
  Vector<T*, OMPallocator<T*>> Adptrs;

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
    handle = batch;

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
  Vector<T, OMPallocator<T>> alpha(batch_count);

  for (int batch = 0; batch < batch_count; batch++)
  {
    alpha[batch] = T(1);
  }

  alpha.updateTo();

  ompBLAS::ger_batched(handle, M, N, alpha.device_data(), Xptrs.device_data(), 1, Yptrs.device_data(), 1,
                       Adptrs.device_data(), M, batch_count);

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

TEST_CASE("OmpBLAS ger", "[OMP]")
{
  const int M           = 137;
  const int N           = 79;
  const int batch_count = 23;

  // Non-batched test
  std::cout << "Testing ger" << std::endl;
  test_ger<float>(M, N);
  test_ger<double>(M, N);
#if defined(QMC_COMPLEX)
  test_ger<std::complex<float>>(N, M);
  test_ger<std::complex<double>>(N, M);
#endif
  // Batched Test
  std::cout << "Testing ger_batched" << std::endl;
  test_ger_batched<float>(M, N, batch_count);
  test_ger_batched<double>(M, N, batch_count);
#if defined(QMC_COMPLEX)
  test_ger_batched<std::complex<float>>(N, M, batch_count);
  test_ger_batched<std::complex<double>>(N, M, batch_count);
#endif
}
} // namespace qmcplusplus
