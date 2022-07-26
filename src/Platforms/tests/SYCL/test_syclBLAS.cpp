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
#include "SYCL/SYCLruntime.hpp"
#include "SYCL/SYCLallocator.hpp"
#include "SYCL/syclBLAS.hpp"
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include "CPU/BLAS.hpp"

namespace qmcplusplus
{
template<typename T, typename Alloc>
void test_gemv(const int M_b, const int N_b, const char trans)
{
  const int M = trans == 'T' ? M_b : N_b;
  const int N = trans == 'T' ? N_b : M_b;

  using vec_t = Vector<T, Alloc>;
  using mat_t = Matrix<T, Alloc>;

  sycl::queue handle = getSYCLDefaultDeviceDefaultQueue();

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
    C[i] = D[i] = T(-0.1);

  A.updateTo();
  B.updateTo();

  T alpha(1);
  T beta(0);

  // in Fortran, B[M][N] is viewed as B^T
  // when trans == 'T', the actual calculation is B * A[N] = C[M]
  // when trans == 'N', the actual calculation is B^T * A[M] = C[N]
  //ompBLAS::gemv(handle, trans, N_b, M_b, alpha, B.device_data(), N_b, A.device_data(), 1, beta, C.device_data(), 1);

  syclBLAS::gemv(handle, trans, N_b, M_b, alpha, B.device_data(), N_b, A.device_data(), 1, beta, C.device_data(), 1)
      .wait();

  C.updateFrom();

  if (trans == 'T')
    BLAS::gemv_trans(M_b, N_b, B.data(), A.data(), D.data());
  else
    BLAS::gemv(M_b, N_b, B.data(), A.data(), D.data());

  for (int index = 0; index < M; index++)
    CHECK(C[index] == D[index]);
}

TEST_CASE("OmpSYCL gemv", "[SYCL]")
{
  const int M           = 137;
  const int N           = 79;
  const int batch_count = 23;

  // Non-batched test
  std::cout << "Testing TRANS gemv" << std::endl;
#if defined(ENABLE_OFFLOAD)
  test_gemv<float, OMPallocator<float>>(M, N, 'T');
  test_gemv<double, OMPallocator<double>>(M, N, 'T');
#endif
}

} // namespace qmcplusplus
