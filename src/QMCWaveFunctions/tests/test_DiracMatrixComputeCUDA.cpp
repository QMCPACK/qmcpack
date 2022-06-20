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

#include <catch.hpp>
#include <algorithm>
#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Fermion/DiracMatrixComputeCUDA.hpp"
#include "makeRngSpdMatrix.hpp"
#include "Utilities/for_testing/checkMatrix.hpp"
#include "Utilities/for_testing/RandomForTest.h"
#include "Platforms/DualAllocatorAliases.hpp"
#include "Platforms/CUDA/CUDALinearAlgebraHandles.h"

// Legacy CPU inversion for temporary testing
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"


namespace qmcplusplus
{
template<typename T>
using OffloadPinnedMatrix = Matrix<T, PinnedDualAllocator<T>>;
template<typename T>
using OffloadPinnedVector = Vector<T, PinnedDualAllocator<T>>;

TEST_CASE("DiracMatrixComputeCUDA_cuBLAS_geam_call", "[wavefunction][fermion]")
{
  OffloadPinnedMatrix<double> mat_a;
  int n = 4;
  mat_a.resize(n, n);
  OffloadPinnedMatrix<double> temp_mat;
  temp_mat.resize(n, n);
  OffloadPinnedMatrix<double> mat_c;
  mat_c.resize(n, n);

  double host_one(1.0);
  double host_zero(0.0);

  std::vector<double> A{2, 5, 8, 7, 5, 2, 2, 8, 7, 5, 6, 6, 5, 4, 4, 8};
  std::copy_n(A.begin(), 16, mat_a.data());
  CUDALinearAlgebraHandles cuda_handles;
  int lda = n;
  cudaCheck(cudaMemcpyAsync((void*)(temp_mat.device_data()), (void*)(mat_a.data()), mat_a.size() * sizeof(double),
                            cudaMemcpyHostToDevice, cuda_handles.hstream));
  cublasErrorCheck(cuBLAS::geam(cuda_handles.h_cublas, CUBLAS_OP_T, CUBLAS_OP_N, n, n, &host_one,
                                temp_mat.device_data(), lda, &host_zero, mat_c.device_data(), lda, mat_a.device_data(),
                                lda),
                   "cuBLAS::geam failed.");
}

TEST_CASE("DiracMatrixComputeCUDA_different_batch_sizes", "[wavefunction][fermion]")
{
  OffloadPinnedMatrix<double> mat_a;
  mat_a.resize(4, 4);
  std::vector<double> A{2, 5, 8, 7, 5, 2, 2, 8, 7, 5, 6, 6, 5, 4, 4, 8};
  std::copy_n(A.data(), 16, mat_a.data());
  OffloadPinnedVector<std::complex<double>> log_values;
  log_values.resize(1);
  OffloadPinnedMatrix<double> inv_mat_a;
  inv_mat_a.resize(4, 4);
  CUDALinearAlgebraHandles cuda_handles;
  DiracMatrixComputeCUDA<double> dmcc;

  dmcc.invert_transpose(cuda_handles, mat_a, inv_mat_a, log_values);


  OffloadPinnedMatrix<double> mat_b;
  mat_b.resize(4, 4);
  double invA[16]{-0.08247423, -0.26804124, 0.26804124, 0.05154639,  0.18556701,  -0.89690722, 0.39690722,  0.13402062,
                  0.24742268,  -0.19587629, 0.19587629, -0.15463918, -0.29896907, 1.27835052,  -0.77835052, 0.06185567};
  std::copy_n(invA, 16, mat_b.data());

  auto check_matrix_result = checkMatrix(inv_mat_a, mat_b);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }

  OffloadPinnedMatrix<double> mat_a2;
  mat_a2.resize(4, 4);
  std::copy_n(A.begin(), 16, mat_a2.data());
  OffloadPinnedMatrix<double> inv_mat_a2;
  inv_mat_a2.resize(4, 4);

  RefVector<const OffloadPinnedMatrix<double>> a_mats{mat_a, mat_a2};
  RefVector<OffloadPinnedMatrix<double>> inv_a_mats{inv_mat_a, inv_mat_a2};

  log_values.resize(2);
  dmcc.mw_invertTranspose(cuda_handles, a_mats, inv_a_mats, log_values);

  check_matrix_result = checkMatrix(inv_mat_a, mat_b);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
  check_matrix_result = checkMatrix(inv_mat_a2, mat_b);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }

  CHECK(log_values[0] == ComplexApprox(std::complex<double>{5.267858159063328, 6.283185307179586}));
  CHECK(log_values[1] == ComplexApprox(std::complex<double>{5.267858159063328, 6.283185307179586}));

  OffloadPinnedMatrix<double> mat_a3;
  mat_a3.resize(4, 4);
  std::copy_n(A.begin(), 16, mat_a3.data());
  OffloadPinnedMatrix<double> inv_mat_a3;
  inv_mat_a3.resize(4, 4);

  a_mats[1] = mat_a3;

  RefVector<const OffloadPinnedMatrix<double>> a_mats3{mat_a, mat_a2, mat_a3};
  RefVector<OffloadPinnedMatrix<double>> inv_a_mats3{inv_mat_a, inv_mat_a2, inv_mat_a3};

  log_values.resize(3);
  dmcc.mw_invertTranspose(cuda_handles, a_mats3, inv_a_mats3, log_values);

  check_matrix_result = checkMatrix(inv_mat_a, mat_b);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
  check_matrix_result = checkMatrix(inv_mat_a2, mat_b);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
  check_matrix_result = checkMatrix(inv_mat_a3, mat_b);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }

  CHECK(log_values[0] == ComplexApprox(std::complex<double>{5.267858159063328, 6.283185307179586}));
  CHECK(log_values[1] == ComplexApprox(std::complex<double>{5.267858159063328, 6.283185307179586}));
  CHECK(log_values[2] == ComplexApprox(std::complex<double>{5.267858159063328, 6.283185307179586}));
}

TEST_CASE("DiracMatrixComputeCUDA_complex_determinants_against_legacy", "[wavefunction][fermion]")
{
  int n = 64;
  CUDALinearAlgebraHandles cuda_handles;

  DiracMatrixComputeCUDA<std::complex<double>> dmcc;

  Matrix<std::complex<double>> mat_spd;
  mat_spd.resize(n, n);
  testing::MakeRngSpdMatrix<std::complex<double>> makeRngSpdMatrix;
  makeRngSpdMatrix(mat_spd);
  // You would hope you could do this
  // OffloadPinnedMatrix<double> mat_a(mat_spd);
  // But you can't
  OffloadPinnedMatrix<std::complex<double>> mat_a(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      mat_a(i, j) = mat_spd(i, j);

  Matrix<std::complex<double>> mat_spd2;
  mat_spd2.resize(n, n);
  makeRngSpdMatrix(mat_spd2);
  // You would hope you could do this
  // OffloadPinnedMatrix<double> mat_a(mat_spd);
  // But you can't
  OffloadPinnedMatrix<std::complex<double>> mat_a2(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      mat_a2(i, j) = mat_spd2(i, j);

  OffloadPinnedVector<std::complex<double>> log_values;
  log_values.resize(2);
  OffloadPinnedMatrix<std::complex<double>> inv_mat_a;
  inv_mat_a.resize(n, n);
  OffloadPinnedMatrix<std::complex<double>> inv_mat_a2;
  inv_mat_a2.resize(n, n);

  RefVector<const OffloadPinnedMatrix<std::complex<double>>> a_mats{mat_a, mat_a2};
  RefVector<OffloadPinnedMatrix<std::complex<double>>> inv_a_mats{inv_mat_a, inv_mat_a2};

  dmcc.mw_invertTranspose(cuda_handles, a_mats, inv_a_mats, log_values);

  DiracMatrix<std::complex<double>> dmat;
  Matrix<std::complex<double>> inv_mat_test(n, n);
  std::complex<double> det_log_value;
  dmat.invert_transpose(mat_spd, inv_mat_test, det_log_value);

  auto check_matrix_result = checkMatrix(inv_mat_a, inv_mat_test);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }

  dmat.invert_transpose(mat_spd2, inv_mat_test, det_log_value);
  check_matrix_result = checkMatrix(inv_mat_a2, inv_mat_test);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
}

TEST_CASE("DiracMatrixComputeCUDA_large_determinants_against_legacy", "[wavefunction][fermion]")
{
  int n = 64;
  CUDALinearAlgebraHandles cuda_handles;
  DiracMatrixComputeCUDA<double> dmcc;

  Matrix<double> mat_spd;
  mat_spd.resize(n, n);
  testing::MakeRngSpdMatrix<double> makeRngSpdMatrix;
  makeRngSpdMatrix(mat_spd);
  // You would hope you could do this
  // OffloadPinnedMatrix<double> mat_a(mat_spd);
  // But you can't
  OffloadPinnedMatrix<double> mat_a(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      mat_a(i, j) = mat_spd(i, j);

  Matrix<double> mat_spd2;
  mat_spd2.resize(n, n);
  makeRngSpdMatrix(mat_spd2);
  // You would hope you could do this
  // OffloadPinnedMatrix<double> mat_a(mat_spd);
  // But you can't
  OffloadPinnedMatrix<double> mat_a2(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      mat_a2(i, j) = mat_spd2(i, j);

  OffloadPinnedVector<std::complex<double>> log_values;
  log_values.resize(2);
  OffloadPinnedMatrix<double> inv_mat_a;
  inv_mat_a.resize(n, n);
  OffloadPinnedMatrix<double> inv_mat_a2;
  inv_mat_a2.resize(n, n);

  RefVector<const OffloadPinnedMatrix<double>> a_mats{mat_a, mat_a2};
  RefVector<OffloadPinnedMatrix<double>> inv_a_mats{inv_mat_a, inv_mat_a2};

  dmcc.mw_invertTranspose(cuda_handles, a_mats, inv_a_mats, log_values);

  DiracMatrix<double> dmat;
  Matrix<double> inv_mat_test(n, n);
  std::complex<double> det_log_value;
  dmat.invert_transpose(mat_spd, inv_mat_test, det_log_value);

  auto check_matrix_result = checkMatrix(inv_mat_a, inv_mat_test);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }

  dmat.invert_transpose(mat_spd2, inv_mat_test, det_log_value);
  check_matrix_result = checkMatrix(inv_mat_a2, inv_mat_test);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
}

} // namespace qmcplusplus
