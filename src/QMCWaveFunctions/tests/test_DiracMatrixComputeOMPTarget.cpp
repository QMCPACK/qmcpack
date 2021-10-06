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
#include "Configuration.h"
#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Fermion/DiracMatrixComputeOMPTarget.hpp"
#include "makeRngSpdMatrix.hpp"
#include "Utilities/Resource.h"
#include "Utilities/for_testing/checkMatrix.hpp"
#include "Utilities/for_testing/RandomForTest.h"
#include "Platforms/PinnedAllocator.h"

// Legacy CPU inversion for temporary testing
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"


namespace qmcplusplus
{
template<typename T>
using OffloadPinnedAllocator = OMPallocator<T, PinnedAlignedAllocator<T>>;

template<typename T>
using OffloadPinnedMatrix = Matrix<T, OffloadPinnedAllocator<T>>;
template<typename T>
using OffloadPinnedVector = Vector<T, OffloadPinnedAllocator<T>>;

TEST_CASE("DiracMatrixComputeOMPTarget_different_batch_sizes", "[wavefunction][fermion]")
{
  OffloadPinnedMatrix<double> mat_a;
  mat_a.resize(4, 4);
  std::vector<double> A{2, 5, 8, 7, 5, 2, 2, 8, 7, 5, 6, 6, 5, 4, 4, 8};
  std::copy_n(A.data(), 16, mat_a.data());
  OffloadPinnedVector<std::complex<double>> log_values;
  log_values.resize(1);
  OffloadPinnedMatrix<double> inv_mat_a;
  inv_mat_a.resize(4, 4);
  DiracMatrixComputeOMPTarget<double> dmc_omp;

  DummyResource dummy;
  std::complex<double> log_value;
  dmc_omp.invert_transpose(dummy, mat_a, inv_mat_a, log_value);
  CHECK(log_value == LogComplexApprox(std::complex<double>{5.267858159063328, 6.283185307179586}));


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
  DummyResource dummy_res;
  dmc_omp.mw_invertTranspose(dummy_res, a_mats, inv_a_mats, log_values);

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
  dmc_omp.mw_invertTranspose(dummy_res, a_mats3, inv_a_mats3, log_values);

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

TEST_CASE("DiracMatrixComputeOMPTarget_large_determinants_against_legacy", "[wavefunction][fermion]")
{
  int n = 64;

  DiracMatrixComputeOMPTarget<double> dmc_omp;

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


  DummyResource dummy_res;
  dmc_omp.mw_invertTranspose(dummy_res, a_mats, inv_a_mats, log_values);

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
