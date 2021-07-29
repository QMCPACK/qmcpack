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
#include <cmath>
#include "Configuration.h"
#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Fermion/DiracMatrixComputeCUDA.hpp"
#include "makeRngSpdMatrix.hpp"
#include "Utilities/for_testing/checkMatrix.hpp"
#include "Utilities/for_testing/RandomForTest.h"
#include "Platforms/PinnedAllocator.h"
#include "Platforms/CUDA/CUDALinearAlgebraHandles.h"
#include "Platforms/tests/CUDA/test_device_value_kernels.hpp"
#include "type_traits/type_tests.hpp"
// Legacy CPU inversion for temporary testing
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"

#ifdef ENABLE_CUDA
#include "DualAllocator.hpp"
#endif

namespace qmcplusplus
{

// The idea is to test both of these when possible without rebuilding. So unlike other parts of the code
// we do not use an ifdef'd type alias to swap the allocator.
#ifdef ENABLE_OFFLOAD
  template<typename T>
  using DualSpacePinnedAllocator = OMPallocator<T, PinnedAlignedAllocator<T>>;
#elif defined(ENABLE_CUDA)
  template<typename T>
  using DualSpacePinnedAllocator = DualAllocator<T, CUDAAllocator<T>, PinnedAlignedAllocator<T>>;
#endif

template<typename T>
using OffloadPinnedMatrix = Matrix<T, DualSpacePinnedAllocator<T>>;
template<typename T>
using OffloadPinnedVector = Vector<T, DualSpacePinnedAllocator<T>>;

template<typename T, IsComplex<T> = true>
void checkLogDetValue(T val_a, T val_b) {
  CHECK(val_a.real() == Approx(val_b.real()));
  using Real = RealAlias<T>;
  constexpr Real two_pi = 2 * M_PI;
  CHECK( std::fmod(val_a.imag(), two_pi) == Approx(std::fmod(val_b.imag(), two_pi)));
}

template<typename T, IsReal<T> = true>
void checkLogDetValue(T val_a, T val_b) {
  CHECK(val_a == ComplexApprox(val_b));
}

template<typename VALUE, typename VALUE_FP>
void test_DiracMatrixComputeCUDA_different_batch_sizes() {
  using FullPrecReal = RealAlias<VALUE_FP>;

  OffloadPinnedMatrix<VALUE> mat_a;
  mat_a.resize(4, 4);
  std::vector<VALUE> A{2, 5, 8, 7, 5, 2, 2, 8, 7, 5, 6, 6, 5, 4, 4, 8};
  std::copy_n(A.data(), 16, mat_a.data());
  OffloadPinnedVector<std::complex<FullPrecReal>> log_values;
  log_values.resize(1);
  OffloadPinnedMatrix<VALUE> inv_mat_a;
  inv_mat_a.resize(4, 4);
  auto cuda_handles = std::make_unique<CUDALinearAlgebraHandles>();
  DiracMatrixComputeCUDA<VALUE_FP> dmcc(cuda_handles->hstream);

  dmcc.invert_transpose(*cuda_handles, mat_a, inv_mat_a, log_values);
  checkLogDetValue(log_values[0], std::complex<FullPrecReal>{5.267858159063328, 6.283185307179586} );

  OffloadPinnedMatrix<VALUE> mat_b;
  mat_b.resize(4, 4);
  double invA[16]{-0.08247423, -0.26804124, 0.26804124, 0.05154639,  0.18556701,  -0.89690722, 0.39690722,  0.13402062,
                  0.24742268,  -0.19587629, 0.19587629, -0.15463918, -0.29896907, 1.27835052,  -0.77835052, 0.06185567};
  std::copy_n(invA, 16, mat_b.data());

  auto check_matrix_result = checkMatrix(inv_mat_a, mat_b);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }

  OffloadPinnedMatrix<VALUE> mat_a2;
  mat_a2.resize(4, 4);
  std::copy_n(A.begin(), 16, mat_a2.data());
  OffloadPinnedMatrix<VALUE> inv_mat_a2;
  inv_mat_a2.resize(4, 4);

  RefVector<OffloadPinnedMatrix<VALUE>> a_mats{mat_a, mat_a2};
  RefVector<OffloadPinnedMatrix<VALUE>> inv_a_mats{inv_mat_a, inv_mat_a2};

  log_values.resize(2);
  dmcc.mw_invertTranspose(*cuda_handles, a_mats, inv_a_mats, log_values, {true, true});

  check_matrix_result = checkMatrix(inv_mat_a, mat_b);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
  check_matrix_result = checkMatrix(inv_mat_a2, mat_b);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }

  checkLogDetValue(log_values[0], std::complex<FullPrecReal>{5.267858159063328, 6.283185307179586});
  checkLogDetValue(log_values[1], std::complex<FullPrecReal>{5.267858159063328, 6.283185307179586});

  OffloadPinnedMatrix<VALUE> mat_a3;
  mat_a3.resize(4, 4);
  // mw_invertTranspose when VALUE == VALUE_FP has the side effect of updating mat_a with LU matrix.
  std::copy_n(A.begin(), 16, mat_a.data());
  std::copy_n(A.begin(), 16, mat_a2.data());
  std::copy_n(A.begin(), 16, mat_a3.data());
  OffloadPinnedMatrix<VALUE> inv_mat_a3;
  inv_mat_a3.resize(4, 4);

  a_mats[1] = mat_a3;

  RefVector<OffloadPinnedMatrix<VALUE>> a_mats3{mat_a, mat_a2, mat_a3};
  RefVector<OffloadPinnedMatrix<VALUE>> inv_a_mats3{inv_mat_a, inv_mat_a2, inv_mat_a3};

  log_values.resize(3);
  dmcc.mw_invertTranspose(*cuda_handles, a_mats3, inv_a_mats3, log_values, {true, true});

  check_matrix_result = checkMatrix(inv_mat_a, mat_b);
  CHECKED_ELSE(check_matrix_result.result) { std::cout << check_matrix_result.result_message << std::endl; }
  check_matrix_result = checkMatrix(inv_mat_a2, mat_b);
  CHECKED_ELSE(check_matrix_result.result) { std::cout << check_matrix_result.result_message << std::endl; }
  check_matrix_result = checkMatrix(inv_mat_a3, mat_b);
  CHECKED_ELSE(check_matrix_result.result) { std::cout << check_matrix_result.result_message << std::endl; }

  checkLogDetValue(log_values[0],std::complex<FullPrecReal>{5.267858159063328, 6.283185307179586});
  checkLogDetValue(log_values[1],std::complex<FullPrecReal>{5.267858159063328, 6.283185307179586});
  checkLogDetValue(log_values[2],std::complex<FullPrecReal>{5.267858159063328, 6.283185307179586});
}

template<typename VALUE, typename VALUE_FP, IsComplex<VALUE> = true>
void test_DiracMatrixComputeCUDA_complex() {
  using FullPrecReal = RealAlias<VALUE_FP>;

  OffloadPinnedMatrix<VALUE> mat_a;
  mat_a.resize(4, 4);
  std::vector<VALUE> A{ {2.0, 0.1}, {5.0, 0.1}, {8.0, 0.5}, {7.0, 1.0}, {5.0, 0.1}, {2.0, 0.2}, {2.0, 0.1}, {8.0, 0.5}, {7.0, 0.2}, {5.0, 1.0}, {6.0, -0.2}, {6.0, -0.2}, {5.0, 0.0}, {4.0, -0.1}, {4.0, -0.6}, {8.0, -2.0} };
  std::copy_n(A.data(), 16, mat_a.data());
  OffloadPinnedVector<std::complex<FullPrecReal>> log_values;
  log_values.resize(1);
  OffloadPinnedMatrix<VALUE> inv_mat_a;
  inv_mat_a.resize(4, 4);
  auto cuda_handles = std::make_unique<CUDALinearAlgebraHandles>();
  DiracMatrixComputeCUDA<VALUE_FP> dmcc(cuda_handles->hstream);

  dmcc.invert_transpose(*cuda_handles, mat_a, inv_mat_a, log_values);

  checkLogDetValue(log_values[0], std::complex<FullPrecReal>{5.6037775791955715, -6.158660333118822});

  OffloadPinnedMatrix<VALUE> mat_b;
  mat_b.resize(4, 4);
  double invA[16]{-0.08247423, -0.26804124, 0.26804124, 0.05154639,  0.18556701,  -0.89690722, 0.39690722,  0.13402062,
                  0.24742268,  -0.19587629, 0.19587629, -0.15463918, -0.29896907, 1.27835052,  -0.77835052, 0.06185567};
  std::copy_n(invA, 16, mat_b.data());

  auto check_matrix_result = checkMatrix(inv_mat_a, mat_b);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }

}


TEST_CASE("DiracMatrixComputeCUDA_cuBLAS_geam_call", "[wavefunction][fermion]")
{
  using Value = QMCTraits::ValueType;

  OffloadPinnedMatrix<Value> mat_a;
  int n = 4;
  mat_a.resize(n, n);
  OffloadPinnedMatrix<Value> temp_mat;
  temp_mat.resize(n, n);
  OffloadPinnedMatrix<Value> mat_c;
  mat_c.resize(n, n);

  Value host_one(1.0);
  Value host_zero(0.0);

  std::vector<Value> A{2, 5, 8, 7, 5, 2, 2, 8, 7, 5, 6, 6, 5, 4, 4, 8};
  std::copy_n(A.begin(), 16, mat_a.data());
  auto cuda_handles = std::make_unique<CUDALinearAlgebraHandles>();
  int lda= n;
  cudaCheck(cudaMemcpyAsync((void*)(temp_mat.device_data()), (void*)(mat_a.data()),
                            mat_a.size() * sizeof(Value), cudaMemcpyHostToDevice, cuda_handles->hstream));
  cublasErrorCheck(cuBLAS::geam(cuda_handles->h_cublas, CUBLAS_OP_T, CUBLAS_OP_N, n, n, &host_one,
                                    temp_mat.device_data(), lda, &host_zero,
                                mat_c.device_data(), lda, mat_a.device_data(), lda),
                   "cuBLAS::geam failed.");
}

TEST_CASE("DiracMatrixComputeCUDA_different_batch_sizes", "[wavefunction][fermion]")
{
  std::cout << "Testing DiracMatrixComputeCUDA<double, double>" << std::endl;
  test_DiracMatrixComputeCUDA_different_batch_sizes<double, double>();
  std::cout << "Testing DiracMatrixComputeCUDA<float, double>" << std::endl;
  test_DiracMatrixComputeCUDA_different_batch_sizes<float, double>();
  std::cout << "Testing DiracMatrixComputeCUDA<std::complex<double>, std::complex<double>>" << std::endl;
  test_DiracMatrixComputeCUDA_different_batch_sizes<std::complex<double>, std::complex<double>>();
  std::cout << "Testing DiracMatrixComputeCUDA<std::complex<float>, std::complex<double>>" << std::endl;
  test_DiracMatrixComputeCUDA_different_batch_sizes<std::complex<float>, std::complex<double>>();
}

TEST_CASE("DiracMatrixComputeCUDA_complex_determinants_against_legacy", "[wavefunction][fermion]")
{
  using Value = std::complex<double>;
  using RealValue = RealAlias<Value>;

  int n = 64;
  auto cuda_handles = std::make_unique<CUDALinearAlgebraHandles>();

  DiracMatrixComputeCUDA<Value> dmcc(cuda_handles->hstream);;

  Matrix<Value> mat_spd;
  mat_spd.resize(n, n);
  testing::MakeRngSpdMatrix<Value> makeRngSpdMatrix;
  makeRngSpdMatrix(mat_spd);
  // You would hope you could do this
  // OffloadPinnedMatrix<double> mat_a(mat_spd);
  // But you can't
  OffloadPinnedMatrix<Value> mat_a(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      mat_a(i, j) = mat_spd(i, j);

  Matrix<Value> mat_spd2;
  mat_spd2.resize(n, n);
  makeRngSpdMatrix(mat_spd2);
  // You would hope you could do this
  // OffloadPinnedMatrix<double> mat_a(mat_spd);
  // But you can't
  OffloadPinnedMatrix<Value> mat_a2(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      mat_a2(i, j) = mat_spd2(i, j);
  
  OffloadPinnedVector<std::complex<RealValue>> log_values;
  log_values.resize(2);
  OffloadPinnedMatrix<Value> inv_mat_a;
  inv_mat_a.resize(n, n);
  OffloadPinnedMatrix<Value> inv_mat_a2;
  inv_mat_a2.resize(n, n);

  RefVector<OffloadPinnedMatrix<Value>> a_mats{mat_a, mat_a2};
  RefVector<OffloadPinnedMatrix<Value>> inv_a_mats{inv_mat_a, inv_mat_a2};

  dmcc.mw_invertTranspose(*cuda_handles, a_mats, inv_a_mats, log_values, {true, true});

  DiracMatrix<Value> dmat;
  Matrix<Value> inv_mat_test(n, n);
  std::complex<RealValue> det_log_value;
  dmat.invert_transpose(mat_spd, inv_mat_test, det_log_value);
  auto check_matrix_result = checkMatrix(inv_mat_a, inv_mat_test);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
  checkLogDetValue(det_log_value, log_values[0]);
  
  dmat.invert_transpose(mat_spd2, inv_mat_test, det_log_value);
  check_matrix_result = checkMatrix(inv_mat_a2, inv_mat_test);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
  checkLogDetValue(det_log_value, log_values[1]);
}

TEST_CASE("DiracMatrixComputeCUDA_large_determinants_against_legacy", "[wavefunction][fermion]")
{
  int n = 64;
  auto cuda_handles = std::make_unique<CUDALinearAlgebraHandles>();

  using Value = double;

  DiracMatrixComputeCUDA<Value> dmcc(cuda_handles->hstream);;

  Matrix<Value> mat_spd;
  mat_spd.resize(n, n);
  testing::MakeRngSpdMatrix<Value> makeRngSpdMatrix;
  makeRngSpdMatrix(mat_spd);
  // You would hope you could do this
  // OffloadPinnedMatrix<Value> mat_a(mat_spd);
  // But you can't
  OffloadPinnedMatrix<Value> mat_a(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      mat_a(i, j) = mat_spd(i, j);

  Matrix<Value> mat_spd2;
  mat_spd2.resize(n, n);
  makeRngSpdMatrix(mat_spd2);
  // You would hope you could do this
  // OffloadPinnedMatrix<Value> mat_a(mat_spd);
  // But you can't
  OffloadPinnedMatrix<Value> mat_a2(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      mat_a2(i, j) = mat_spd2(i, j);
  
  OffloadPinnedVector<std::complex<Value>> log_values;
  log_values.resize(2);
  OffloadPinnedMatrix<Value> inv_mat_a;
  inv_mat_a.resize(n, n);
  OffloadPinnedMatrix<Value> inv_mat_a2;
  inv_mat_a2.resize(n, n);

  RefVector<OffloadPinnedMatrix<Value>> a_mats{mat_a, mat_a2};
  RefVector<OffloadPinnedMatrix<Value>> inv_a_mats{inv_mat_a, inv_mat_a2};

  dmcc.mw_invertTranspose(*cuda_handles, a_mats, inv_a_mats, log_values, {true, true});

  DiracMatrix<Value> dmat;
  Matrix<Value> inv_mat_test(n, n);
  std::complex<double> det_log_value;
  dmat.invert_transpose(mat_spd, inv_mat_test, det_log_value);
  
  auto check_matrix_result = checkMatrix(inv_mat_a, inv_mat_test);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }

  dmat.invert_transpose(mat_spd2, inv_mat_test, det_log_value);
  check_matrix_result = checkMatrix(inv_mat_a2, inv_mat_test);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
}

} // namespace qmcplusplus
