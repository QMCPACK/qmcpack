//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Configuration.h"
#include "Common/Queue.hpp"
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"
#include "QMCWaveFunctions/Fermion/AccelEngine.hpp"
#include "checkMatrix.hpp"
#include "createTestMatrix.h"

namespace qmcplusplus
{

template<PlatformKind PL, typename T>
using DeviceAllocator = typename compute::MemManage<PL>::template DeviceAllocator<T>;

template<PlatformKind PL, typename T, typename FP_T>
using SolverInverter = typename AccelEngine<PL, T, FP_T>::Inverter;

template<PlatformKind PL, typename T, typename FP_T>
void test_solver(const size_t N)
{
  using Value    = T;
  using LogValue = std::complex<double>;

  SolverInverter<PL, T, FP_T> solver;
  compute::Queue<PL> queue;

  Matrix<Value> m(N, N);
  Matrix<Value> m_invT(N, N);
  Matrix<Value, DeviceAllocator<PL, Value>> m_invGPU(N, N);
  LogValue log_value;

  SECTION("identity")
  {
    fillIdentityMatrix(m);

    solver.invert_transpose(m, m_invT, m_invGPU, log_value, queue.getNative());
    CHECK(log_value == LogComplexApprox(0.0));

    Matrix<Value> eye;
    eye.resize(3, 3);
    fillIdentityMatrix(eye);

    CheckMatrixResult check_matrix_result = checkMatrix(m_invT, eye);
    CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
  }

  SECTION("3x3 matrix")
  {
    Matrix<Value> a(N, N);
    Matrix<Value> a_T(N, N);
    TestMatrix1::fillInput(a);

    simd::transpose(a.data(), a.rows(), a.cols(), a_T.data(), a_T.rows(), a_T.cols());
    solver.invert_transpose(a_T, m_invT, m_invGPU, log_value, queue.getNative());
    CHECK(log_value == LogComplexApprox(TestMatrix1::logDet()));

    Matrix<Value> b(3, 3);

    TestMatrix1::fillInverse(b);

    auto check_matrix_result = checkMatrix(m_invT, b);
    CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
  }
}

TEMPLATE_TEST_CASE("cuSolverInverter", "[wavefunction][fermion]", double, float)
{
  // TestType is defined by Catch. It is the type in each instantiation of the templated test case.
#ifdef QMC_COMPLEX
  using FullPrecValue = std::complex<double>;
  using Value         = std::complex<TestType>;
#else
  using FullPrecValue = double;
  using Value         = TestType;
#endif
#if defined(ENABLE_CUDA)
  test_solver<PlatformKind::CUDA, Value, FullPrecValue>(3);
#endif
#if defined(ENABLE_SYCL)
  test_solver<PlatformKind::SYCL, Value, FullPrecValue>(3);
#endif
}

} // namespace qmcplusplus
