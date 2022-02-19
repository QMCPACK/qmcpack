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
#include "Platforms/CPU/SIMD/aligned_allocator.hpp"
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"
#include "QMCWaveFunctions/Fermion/cuSolverInverter.hpp"
#include "checkMatrix.hpp"
#include "createTestMatrix.h"

namespace qmcplusplus
{

using LogValueType = std::complex<double>;

TEMPLATE_TEST_CASE("cuSolverInverter", "[wavefunction][fermion]", double, float)
{
#ifdef QMC_COMPLEX
  using FullPrecValueType = std::complex<double>;
  using ValueType         = std::complex<TestType>;
#else
  using FullPrecValueType = double;
  using ValueType         = TestType;
#endif
  cuSolverInverter<FullPrecValueType> solver;
  const int N = 3;

  Matrix<ValueType> m(N, N);
  Matrix<ValueType> m_invT(N, N);
  Matrix<ValueType, CUDAAllocator<ValueType>> m_invGPU(N, N);
  LogValueType log_value;

  SECTION("identity")
  {
    fillIdentityMatrix(m);

    solver.invert_transpose(m, m_invT, m_invGPU, log_value);
    REQUIRE(log_value == LogComplexApprox(0.0));

    Matrix<ValueType> eye;
    eye.resize(3, 3);
    fillIdentityMatrix(eye);

    CheckMatrixResult check_matrix_result = checkMatrix(m_invT, eye);
    CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
  }

  SECTION("3x3 matrix")
  {
    Matrix<ValueType> a(N, N);
    Matrix<ValueType> a_T(N, N);
    fillTestMatrix1_Input(a);

    simd::transpose(a.data(), a.rows(), a.cols(), a_T.data(), a_T.rows(), a_T.cols());
    solver.invert_transpose(a_T, m_invT, m_invGPU, log_value);
    REQUIRE(log_value == LogComplexApprox(testMatrix1_logDet()));

    Matrix<ValueType> b(3, 3);

    fillTestMatrix1_Inverse(b);

    auto check_matrix_result = checkMatrix(m_invT, b);
    CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
  }
}

} // namespace qmcplusplus
