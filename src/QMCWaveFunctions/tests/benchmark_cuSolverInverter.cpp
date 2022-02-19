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
#include "Utilities/for_testing/checkMatrix.hpp"
#include "Containers/tests/makeRngSpdMatrix.hpp"

namespace qmcplusplus
{

TEST_CASE("cuSolverInverter_bench", "[wavefunction][benchmark]")
{
#ifdef QMC_COMPLEX
  using FullPrecValueType = std::complex<double>;
#else
  using FullPrecValueType = double;
#endif

  cuSolverInverter<FullPrecValueType> solver;

  int N = 1024;

  Matrix<FullPrecValueType> m(N, N);
  Matrix<FullPrecValueType> m_invT(N, N);
  Matrix<FullPrecValueType> m_invT_CPU(N, N);
  Matrix<FullPrecValueType, CUDAAllocator<FullPrecValueType>> m_invGPU;
  std::complex<double> log_value;
  m.resize(N, N);
  m_invT.resize(N, N);
  m_invT_CPU.resize(N, N);
  m_invGPU.resize(N, N);

  testing::MakeRngSpdMatrix<FullPrecValueType> makeRngSpdMatrix{};
  makeRngSpdMatrix(m);

  BENCHMARK("cuSolverInverter") { solver.invert_transpose(m, m_invT, m_invGPU, log_value); };

  DiracMatrix<FullPrecValueType> dmat;
  BENCHMARK("CPU") { dmat.invert_transpose(m, m_invT_CPU, log_value); };

  auto check_matrix_result = checkMatrix(m_invT, m_invT_CPU);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
}

} // namespace qmcplusplus
