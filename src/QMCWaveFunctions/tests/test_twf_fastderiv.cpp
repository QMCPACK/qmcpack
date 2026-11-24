//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//////////////////////////////////////////////////////////////////////////////////////

/** @file test_twf_fastderiv.cpp
 *  @brief Main unit tests of the TWFFastDerivWrapper API.
 *
 *  Tests here call TWFFastDerivWrapper member functions directly on small hand-built
 *  matrices and check their behavior.
 *
 *  This is in contrast to test_twf_fastderiv_raii.cpp, which only tests how the wrapper's
 *  resources are created, acquired and released when it is managed by TrialWaveFunction.
 */

#include <catch2/catch_test_macros.hpp>

#include "QMCWaveFunctions/TWFFastDerivWrapper.h"

namespace qmcplusplus
{
// With no virtual orbitals (norb == nelec) buildIntermediates must not issue a GEMM with
// M=0/LDC=0, which violates LDC >= max(1,M) ("Parameter 13 was incorrect on entry to ZGEMM").
TEST_CASE("TWFFastDerivWrapper buildIntermediates without virtual orbitals", "[wavefunction][.wip]")
{
  using ValueMatrix = TWFFastDerivWrapper::ValueMatrix;

  const int n = 2;
  std::vector<ValueMatrix> Minv(1, ValueMatrix(n, n)), B(1, ValueMatrix(n, n)), M(1, ValueMatrix(n, n)),
      X(1, ValueMatrix(n, n)), Minv_B(1, ValueMatrix(n, n)), Minv_Mv(1, ValueMatrix(n, 0));
  Minv[0] = B[0] = M[0] = 1.0;

  TWFFastDerivWrapper twf;
  twf.buildIntermediates(Minv, B, M, X, Minv_B, Minv_Mv);
}
} // namespace qmcplusplus
