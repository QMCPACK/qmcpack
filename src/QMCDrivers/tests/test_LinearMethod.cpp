//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com
//
// File created by: Mark Dewing, markdewing@gmail.com
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "QMCDrivers/WFOpt/LinearMethod.h"
#include "Utilities/RuntimeOptions.h"
#include <random>


namespace qmcplusplus
{

using Real = QMCTraits::RealType;

TEST_CASE("selectEigenvalues", "[drivers]")
{
  LinearMethod lm;
  const int Ne = 4; // Number of eigenvalues
  const int Nv = 4; // Size of eigenvectors

  // For direct solvers, Ne == Nv, but for iterative solvers we may have Ne <= Nv;

  Matrix<Real> evecs(Ne, Nv);
  evecs       = 1.0;
  evecs(0, 0) = 2.0;
  evecs(1, 0) = 3.0;
  std::vector<Real> evals{-1.2, -1.5}; // size Ne
  std::vector<Real> selected_evec(Nv);
  Real zerozero      = -1.0;
  Real selected_eval = lm.selectEigenvalue(evals, evecs, zerozero, selected_evec);

  // Selects the closest to zerozero - 2.0
  CHECK(selected_eval == Approx(-1.5));

  CHECK(selected_evec[0] == Approx(1.0));
  CHECK(selected_evec[1] == Approx(1.0 / 3.0));
}

} // namespace qmcplusplus
