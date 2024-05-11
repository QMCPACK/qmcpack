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

TEST_CASE("solveGeneralizedEigenvalues", "[drivers]")
{
  // Eigenvalues and eigenvectors from gen_eigenval.py
  LinearMethod lm;
  const int N = 2;
  Matrix<Real> Ovlp(N, N);
  Matrix<Real> Ham(N, N);

  Ovlp(0, 0) = 1.0;
  Ovlp(0, 1) = 0.1;
  Ovlp(1, 0) = 0.15;
  Ovlp(1, 1) = 1.0;

  Ham(0, 0) = -4.0;
  Ham(0, 1) = 0.2;
  Ham(1, 0) = 0.3;
  Ham(1, 1) = 1.0;
  std::vector<Real> ev(N);
  Matrix<Real> evec(N, N);
  lm.solveGeneralizedEigenvalues(Ham, Ovlp, ev, evec);
  CHECK(ev[0] == Approx(-4.10958));
  CHECK(ev[1] == Approx(1.00298));
  app_log() << "ev = " << ev[0] << " " << ev[1] << std::endl;

  app_log() << " evec " << evec(0, 0) << " " << evec(0, 1) << std::endl;
  app_log() << " evec " << evec(1, 0) << " " << evec(1, 1) << std::endl;

  CHECK(evec(0, 0) == Approx(-1.0));
  CHECK(evec(0, 1) == Approx(0.17935662));
  CHECK(evec(1, 0) == Approx(0.01992851));
  CHECK(evec(1, 1) == Approx(1.0));
}

TEST_CASE("solveGeneralizedEigenvaluesInv", "[drivers]")
{
  LinearMethod lm;
  const int N = 2;
  Matrix<Real> Ovlp(N, N);
  Matrix<Real> Ham(N, N);

  Ovlp(0, 0) = 1.0;
  Ovlp(0, 1) = 0.1;
  Ovlp(1, 0) = 0.15;
  Ovlp(1, 1) = 1.0;

  Ham(0, 0) = -4.0;
  Ham(0, 1) = 0.2;
  Ham(1, 0) = 0.3;
  Ham(1, 1) = 1.0;
  std::vector<Real> ev(N);
  Matrix<Real> evec(N, N);
  lm.solveGeneralizedEigenvalues_Inv(Ham, Ovlp, ev, evec);
  CHECK(ev[0] == Approx(-4.10958));
  CHECK(ev[1] == Approx(1.00298));
  app_log() << "ev = " << ev[0] << " " << ev[1] << std::endl;

  app_log() << " evec " << evec(0, 0) << " " << evec(0, 1) << std::endl;
  app_log() << " evec " << evec(1, 0) << " " << evec(1, 1) << std::endl;

  CHECK(evec(0, 0) == Approx(-0.98429354));
  CHECK(evec(0, 1) == Approx(0.17653957));
  CHECK(evec(1, 0) == Approx(-0.01992456));
  CHECK(evec(1, 1) == Approx(-0.99980149));
}


TEST_CASE("solveGeneralizedEigenvaluesCompare", "[drivers]")
{
  LinearMethod lm;
  const int N = 4;
  Matrix<Real> Ovlp(N, N);
  Matrix<Real> Ham(N, N);

  std::mt19937 mt(100);
  std::uniform_real_distribution<Real> rnd(0.0, 1.0);


  for (int i = 0; i < N; i++)
  {
    Ovlp(i, i) = 1.0;
  }

  for (int i = 0; i < N; i++)
  {
    Ovlp(i, i) = 1.0;
    for (int j = 0; j < i; j++)
    {
      Ovlp(i, j) = rnd(mt);
      Ovlp(j, i) = Ovlp(i, j) + 0.1;
    }
  }


  for (int i = 0; i < N; i++)
  {
    Ham(i, i) = 1.0;
    for (int j = 0; j < i; j++)
    {
      Ham(i, j) = rnd(mt);
      Ham(j, i) = Ham(i, j) - 0.1;
    }
  }

  Matrix<Real> Ovlp_copy(N, N);
  Matrix<Real> Ham_copy(N, N);

  Ovlp_copy = Ovlp;
  Ham_copy  = Ham;

  std::vector<Real> ev1(N);
  Matrix<Real> evec1(N, N);
  lm.solveGeneralizedEigenvalues(Ham_copy, Ovlp_copy, ev1, evec1);

  Ovlp_copy = Ovlp;
  Ham_copy  = Ham;

  std::vector<Real> ev2(N);
  Matrix<Real> evec2(N, N);
  lm.solveGeneralizedEigenvalues_Inv(Ham_copy, Ovlp_copy, ev2, evec2);

  for (int i = 0; i < N; i++)
  {
    CHECK(ev1[i] == Approx(ev2[i]));
  }
}

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
