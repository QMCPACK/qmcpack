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


TEST_CASE("solveGeneralizedEigenvalues", "[drivers]")
{
  // Eigenvalues and eigenvectors from gen_eigenval.py
  LinearMethod lm;
  const int N = 2;
  Matrix<double> Ovlp(N, N);
  Matrix<double> Ham(N, N);

  Ovlp(0, 0) = 1.0;
  Ovlp(0, 1) = 0.1;
  Ovlp(1, 0) = 0.15;
  Ovlp(1, 1) = 1.0;

  Ham(0, 0) = -4.0;
  Ham(0, 1) = 0.2;
  Ham(1, 0) = 0.3;
  Ham(1, 1) = 1.0;
  std::vector<double> ev(N);
  Matrix<double> evec(N, N);
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
  Matrix<double> Ovlp(N, N);
  Matrix<double> Ham(N, N);

  Ovlp(0, 0) = 1.0;
  Ovlp(0, 1) = 0.1;
  Ovlp(1, 0) = 0.15;
  Ovlp(1, 1) = 1.0;

  Ham(0, 0) = -4.0;
  Ham(0, 1) = 0.2;
  Ham(1, 0) = 0.3;
  Ham(1, 1) = 1.0;
  std::vector<double> ev(N);
  Matrix<double> evec(N, N);
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
  Matrix<double> Ovlp(N, N);
  Matrix<double> Ham(N, N);

  std::mt19937 mt(100);
  std::uniform_real_distribution<double> rnd(0.0, 1.0);


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

  Matrix<double> Ovlp_copy(N, N);
  Matrix<double> Ham_copy(N, N);

  Ovlp_copy = Ovlp;
  Ham_copy  = Ham;

  std::vector<double> ev1(N);
  Matrix<double> evec1(N, N);
  lm.solveGeneralizedEigenvalues(Ham_copy, Ovlp_copy, ev1, evec1);

  Ovlp_copy = Ovlp;
  Ham_copy  = Ham;

  std::vector<double> ev2(N);
  Matrix<double> evec2(N, N);
  lm.solveGeneralizedEigenvalues_Inv(Ham_copy, Ovlp_copy, ev2, evec2);

  for (int i = 0; i < N; i++)
  {
    CHECK(ev1[i] == Approx(ev2[i]));
  }
}

} // namespace qmcplusplus
