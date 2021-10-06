//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"
#include <limits>
#include "QMCHamiltonians/NonLocalTOperator.h"

namespace qmcplusplus
{

TEST_CASE("NonLocalTOperator", "[hamiltonian]")
{
  using PosType = QMCTraits::PosType;
  NonLocalTOperator t_op;
  t_op.thingsThatShouldBeInMyConstructor("v0", 1.0, 0.0, 0.0);

  std::vector<NonLocalData> Txy;
  Txy.emplace_back(0, 0.4, PosType(0.1, 0.2, 0.3));
  Txy.emplace_back(1, -0.4, PosType(0.2, 0.3, 0.1));
  Txy.emplace_back(1, -0.3, PosType(0.2, 0.1, 0.3));
  Txy.emplace_back(2, -0.2, PosType(0.3, 0.1, 0.2));
  Txy.emplace_back(2, -0.1, PosType(0.3, 0.1, 0.2));
  Txy.emplace_back(3, 0.0, PosType(0.3, 0.2, 0.1));

  auto select0 = t_op.selectMove(0.0, Txy);
  CHECK(select0 == nullptr);

  auto select1 = t_op.selectMove(0.4, Txy);
  CHECK(select1 == nullptr);

  auto select2 = t_op.selectMove(0.5, Txy);
  REQUIRE(select2 != nullptr);
  CHECK(select2->PID == 1);

  auto select3 = t_op.selectMove(0.6, Txy);
  REQUIRE(select3 != nullptr);
  CHECK(select3->PID == 1);

  auto select4 = t_op.selectMove(0.85, Txy);
  REQUIRE(select4 != nullptr);
  CHECK(select4->PID == 2);

  auto select5 = t_op.selectMove(0.9, Txy);
  REQUIRE(select5 != nullptr);
  CHECK(select5->PID == 2);

  auto select6 = t_op.selectMove(float(1) - std::numeric_limits<float>::epsilon(), Txy);
  REQUIRE(select6 != nullptr);
  CHECK(select6->PID == 2);

  t_op.groupByElectron(4, Txy);

  auto select7 = t_op.selectMove(0.7, 1);
  REQUIRE(select7 != nullptr);
  CHECK(select7->Weight == Approx(-0.4));

  auto select8 = t_op.selectMove(0.7, 2);
  REQUIRE(select8 == nullptr);

  auto select9 = t_op.selectMove(0.8, 2);
  REQUIRE(select9 != nullptr);
  CHECK(select9->Weight == Approx(-0.2));
}
}
