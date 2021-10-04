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
  Txy.emplace_back(1, -0.7, PosType(0.2, 0.3, 0.1));
  Txy.emplace_back(2, -0.3, PosType(0.3, 0.1, 0.2));
  Txy.emplace_back(3, 0.0, PosType(0.3, 0.2, 0.1));

  auto Txy0 = Txy;
  auto select0 = t_op.selectMove(0.0, Txy0);
  CHECK(select0 == nullptr);

  auto Txy1 = Txy;
  auto select1 = t_op.selectMove(0.4, Txy1);
  CHECK(select1 == nullptr);

  auto Txy2 = Txy;
  auto select2 = t_op.selectMove(0.5, Txy2);
  REQUIRE(select2 != nullptr);
  CHECK(select2->PID == 1);

  auto Txy3 = Txy;
  auto select3 = t_op.selectMove(0.6, Txy3);
  REQUIRE(select3 != nullptr);
  CHECK(select3->PID == 1);

  auto Txy4 = Txy;
  auto select4 = t_op.selectMove(0.85, Txy4);
  REQUIRE(select4 != nullptr);
  CHECK(select4->PID == 2);

  auto Txy5 = Txy;
  auto select5 = t_op.selectMove(0.9, Txy5);
  REQUIRE(select5 != nullptr);
  CHECK(select5->PID == 2);

  auto Txy6 = Txy;
  auto select6 = t_op.selectMove(float(1) - std::numeric_limits<float>::epsilon(), Txy6);
  REQUIRE(select6 != nullptr);
  CHECK(select6->PID == 2);
}
}
