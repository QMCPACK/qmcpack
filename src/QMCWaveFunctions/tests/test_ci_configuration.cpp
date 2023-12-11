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

#include "Fermion/ci_configuration2.h"

namespace qmcplusplus
{
TEST_CASE("ci_configuration2", "[wavefunction]")
{
  const int n_states = 6;
  std::vector<size_t> ref{0, 1, 2, 3};
  ci_configuration2 ref_state(ref);

  size_t n_excited;
  double sign = 0.0;
  std::vector<size_t> pos(n_states);
  std::vector<size_t> occupied(n_states);
  std::vector<size_t> unoccupied(n_states);

  std::vector<size_t> ext1{0, 1, 4, 5};
  ci_configuration2 ext1_state(ext1);

  n_excited = ext1_state.calculateNumOfExcitations(ref_state);
  REQUIRE(n_excited == 2);

  sign = ext1_state.calculateExcitations(ref_state, n_excited, pos, occupied, unoccupied);

  REQUIRE(n_excited == 2);
  CHECK(sign == 1.0);
  CHECK(pos[0] == 2);
  CHECK(pos[1] == 3);
  CHECK(occupied[0] == 4);
  CHECK(occupied[1] == 5);
  CHECK(unoccupied[0] == 2);
  CHECK(unoccupied[1] == 3);

  std::vector<size_t> ext2{0, 1, 2, 6};
  ci_configuration2 ext2_state(ext2);

  n_excited = ext2_state.calculateNumOfExcitations(ref_state);
  REQUIRE(n_excited == 1);

  sign = ext2_state.calculateExcitations(ref_state, n_excited, pos, occupied, unoccupied);

  REQUIRE(n_excited == 1);
  CHECK(sign == 1.0);
  CHECK(pos[0] == 3);
  CHECK(occupied[0] == 6);
  CHECK(unoccupied[0] == 3);

  std::vector<size_t> ref2{0, 1};
  ci_configuration2 ref2_state(ref2);
  std::vector<size_t> ext3{1, 6};
  ci_configuration2 ext3_state(ext3);

  sign = ext3_state.calculateExcitations(ref2_state, n_excited, pos, occupied, unoccupied);
  REQUIRE(n_excited == 1);
  CHECK(sign == -1.0);

  sign = ref2_state.calculateExcitations(ext3_state, n_excited, pos, occupied, unoccupied);
  REQUIRE(n_excited == 1);
  CHECK(sign == -1.0);
}
} // namespace qmcplusplus
