//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include <stdio.h>
#include <string>

#include "catch.hpp"

#include "Configuration.h"
#include "OhmmsPETE/TinyVector.h"
#include "Lattice/CrystalLattice.h"

namespace qmcplusplus
{
using vec_t = TinyVector<double, 3>;

/** Lattice is defined but Open BC is also used.
 */
TEST_CASE("Crystal_lattice_periodic_bulk", "[lattice]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> Lattice;
  Lattice.BoxBConds = false; // Open BC
  Lattice.R.diagonal(0.4);
  Lattice.reset();

  CHECK(Lattice.Volume == Approx(0.4 * 0.4 * 0.4));

  vec_t v3(0.6, 1.2, -1.7);
  REQUIRE(Lattice.isValid(v3) == false);
  REQUIRE(Lattice.outOfBound(v3) == true);

  vec_t v4(0.45, 0.2, 0.1);
  REQUIRE(Lattice.isValid(v4) == true);
  REQUIRE(Lattice.outOfBound(v4) == false);
}

} // namespace qmcplusplus
