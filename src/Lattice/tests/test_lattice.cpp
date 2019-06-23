//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include <stdio.h>
#include <string>

#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/TinyVector.h"
#include "Lattice/CrystalLattice.h"
#include "Lattice/ParticleBConds.h"
#include "Configuration.h"

using std::string;

namespace qmcplusplus
{
typedef TinyVector<double, 3> vec_t;

/** Lattice is defined but Open BC is also used.
 */
TEST_CASE("Crystal_lattice_periodic_bulk", "[lattice]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> Lattice;
  Lattice.BoxBConds = false; // Open BC
  Lattice.R.diagonal(0.4);
  Lattice.reset();

  REQUIRE(Lattice.Volume == Approx(0.4 * 0.4 * 0.4));

  vec_t v3(0.6, 1.2, -1.7);
  REQUIRE(Lattice.isValid(v3) == false);
  REQUIRE(Lattice.outOfBound(v3) == true);

  vec_t v4(0.45, 0.2, 0.1);
  REQUIRE(Lattice.isValid(v4) == true);
  REQUIRE(Lattice.outOfBound(v4) == false);
}

} // namespace qmcplusplus
