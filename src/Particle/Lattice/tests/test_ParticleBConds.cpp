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
using vec_t = TinyVector<OHMMS_PRECISION, 3>;

TEST_CASE("open_bconds", "[lattice]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> Lattice;
  DTD_BConds<OHMMS_PRECISION, 3, SUPERCELL_OPEN> bcond(Lattice);

  vec_t v(3.0, 4.0, 5.0);

  OHMMS_PRECISION r2 = bcond.apply_bc(v);
  CHECK(Approx(r2) == 50.0);


  std::vector<vec_t> disps(1);
  disps[0] = v;
  std::vector<OHMMS_PRECISION> r(1), rinv(1), rr(1);

  bcond.apply_bc(disps, r, rinv);

  CHECK(Approx(r[0]) == std::sqrt(50.0));
  CHECK(Approx(rinv[0]) == 1.0 / std::sqrt(50.0));

  r[0] = 0.0;
  bcond.apply_bc(disps, r);
  CHECK(Approx(r[0]) == std::sqrt(50.0));

  bcond.evaluate_rsquared(disps.data(), rr.data(), disps.size());
  CHECK(Approx(rr[0]) == 50.0);
}

/** Lattice is defined but Open BC is also used.
 */
TEST_CASE("periodic_bulk_bconds", "[lattice]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> Lattice;
  Lattice.BoxBConds = false; // Open BC
  Lattice.R.diagonal(0.4);
  Lattice.reset();

  CHECK(Lattice.Volume == Approx(0.4 * 0.4 * 0.4));

  DTD_BConds<OHMMS_PRECISION, 3, SUPERCELL_BULK> bcond(Lattice);

  vec_t v1(0.0, 0.0, 0.0);

  OHMMS_PRECISION r2 = bcond.apply_bc(v1);
  REQUIRE(r2 == 0.0);

  vec_t v2(0.5, 0.0, 0.0);
  r2 = bcond.apply_bc(v2);
  CHECK(r2 == Approx(0.01));
}

TEST_CASE("uniform 3D Lattice layout", "[lattice]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> Lattice;
  Lattice.BoxBConds = true; // periodic

  Lattice.R.diagonal(1.0);
  Lattice.reset();

  REQUIRE(Lattice.R(0, 0) == 1.0);
  REQUIRE(Lattice.R(0, 1) == 0.0);
  REQUIRE(Lattice.R(0, 2) == 0.0);
  REQUIRE(Lattice.R(1, 0) == 0.0);
  REQUIRE(Lattice.R(1, 1) == 1.0);
  REQUIRE(Lattice.R(1, 2) == 0.0);
  REQUIRE(Lattice.R(2, 0) == 0.0);
  REQUIRE(Lattice.R(2, 1) == 0.0);
  REQUIRE(Lattice.R(2, 2) == 1.0);
}

} // namespace qmcplusplus
