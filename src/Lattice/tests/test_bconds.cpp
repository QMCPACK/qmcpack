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


#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/TinyVector.h"
#include "Lattice/CrystalLattice.h"
#include "Lattice/ParticleBConds.h"
#include "Lattice/Uniform3DGridLayout.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{

typedef TinyVector<double ,3> vec_t;

TEST_CASE("open_bconds", "[lattice]")
{

  typedef DTD_BConds<double, 3, SUPERCELL_OPEN> bcond_t;
  // Following give a compile error
  // test_bconds.cpp:27:23: error: request for member ‘apply_bc’ in ‘qmcplusplus::bcond’, which is of non-class type ‘qmcplusplus::DTD_BConds<double, 3u, 0>(qmcplusplus::CrystalLattice<double, 3u>)’
  //DTD_BConds<double, 3, SUPERCELL_OPEN> bcond(CrystalLattice<double, 3>);
  //double r2 = bcond.apply_bc(v);

  bcond_t *bcond = new bcond_t(CrystalLattice<double, 3>());

  vec_t  v(3.0, 4.0, 5.0);

  double r2 = bcond->apply_bc(v);
  REQUIRE(Approx(r2) == 50.0);


  std::vector<vec_t> disps(1);
  disps[0] = v;
  std::vector<double> r(1), rinv(1), rr(1);

  bcond->apply_bc(disps, r, rinv);

  REQUIRE(Approx(r[0]) == std::sqrt(50.0));
  REQUIRE(Approx(rinv[0]) == 1.0/std::sqrt(50.0));

  r[0] = 0.0;
  bcond->apply_bc(disps, r);
  REQUIRE(Approx(r[0]) == std::sqrt(50.0));

  bcond->evaluate_rsquared(disps.data(), rr.data(), disps.size());
  REQUIRE(Approx(rr[0]) == 50.0);
}

TEST_CASE("periodic_bulk_bconds", "[lattice]")
{

  typedef DTD_BConds<double, 3, SUPERCELL_BULK> bcond_t;

  CrystalLattice<double, 3>* cl = new CrystalLattice<double, 3>();
  std::vector<string> argv;
  argv.push_back("cubic");
  argv.push_back("0.4");
    

  cl->set(argv);

  REQUIRE(cl->Volume == Approx(0.4*0.4*0.4));


  bcond_t *bcond = new bcond_t(*cl);

  vec_t v1(0.0, 0.0, 0.0);

  double r2 = bcond->apply_bc(v1);
  REQUIRE(r2 == 0.0);

  vec_t v2(0.5, 0.0, 0.0);
  r2 = bcond->apply_bc(v2);
  REQUIRE(r2 == Approx(0.01));

}

TEST_CASE("uniform 3D grid layout", "[lattice]")
{
  Uniform3DGridLayout grid;
  grid.BoxBConds = true; // periodic

  grid.R.diagonal(1.0);
  grid.reset();
}

}
