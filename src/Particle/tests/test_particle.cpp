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


#include "Message/catch_mpi_main.hpp"

#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/TinyVector.h"
#include "Lattice/CrystalLattice.h"
#include "Lattice/ParticleBConds.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "Particle/SymmetricDistanceTableData.h"


#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("symmetric_distance_table OpenBC", "[particle]")
{
  OHMMS::Controller->initialize(0, NULL);

  ParticleSet source;

  source.setName("electrons");

  source.create(2);
  source.R[0][0] = 0.0;
  source.R[0][1] = 1.0;
  source.R[0][2] = 2.0;
  source.R[1][0] = 1.1;
  source.R[1][1] = 1.0;
  source.R[1][2] = 3.2;

  const int TableID = source.addTable(source, DT_SOA);
  source.update();
  const auto& d_aa = source.getDistTable(TableID);

  REQUIRE(d_aa.Distances[0][1] == Approx(1.62788206));
  REQUIRE(d_aa.Distances[1][0] == Approx(1.62788206));
  REQUIRE(d_aa.Displacements[0][1][0] == Approx(1.1));
  REQUIRE(d_aa.Displacements[0][1][1] == Approx(0.0));
  REQUIRE(d_aa.Displacements[0][1][2] == Approx(1.2));
  REQUIRE(d_aa.Displacements[1][0][0] == Approx(-1.1));
  REQUIRE(d_aa.Displacements[1][0][1] == Approx(0.0));
  REQUIRE(d_aa.Displacements[1][0][2] == Approx(-1.2));
}

TEST_CASE("symmetric_distance_table PBC", "[particle]")
{
  OHMMS::Controller->initialize(0, NULL);

  ParticleSet source;

  Uniform3DGridLayout grid;
  grid.BoxBConds = true; // periodic
  grid.R = ParticleSet::Tensor_t(6.74632230, 6.74632230, 0.00000000, 0.00000000, 3.37316115, 3.37316115, 3.37316115,
                                 0.00000000, 3.37316115);
  grid.reset();

  source.setName("electrons");
  source.Lattice.copy(grid);

  source.create(4);
  source.R[0] = ParticleSet::PosType(0.00000000, 0.00000000, 0.00000000);
  source.R[1] = ParticleSet::PosType(1.68658058, 1.68658058, 1.68658058);
  source.R[2] = ParticleSet::PosType(3.37316115, 3.37316115, 0.00000000);
  source.R[3] = ParticleSet::PosType(5.05974172, 5.05974172, 1.68658058);

  const int TableID = source.addTable(source, DT_SOA);
  source.update();
  const auto& d_aa = source.getDistTable(TableID);

  REQUIRE(d_aa.Distances[1][2] == Approx(2.9212432441));
  REQUIRE(d_aa.Distances[2][1] == Approx(2.9212432441));
  REQUIRE(d_aa.Displacements[1][2][0] == Approx(1.68658057));
  REQUIRE(d_aa.Displacements[1][2][1] == Approx(1.68658057));
  REQUIRE(d_aa.Displacements[1][2][2] == Approx(-1.68658058));
  REQUIRE(d_aa.Displacements[2][1][0] == Approx(-1.68658057));
  REQUIRE(d_aa.Displacements[2][1][1] == Approx(-1.68658057));
  REQUIRE(d_aa.Displacements[2][1][2] == Approx(1.68658057));
}

TEST_CASE("particle set lattice with vacuum", "[particle]")
{
  OHMMS::Controller->initialize(0, NULL);

  typedef SymmetricDTD<double, 3, SUPERCELL_BULK> sym_dtd_t;
  ParticleSet source;

  Uniform3DGridLayout grid;
  // PPP case
  grid.BoxBConds = true;
  grid.R(0)      = 1.0;
  grid.R(1)      = 2.0;
  grid.R(2)      = 3.0;

  grid.R(3) = 0.0;
  grid.R(4) = 1.0;
  grid.R(5) = 0.0;

  grid.R(6) = 0.0;
  grid.R(7) = 0.0;
  grid.R(8) = 1.0;

  grid.VacuumScale = 2.0;
  grid.reset();

  source.setName("electrons");
  source.Lattice.copy(grid);
  source.createSK();

  REQUIRE(source.LRBox.R(0, 0) == 1.0);
  REQUIRE(source.LRBox.R(0, 1) == 2.0);
  REQUIRE(source.LRBox.R(0, 2) == 3.0);

  // PPN case
  grid.BoxBConds[2] = false;
  grid.reset();
  source.Lattice.copy(grid);
  source.createSK();

  REQUIRE(source.LRBox.R(2, 0) == 0.0);
  REQUIRE(source.LRBox.R(2, 1) == 0.0);
  REQUIRE(source.LRBox.R(2, 2) == 2.0);

  // PNN case
  grid.BoxBConds[1] = false;
  grid.reset();
  source.Lattice.copy(grid);
  source.createSK();

  REQUIRE(source.LRBox.R(0, 0) == 1.0);
  REQUIRE(source.LRBox.R(0, 1) == 2.0);
  REQUIRE(source.LRBox.R(0, 2) == 3.0);
  REQUIRE(source.LRBox.R(1, 0) == 0.0);
  REQUIRE(source.LRBox.R(1, 1) == 2.0);
  REQUIRE(source.LRBox.R(1, 2) == 0.0);
  REQUIRE(source.LRBox.R(2, 0) == 0.0);
  REQUIRE(source.LRBox.R(2, 1) == 0.0);
  REQUIRE(source.LRBox.R(2, 2) == 2.0);
}

} // namespace qmcplusplus
