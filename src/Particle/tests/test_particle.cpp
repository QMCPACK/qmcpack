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
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Particle/SymmetricDistanceTableData.h"



#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{


TEST_CASE("symmetric_distance_table", "[particle]")
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

  source.addTable(source,DT_SOA_PREFERRED);
  source.update();

  DistanceTableData *dist2 = createDistanceTable(source,DT_AOS);
}

TEST_CASE("particle set lattice with vacuum", "[particle]")
{

  OHMMS::Controller->initialize(0, NULL);

  typedef SymmetricDTD<double, 3, SUPERCELL_BULK> sym_dtd_t;
  ParticleSet source;

  Uniform3DGridLayout grid;
  // PPP case
  grid.BoxBConds = true;
  grid.R(0) = 1.0;
  grid.R(1) = 2.0;
  grid.R(2) = 3.0;

  grid.R(3) = 0.0;
  grid.R(4) = 1.0;
  grid.R(5) = 0.0;

  grid.R(6) = 0.0;
  grid.R(7) = 0.0;
  grid.R(8) = 1.0;

  grid.VacuumScale=2.0;
  grid.reset();

  source.setName("electrons");
  source.Lattice.copy(grid);
  source.createSK();

  REQUIRE( source.LRBox.R(0,0) == 1.0 );
  REQUIRE( source.LRBox.R(0,1) == 2.0 );
  REQUIRE( source.LRBox.R(0,2) == 3.0 );

  // PPN case
  grid.BoxBConds[2] = false;
  grid.reset();
  source.Lattice.copy(grid);
  source.createSK();

  REQUIRE( source.LRBox.R(2,0) == 0.0 );
  REQUIRE( source.LRBox.R(2,1) == 0.0 );
  REQUIRE( source.LRBox.R(2,2) == 2.0 );

  // PNN case
  grid.BoxBConds[1] = false;
  grid.reset();
  source.Lattice.copy(grid);
  source.createSK();

  REQUIRE( source.LRBox.R(0,0) ==  1.0 );
  REQUIRE( source.LRBox.R(0,1) ==  2.0 );
  REQUIRE( source.LRBox.R(0,2) ==  3.0 );
  REQUIRE( source.LRBox.R(1,0) ==  0.0 );
  REQUIRE( source.LRBox.R(1,1) ==  2.0 );
  REQUIRE( source.LRBox.R(1,2) ==  0.0 );
  REQUIRE( source.LRBox.R(2,0) ==  0.0 );
  REQUIRE( source.LRBox.R(2,1) ==  0.0 );
  REQUIRE( source.LRBox.R(2,2) ==  2.0 );
}

}
