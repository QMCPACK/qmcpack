
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

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


TEST_CASE("particle set open bc", "[particle]")
{

  OHMMS::Controller->initialize(0, NULL);
  OhmmsInfo("testlogfile");

  typedef SymmetricDTD<double, 3, SUPERCELL_OPEN> sym_dtd_t;
  ParticleSet source;

  source.setName("electrons");

  source.R.resize(2);
  source.R[0][0] = 0.0;
  source.R[0][1] = 1.0;
  source.R[0][2] = 2.0;
  source.R[1][0] = 1.1;
  source.R[1][1] = 1.0;
  source.R[1][2] = 3.2;

  sym_dtd_t dist(source, source);

  dist.evaluate(source);
  source.addTable(source);
  source.update();

  DistanceTableData *dist2 = createDistanceTable(source);
}

TEST_CASE("particle set lattice", "[particle]")
{

  OHMMS::Controller->initialize(0, NULL);
  OhmmsInfo("testlogfile");

  typedef SymmetricDTD<double, 3, SUPERCELL_BULK> sym_dtd_t;
  ParticleSet source;

  Uniform3DGridLayout grid;
  grid.BoxBConds = true; // periodic
  grid.R.diagonal(1.0);
  grid.reset();

  source.setName("electrons");
  source.Lattice.copy(grid);

  source.R.resize(2);
  source.R[0][0] = 0.0;
  source.R[0][1] = 1.0;
  source.R[0][2] = 2.0;
  source.R[1][0] = 1.1;
  source.R[1][1] = 1.0;
  source.R[1][2] = 3.2;

  sym_dtd_t dist(source, source);

  dist.evaluate(source);
  source.addTable(source);
  source.update();

  DistanceTableData *dist2 = createDistanceTable(source);
}

}
