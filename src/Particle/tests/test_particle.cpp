
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
  OhmmsInfo("testlogfile");

  typedef SymmetricDTD<double, 3, SUPERCELL_OPEN> sym_dtd_t;
  ParticleSet source;

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
