//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "Message/catch_mpi_main.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Lattice/ParticleBConds.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "Particle/SymmetricDistanceTableData.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCHamiltonians/BareKineticEnergy.h"


#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{

TEST_CASE("Bare Kinetic Energy", "[hamiltonian]")
{

  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  ParticleSet ions;
  ParticleSet elec;

  ions.setName("ion");
  ions.create(1);
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;

  elec.setName("elec");
  elec.create(2);
  elec.R[0][0] = 0.0;
  elec.R[0][1] = 1.0;
  elec.R[0][2] = 0.0;
  elec.R[1][0] = 1.0;
  elec.R[1][1] = 1.0;
  elec.R[1][2] = 0.0;

  SpeciesSet &tspecies =  elec.getSpeciesSet();
  int upIdx = tspecies.addSpecies("u");
  int downIdx = tspecies.addSpecies("d");
  //int chargeIdx = tspecies.addAttribute("charge");
  int massIdx = tspecies.addAttribute("mass");
  //tspecies(chargeIdx, upIdx) = -1;
  //tspecies(chargeIdx, downIdx) = -1;
  tspecies(massIdx, upIdx) = 1.0;
  tspecies(massIdx, downIdx) = 1.0;

  elec.addTable(ions,DT_AOS);
  elec.update();


  ParticleSetPool ptcl = ParticleSetPool(c);


  const char *particles = \
"<tmp> \
</tmp> \
";
  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr h1 = xmlFirstElementChild(root);

  // This constructor has compile errors (Ps member not initialized)
  //BareKineticEnergy<double> bare_ke;

  BareKineticEnergy<double> bare_ke(elec);
  bare_ke.put(h1);

  elec.L[0] = 1.0;
  elec.L[1] = 0.0;
  double v = bare_ke.evaluate(elec);
  REQUIRE(v == -0.5);

  elec.L[0] = 0.0;
  elec.L[1] = 0.0;
  elec.G[0][0] = 1.0;
  elec.G[0][1] = 0.0;
  elec.G[0][2] = 0.0;
  elec.G[1][0] = 0.0;
  elec.G[1][1] = 0.0;
  elec.G[1][2] = 0.0;
  v = bare_ke.evaluate(elec);
  REQUIRE(v == -0.5);
}
}

