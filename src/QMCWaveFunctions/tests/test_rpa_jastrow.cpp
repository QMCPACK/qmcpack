//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Lattice/ParticleBConds.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "Particle/SymmetricDistanceTableData.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/Jastrow/RPAJastrow.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "ParticleIO/ParticleLayoutIO.h"


#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("RPA Jastrow", "[wavefunction]")
{
  Communicate* c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  ParticleSet ions_;
  ParticleSet elec_;

  ions_.setName("ion");
  ions_.create(2);
  ions_.R[0][0] = 2.0;
  ions_.R[0][1] = 0.0;
  ions_.R[0][2] = 0.0;
  ions_.R[1][0] = -2.0;
  ions_.R[1][1] = 0.0;
  ions_.R[1][2] = 0.0;
  SpeciesSet& source_species(ions_.getSpeciesSet());
  source_species.addSpecies("O");
  ions_.RSoA = ions_.R;
  //ions_.resetGroups();

  elec_.setName("elec");
  std::vector<int> ud(2);
  ud[0] = ud[1] = 2;
  elec_.create(ud);
  elec_.R[0][0] = 1.00;
  elec_.R[0][1] = 0.0;
  elec_.R[0][2] = 0.0;
  elec_.R[1][0] = 0.0;
  elec_.R[1][1] = 0.0;
  elec_.R[1][2] = 0.0;
  elec_.R[2][0] = -1.00;
  elec_.R[2][1] = 0.0;
  elec_.R[2][2] = 0.0;
  elec_.R[3][0] = 0.0;
  elec_.R[3][1] = 0.0;
  elec_.R[3][2] = 2.0;

  SpeciesSet& target_species(elec_.getSpeciesSet());
  int upIdx                          = target_species.addSpecies("u");
  int downIdx                        = target_species.addSpecies("d");
  int chargeIdx                      = target_species.addAttribute("charge");
  target_species(chargeIdx, upIdx)   = -1;
  target_species(chargeIdx, downIdx) = -1;

  // initialize simulationcell for kvectors
  const char* xmltext = "<tmp> \
  <simulationcell>\
     <parameter name=\"lattice\" units=\"bohr\">\
              6.00000000        0.00000000        0.00000000\
              0.00000000        6.00000000        0.00000000\
              0.00000000        0.00000000        6.00000000\
     </parameter>\
     <parameter name=\"bconds\">\
        p p p\
     </parameter>\
     <parameter name=\"LR_dim_cutoff\"       >    15                 </parameter>\
  </simulationcell>\
</tmp> ";
  Libxml2Document doc;
  bool okay = doc.parseFromString(xmltext);
  REQUIRE(okay);

  xmlNodePtr root  = doc.getRoot();
  xmlNodePtr part1 = xmlFirstElementChild(root);

  // read lattice
  ParticleSet::ParticleLayout_t* SimulationCell = new ParticleSet::ParticleLayout_t;
  LatticeParser lp(*SimulationCell);
  lp.put(part1);
  SimulationCell->print(app_log(), 0);
  elec_.Lattice = *SimulationCell;
  // initialize SK
  elec_.createSK();

  TrialWaveFunction psi = TrialWaveFunction(c);

  xmltext = "<tmp> \
  <jastrow name=\"Jee\" type=\"Two-Body\" function=\"rpa\"/>\
</tmp> ";
  okay = doc.parseFromString(xmltext);
  REQUIRE(okay);

  root = doc.getRoot();

  xmlNodePtr jas_node = xmlFirstElementChild(root);
  bool is_manager=false;
  RPAJastrow* jas = new RPAJastrow(elec_, is_manager);
  jas->put(root);

  psi.addComponent(jas, "Jee");

  // update all distance tables
  elec_.update();

  double logpsi = psi.evaluateLog(elec_);
  REQUIRE(logpsi == Approx(-1.3327837613)); // note: number not validated
}
} // namespace qmcplusplus
