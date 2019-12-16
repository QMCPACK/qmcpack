//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Yubo "Paul Yang", yubo.paul.yang@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Yubo "Paul Yang", yubo.paul.yang@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Lattice/ParticleBConds.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Jastrow/kSpaceJastrow.h"
#include "QMCWaveFunctions/Jastrow/kSpaceJastrowBuilder.h"
#include "ParticleIO/ParticleLayoutIO.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("kspace jastrow", "[wavefunction]")
{
  Communicate* c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  ParticleSet ions_;
  ParticleSet elec_;

  ions_.setName("ion");
  ions_.create(1);
  ions_.R[0][0] = 0.0;
  ions_.R[0][1] = 0.0;
  ions_.R[0][2] = 0.0;

  elec_.setName("elec");
  std::vector<int> ud(2);
  ud[0] = 2;
  ud[1] = 0;
  elec_.create(ud);
  elec_.R[0][0] = -0.28;
  elec_.R[0][1] = 0.0225;
  elec_.R[0][2] = -2.709;
  elec_.R[1][0] = -1.08389;
  elec_.R[1][1] = 1.9679;
  elec_.R[1][2] = -0.0128914;

  SpeciesSet& tspecies         = elec_.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;
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
  // kSpaceJastrow::setupGvecs requires ions to have lattice
  ions_.Lattice = *SimulationCell;
  // initialize SK
  elec_.createSK();

  const char* particles = "<tmp> \
<jastrow name=\"Jk\" type=\"kSpace\" source=\"ion\"> \
  <correlation kc=\"1.5\" type=\"Two-Body\" symmetry=\"isotropic\"> \
    <coefficients id=\"cG2\" type=\"Array\"> \
      -100. -50. \
    </coefficients> \
  </correlation> \
</jastrow> \
</tmp> \
";
  okay = doc.parseFromString(particles);
  REQUIRE(okay);

  root = doc.getRoot();
  xmlNodePtr jas1 = xmlFirstElementChild(root);

  kSpaceJastrowBuilder jastrow(c, elec_, ions_);
  std::unique_ptr<WaveFunctionComponent> jas(jastrow.buildComponent(jas1));

  // update all distance tables
  elec_.update();

  double logpsi_real = std::real(jas->evaluateLog(elec_, elec_.G, elec_.L));
  REQUIRE(logpsi_real == Approx(-4.4088303951)); // !!!! value not checked
}
} // namespace qmcplusplus
