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
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/LatticeGaussianProduct.h"
#include "QMCWaveFunctions/LatticeGaussianProductBuilder.h"
#include "ParticleIO/ParticleLayoutIO.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
using RealType = QMCTraits::RealType;
TEST_CASE("lattice gaussian", "[wavefunction]")
{
  Communicate* c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  ParticleSet ions_;
  ParticleSet elec_;

  ions_.setName("ion");
  ions_.create(2);
  ions_.R[0][0] = 0.0;
  ions_.R[0][1] = 0.0;
  ions_.R[0][2] = 0.0;
  ions_.R[1][0] = 0.0;
  ions_.R[1][1] = 0.0;
  ions_.R[1][2] = 0.0;

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

  std::map<string, ParticleSet*> pp;
  pp["ion"] = &ions_;

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

  // read lattice
  xmlNodePtr root  = doc.getRoot();
  xmlNodePtr part1 = xmlFirstElementChild(root);
  ParticleSet::ParticleLayout_t* SimulationCell = new ParticleSet::ParticleLayout_t;
  LatticeParser lp(*SimulationCell);
  lp.put(part1);
  SimulationCell->print(app_log(), 0);
  elec_.Lattice = *SimulationCell;
  // initialize SK
  elec_.createSK();

  const char* particles = "<tmp> \
  <ionwf name=\"ionwf\" source=\"ion\" width=\"0.5 0.5\"/> \
</tmp> \
";
  okay = doc.parseFromString(particles);
  REQUIRE(okay);

  root = doc.getRoot();

  xmlNodePtr jas1 = xmlFirstElementChild(root);

  LatticeGaussianProductBuilder jastrow(c, elec_, pp);
  std::unique_ptr<LatticeGaussianProduct> LGP(dynamic_cast<LatticeGaussianProduct*>(jastrow.buildComponent(jas1)));
  double width = 0.5;
  double alpha = 1./(2*width*width);
  // check initialization. Nope, cannot access Psi.Z
  for (int i=0; i<2; i++)
  {
    REQUIRE(LGP->ParticleAlpha[i] == Approx(alpha));
  }

  // update all distance tables
  ions_.update();
  elec_.update();

  RealType logpsi = LGP->evaluateLog(elec_, elec_.G, elec_.L);
  // check answer
  RealType r2 = Dot(elec_.R, elec_.R);
  double wfval = std::exp(-alpha*r2);
  REQUIRE(logpsi == Approx(std::log(wfval)));
}
} // namespace qmcplusplus
