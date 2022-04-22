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
#include "Particle/ParticleSet.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/LatticeGaussianProduct.h"
#include "QMCWaveFunctions/LatticeGaussianProductBuilder.h"
#include "ParticleIO/LatticeIO.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
using RealType     = QMCTraits::RealType;
using LogValueType = std::complex<QMCTraits::QTFull::RealType>;

TEST_CASE("lattice gaussian", "[wavefunction]")
{
  Communicate* c;
  c = OHMMS::Controller;

  ParticleSet::ParticleLayout lattice;
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

  LatticeParser lp(lattice);
  lp.put(part1);
  lattice.print(app_log(), 0);
  const SimulationCell simulation_cell(lattice);
  auto ions_ptr = std::make_unique<ParticleSet>(simulation_cell);
  auto elec_ptr = std::make_unique<ParticleSet>(simulation_cell);
  auto &ions(*ions_ptr), elec(*elec_ptr);

  ions.setName("ion");
  ions.create({2});
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.R[1][0] = 0.0;
  ions.R[1][1] = 0.0;
  ions.R[1][2] = 0.0;

  elec.setName("elec");
  elec.create({2, 0});
  elec.R[0][0] = -0.28;
  elec.R[0][1] = 0.0225;
  elec.R[0][2] = -2.709;
  elec.R[1][0] = -1.08389;
  elec.R[1][1] = 1.9679;
  elec.R[1][2] = -0.0128914;

  std::map<string, const std::unique_ptr<ParticleSet>> pp;
  pp.emplace(ions_ptr->getName(), std::move(ions_ptr));
  pp.emplace(elec_ptr->getName(), std::move(elec_ptr));

  SpeciesSet& tspecies         = elec.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;

  // initialize SK
  elec.createSK();

  const char* particles = "<tmp> \
  <ionwf name=\"ionwf\" source=\"ion\" width=\"0.5 0.5\"/> \
</tmp> \
";
  okay                  = doc.parseFromString(particles);
  REQUIRE(okay);

  root = doc.getRoot();

  xmlNodePtr jas1 = xmlFirstElementChild(root);

  LatticeGaussianProductBuilder jastrow(c, elec, pp);
  auto LGP_uptr = jastrow.buildComponent(jas1);
  auto LGP      = dynamic_cast<LatticeGaussianProduct*>(LGP_uptr.get());
  double width  = 0.5;
  double alpha  = 1. / (2 * width * width);
  // check initialization. Nope, cannot access Psi.Z
  for (int i = 0; i < 2; i++)
  {
    REQUIRE(LGP->ParticleAlpha[i] == Approx(alpha));
  }

  // update all distance tables
  ions.update();
  elec.update();

  LogValueType logpsi = LGP->evaluateLog(elec, elec.G, elec.L);
  // check answer
  RealType r2  = Dot(elec.R, elec.R);
  double wfval = std::exp(-alpha * r2);
  REQUIRE(logpsi == ComplexApprox(std::log(wfval)));
}
} // namespace qmcplusplus
