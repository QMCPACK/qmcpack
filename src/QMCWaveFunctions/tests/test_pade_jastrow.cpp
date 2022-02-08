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


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Jastrow/PadeFunctors.h"
#include "QMCWaveFunctions/Jastrow/RadialJastrowBuilder.h"
#include "QMCWaveFunctions/WaveFunctionFactory.h"


#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("Pade functor", "[wavefunction]")
{
  double A = -0.25;
  double B = 0.1;
  PadeFunctor<double> pf;
  pf.B0 = B;
  pf.setCusp(A);

  double r = 1.2;
  double u = pf.evaluate(r);
  REQUIRE(u == Approx(2.232142857142857));
}

TEST_CASE("Pade2 functor", "[wavefunction]")
{
  double A = 0.8;
  double B = 5.0;
  double C = -0.1;
  Pade2ndOrderFunctor<double> pf2;
  pf2.A = A;
  pf2.B = B;
  pf2.C = C;
  pf2.reset();

  double r = 1.2;
  double u = pf2.evaluate(r);
  REQUIRE(u == Approx(0.11657142857142856));
}


TEST_CASE("Pade Jastrow", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  const SimulationCell simulation_cell;
  ParticleSet ions_(simulation_cell);
  ParticleSet elec_(simulation_cell);

  // Need 1 electron and 1 proton, somehow
  ions_.setName("ion");
  ions_.create({1});
  ions_.R[0][0] = 0.0;
  ions_.R[0][1] = 0.0;
  ions_.R[0][2] = 0.0;

  elec_.setName("elec");
  elec_.create({2,0});
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
  int massIdx                  = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;
  tspecies(massIdx, upIdx)     = 1;
  tspecies(massIdx, downIdx)   = 1;

  const char* particles = "<tmp> \
<jastrow name=\"Jee\" type=\"Two-Body\" function=\"pade\"> \
  <correlation speciesA=\"u\" speciesB=\"u\"> \
        <var id=\"juu_b\" name=\"B\">0.1</var> \
  </correlation> \
</jastrow> \
</tmp> \
";
  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr jas1 = xmlFirstElementChild(root);

  // cusp = -0.25
  // r_ee = 3.42050023755
  RadialJastrowBuilder jastrow(c, elec_);
  std::unique_ptr<WaveFunctionComponent> jas(jastrow.buildComponent(jas1));

  // update all distance tables
  elec_.update();

  double logpsi_real = std::real(jas->evaluateLog(elec_, elec_.G, elec_.L));
  REQUIRE(logpsi_real == Approx(-1.862821769493147));
}

TEST_CASE("Pade2 Jastrow", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  ParticleSetPool ptcl = ParticleSetPool(c);
  auto& simulation_cell(ptcl.getSimulationCell());
  auto ions_uptr = std::make_unique<ParticleSet>(simulation_cell);
  auto elec_uptr = std::make_unique<ParticleSet>(simulation_cell);
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion0");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({1});
  ions_.R[0]                 = {0.0, 0.0, 0.0};
  SpeciesSet& ispecies       = ions_.getSpeciesSet();
  int HIdx                   = ispecies.addSpecies("H");
  int ichargeIdx             = ispecies.addAttribute("charge");
  ispecies(ichargeIdx, HIdx) = 1.0;

  elec_.setName("e");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({1, 1});
  elec_.R[0] = {0.5, 0.5, 0.5};
  elec_.R[1] = {-0.5, -0.5, -0.5};

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int downIdx                = tspecies.addSpecies("d");
  int massIdx                = tspecies.addAttribute("mass");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(massIdx, upIdx)   = 1.0;
  tspecies(massIdx, downIdx) = 1.0;
  tspecies(chargeIdx, upIdx) = -1.0;
  tspecies(massIdx, downIdx) = -1.0;
  // Necessary to set mass
  elec_.resetGroups();

  ions_.update();
  elec_.addTable(elec_);
  elec_.addTable(ions_);
  elec_.update();

  const char* jasxml = "<wavefunction name=\"psi0\" target=\"e\"> \
  <jastrow name=\"J1\" type=\"One-Body\" function=\"pade2\" print=\"yes\" source=\"ion0\"> \
    <correlation elementType=\"H\"> \
        <var id=\"J1H_A\" name=\"A\">0.8</var> \
        <var id=\"J1H_B\" name=\"B\">5.0</var> \
        <var id=\"J1H_C\" name=\"C\">-0.1</var> \
    </correlation> \
  </jastrow> \
</wavefunction> \
";
  Libxml2Document doc;
  bool okay = doc.parseFromString(jasxml);
  REQUIRE(okay);

  xmlNodePtr jas1 = doc.getRoot();

  // update all distance tables
  elec_.update();
  WaveFunctionFactory wf_factory("psi0", elec_, ptcl.getPool(), c);
  wf_factory.put(jas1);
  auto& twf(*wf_factory.getTWF());
  twf.setMassTerm(elec_);
  twf.evaluateLog(elec_);
  twf.prepareGroup(elec_, 0);

  auto& twf_component_list = twf.getOrbitals();

  opt_variables_type active;
  twf.checkInVariables(active);
  active.removeInactive();
  int nparam = active.size_of_active();
  REQUIRE(nparam == 3);

  using ValueType = QMCTraits::ValueType;
  std::vector<ValueType> dlogpsi(nparam);
  std::vector<ValueType> dhpsioverpsi(nparam);
  //twf.evaluateDerivatives(elec_, active, dlogpsi, dhpsioverpsi);
  twf_component_list[0]->evaluateDerivatives(elec_, active, dlogpsi, dhpsioverpsi);

  // Numbers not validated
  std::vector<ValueType> expected_dlogpsi      = {-0.3249548841, 0.0376658437, -0.2814191847};
  std::vector<ValueType> expected_dhpsioverpsi = {0.0146266746, 0.0031788682, 0.4554097531};
  for (int i = 0; i < nparam; i++)
  {
    CHECK(dlogpsi[i] == ValueApprox(expected_dlogpsi[i]));
    CHECK(dhpsioverpsi[i] == ValueApprox(expected_dhpsioverpsi[i]));
  }
}
} // namespace qmcplusplus
