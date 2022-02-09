//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Shiv Upadhyay, shivnupadhyay@gmail.com, University of Pittsburgh
//
// File created by: Shiv Upadhyay, shivnupadhyay@gmail.com, University of Pittsburgh
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/Jastrow/J1OrbitalSoA.h"
#include "QMCWaveFunctions/Jastrow/RadialJastrowBuilder.h"
#include "QMCWaveFunctions/WaveFunctionFactory.h"

namespace qmcplusplus
{
TEST_CASE("J1 evaluate derivatives Jastrow", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  ParticleSetPool ptcl = ParticleSetPool(c);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
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

  ions_.get(app_log());
  elec_.get(app_log());

  const char* jasxml = "<wavefunction name=\"psi0\" target=\"e\"> \
  <jastrow name=\"J1\" type=\"One-Body\" function=\"Bspline\" print=\"yes\" source=\"ion0\"> \
    <correlation elementType=\"H\" cusp=\"0.0\" size=\"2\" rcut=\"5.0\"> \
      <coefficients id=\"J1H\" type=\"Array\"> 0.5 0.1 </coefficients> \
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
  REQUIRE(nparam == 2);

  using ValueType = QMCTraits::ValueType;
  std::vector<ValueType> dlogpsi(nparam);
  std::vector<ValueType> dhpsioverpsi(nparam);
  //twf.evaluateDerivatives(elec_, active, dlogpsi, dhpsioverpsi);
  twf_component_list[0]->evaluateDerivatives(elec_, active, dlogpsi, dhpsioverpsi);

  // Numbers not validated
  std::vector<ValueType> expected_dlogpsi      = {-0.9336294487, -1.0196051794};
  std::vector<ValueType> expected_dhpsioverpsi = {-1.1596433096, 0.7595492539};
  for (int i = 0; i < nparam; i++)
  {
    CHECK(dlogpsi[i] == ValueApprox(expected_dlogpsi[i]));
    CHECK(dhpsioverpsi[i] == ValueApprox(expected_dhpsioverpsi[i]));
  }
}

TEST_CASE("J1 evaluate derivatives Jastrow with two species", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;
  ParticleSetPool ptcl = ParticleSetPool(c);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion0");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({1, 1});
  ions_.R[0]                 = {0.0, 0.0, 1.0};
  ions_.R[1]                 = {0.0, 0.0, 0.0};
  SpeciesSet& ispecies       = ions_.getSpeciesSet();
  int OIdx                   = ispecies.addSpecies("O");
  int HIdx                   = ispecies.addSpecies("H");
  int imassIdx               = ispecies.addAttribute("mass");
  int ichargeIdx             = ispecies.addAttribute("charge");
  ispecies(ichargeIdx, HIdx) = 1.0;
  ispecies(ichargeIdx, OIdx) = 8.0;
  ispecies(imassIdx, HIdx)   = 1.0;
  ispecies(imassIdx, OIdx)   = 16.0;

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

  ions_.get(app_log());
  elec_.get(app_log());

  const char* jasxml = "<wavefunction name=\"psi0\" target=\"e\"> \
  <jastrow name=\"J1\" type=\"One-Body\" function=\"Bspline\" print=\"yes\" source=\"ion0\"> \
    <correlation elementType=\"H\" cusp=\"0.0\" size=\"2\" rcut=\"5.0\"> \
      <coefficients id=\"J1H\" type=\"Array\"> 0.5 0.1 </coefficients> \
    </correlation> \
    <correlation elementType=\"O\" cusp=\"0.0\" size=\"2\" rcut=\"5.0\"> \
      <coefficients id=\"J1O\" type=\"Array\"> 0.2 0.1 </coefficients> \
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
  REQUIRE(nparam == 4);

  using ValueType = QMCTraits::ValueType;
  std::vector<ValueType> dlogpsi(nparam);
  std::vector<ValueType> dhpsioverpsi(nparam);
  //twf.evaluateDerivatives(elec_, active, dlogpsi, dhpsioverpsi);
  twf_component_list[0]->evaluateDerivatives(elec_, active, dlogpsi, dhpsioverpsi);

  // Numbers not validated independently
  std::vector<ValueType> expected_dlogpsi      = {-0.6360001724, -1.1764442146, -0.9336294487, -1.0196051794};
  std::vector<ValueType> expected_dhpsioverpsi = {-0.6225838942, 0.0099980417, -1.1853702074, 0.7798000176};
  for (int i = 0; i < nparam; i++)
  {
    CHECK(dlogpsi[i] == ValueApprox(expected_dlogpsi[i]));
    CHECK(dhpsioverpsi[i] == ValueApprox(expected_dhpsioverpsi[i]));
  }
}

TEST_CASE("J1 evaluate derivatives Jastrow with two species one without Jastrow", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;
  ParticleSetPool ptcl = ParticleSetPool(c);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion0");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({1, 1});
  ions_.R[0]                 = {0.0, 0.0, 1.0};
  ions_.R[1]                 = {0.0, 0.0, 0.0};
  SpeciesSet& ispecies       = ions_.getSpeciesSet();
  int OIdx                   = ispecies.addSpecies("O");
  int HIdx                   = ispecies.addSpecies("H");
  int imassIdx               = ispecies.addAttribute("mass");
  int ichargeIdx             = ispecies.addAttribute("charge");
  ispecies(ichargeIdx, HIdx) = 1.0;
  ispecies(ichargeIdx, OIdx) = 8.0;
  ispecies(imassIdx, HIdx)   = 1.0;
  ispecies(imassIdx, OIdx)   = 16.0;

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

  ions_.get(app_log());
  elec_.get(app_log());

  const char* jasxml = "<wavefunction name=\"psi0\" target=\"e\"> \
  <jastrow name=\"J1\" type=\"One-Body\" function=\"Bspline\" print=\"yes\" source=\"ion0\"> \
    <correlation elementType=\"H\" cusp=\"0.0\" size=\"2\" rcut=\"5.0\"> \
      <coefficients id=\"J1H\" type=\"Array\"> 0.5 0.1 </coefficients> \
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
  REQUIRE(nparam == 2);

  using ValueType = QMCTraits::ValueType;
  std::vector<ValueType> dlogpsi(nparam);
  std::vector<ValueType> dhpsioverpsi(nparam);
  //twf.evaluateDerivatives(elec_, active, dlogpsi, dhpsioverpsi);
  twf_component_list[0]->evaluateDerivatives(elec_, active, dlogpsi, dhpsioverpsi);

  // Numbers not validated independently
  std::vector<ValueType> expected_dlogpsi      = {-0.9336294487, -1.0196051794};
  std::vector<ValueType> expected_dhpsioverpsi = {-1.1596433096, 0.7595492539};
  for (int i = 0; i < nparam; i++)
  {
    CHECK(dlogpsi[i] == ValueApprox(expected_dlogpsi[i]));
    CHECK(dhpsioverpsi[i] == ValueApprox(expected_dhpsioverpsi[i]));
  }
}
} // namespace qmcplusplus
