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
#include "QMCWaveFunctions/Jastrow/J1Spin.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "QMCWaveFunctions/Jastrow/RadialJastrowBuilder.h"
#include "QMCWaveFunctions/WaveFunctionFactory.h"
#include "Utilities/RuntimeOptions.h"

namespace qmcplusplus
{
using RealType     = WaveFunctionComponent::RealType;
using LogValueType = WaveFunctionComponent::LogValueType;
using ValueType    = QMCTraits::ValueType;

TEST_CASE("J1 spin evaluate derivatives Jastrow", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  ParticleSetPool ptcl = ParticleSetPool(c);
  auto ions_uptr       = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr       = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
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

  const char* jasxml = R"(<wavefunction name="psi0" target="e">
<jastrow name="J1" type="One-Body" function="Bspline" print="yes" source="ion0" spin="yes">
  <correlation speciesA="H" speciesB="u" cusp="0.0" size="2" rcut="5.0">
    <coefficients id="J1uH" type="Array"> 0.5 0.1 </coefficients>
  </correlation>
  <correlation speciesA="H" speciesB="d" cusp="0.0" size="2" rcut="5.0">
    <coefficients id="J1dH" type="Array"> 0.5 0.1 </coefficients>
  </correlation>
</jastrow>
</wavefunction>
)";
  Libxml2Document doc;
  bool okay = doc.parseFromString(jasxml);
  REQUIRE(okay);
  xmlNodePtr jas1 = doc.getRoot();
  WaveFunctionFactory wf_factory(elec_, ptcl.getPool(), c);
  RuntimeOptions runtime_options;
  auto twf_ptr = wf_factory.buildTWF(jas1, runtime_options);
  auto& twf(*twf_ptr);
  twf.setMassTerm(elec_);
  auto& twf_component_list = twf.getOrbitals();
  auto cloned_j1spin       = twf_component_list[0]->makeClone(elec_);

  opt_variables_type active;
  twf.checkInVariables(active);
  active.removeInactive();
  int nparam = active.size_of_active();
  REQUIRE(nparam == 4);

  // check logs
  //evaluateLog += into G + L so reset
  elec_.G          = 0.0;
  elec_.L          = 0.0;
  LogValueType log = twf_component_list[0]->evaluateLog(elec_, elec_.G, elec_.L);
  LogValueType expected_log{-0.568775, 0.0};
  CHECK(log == LogComplexApprox(expected_log));
  //evaluateLog += into G + L so reset
  elec_.G                 = 0.0;
  elec_.L                 = 0.0;
  LogValueType cloned_log = cloned_j1spin->evaluateLog(elec_, elec_.G, elec_.L);
  CHECK(cloned_log == LogComplexApprox(expected_log));

  // check derivatives
  twf.evaluateLog(elec_);
  Vector<ValueType> dlogpsi(nparam);
  Vector<ValueType> dhpsioverpsi(nparam);
  Vector<ValueType> cloned_dlogpsi(nparam);
  Vector<ValueType> cloned_dhpsioverpsi(nparam);

  twf_component_list[0]->evaluateDerivatives(elec_, active, dlogpsi, dhpsioverpsi);
  cloned_j1spin->evaluateDerivatives(elec_, active, cloned_dlogpsi, cloned_dhpsioverpsi);
  // Numbers not validated
  std::vector<ValueType> expected_dlogpsi      = {-0.46681472435, -0.5098025897, -0.46681472435, -0.5098025897};
  std::vector<ValueType> expected_dhpsioverpsi = {-0.5798216548, 0.37977462695, -0.5798216548, 0.37977462695};
  for (int i = 0; i < nparam; i++)
  {
    CHECK(dlogpsi[i] == ValueApprox(expected_dlogpsi[i]));
    CHECK(cloned_dlogpsi[i] == ValueApprox(expected_dlogpsi[i]));
    CHECK(dhpsioverpsi[i] == ValueApprox(expected_dhpsioverpsi[i]));
    CHECK(cloned_dhpsioverpsi[i] == ValueApprox(expected_dhpsioverpsi[i]));
  }
}

TEST_CASE("J1 spin evaluate derivatives multiparticle Jastrow", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  ParticleSetPool ptcl = ParticleSetPool(c);
  auto ions_uptr       = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr       = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion0");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({2});
  ions_.R[0]                 = {-1.0, 0.0, 0.0};
  ions_.R[1]                 = { 1.0, 0.0, 0.0};
  SpeciesSet& ispecies       = ions_.getSpeciesSet();
  int BeIdx                   = ispecies.addSpecies("Be");
  int ichargeIdx             = ispecies.addAttribute("charge");
  ispecies(ichargeIdx, BeIdx) = 4.0;

  elec_.setName("e");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({4, 4, 1});
  elec_.R[0] = { 0.5,  0.5,  0.5};
  elec_.R[1] = {-0.5,  0.5,  0.5};
  elec_.R[2] = { 0.5, -0.5,  0.5};
  elec_.R[3] = { 0.5,  0.5, -0.5};
  elec_.R[4] = {-0.5, -0.5,  0.5};
  elec_.R[5] = { 0.5, -0.5, -0.5};
  elec_.R[6] = {-0.5,  0.5, -0.5};
  elec_.R[7] = {-0.5, -0.5, -0.5};
  elec_.R[8] = { 1.5,  1.5,  1.5};

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int downIdx                = tspecies.addSpecies("d");
  int posIdx                = tspecies.addSpecies("p");
  int massIdx                = tspecies.addAttribute("mass");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(massIdx, upIdx)   = 1.0;
  tspecies(massIdx, downIdx) = 1.0;
  tspecies(massIdx, posIdx) = 1.0;
  tspecies(chargeIdx, upIdx) = -1.0;
  tspecies(massIdx, downIdx) = -1.0;
  tspecies(massIdx, posIdx) =  1.0;
  // Necessary to set mass
  elec_.resetGroups();

  ions_.update();
  elec_.addTable(elec_);
  elec_.addTable(ions_);
  elec_.update();

  const char* jasxml = R"(<wavefunction name="psi0" target="e">
<jastrow name="J1" type="One-Body" function="Bspline" print="yes" source="ion0" spin="yes">
  <correlation speciesA="Be" speciesB="u" cusp="0.0" size="2" rcut="5.0">
    <coefficients id="J1uH" type="Array"> 0.5 0.1 </coefficients>
  </correlation>
  <correlation speciesA="Be" speciesB="d" cusp="0.0" size="2" rcut="5.0">
    <coefficients id="J1dH" type="Array"> 0.5 0.1 </coefficients>
  </correlation>
  <correlation speciesA="Be" speciesB="p" cusp="0.0" size="2" rcut="5.0">
    <coefficients id="J1pH" type="Array"> 0.5 0.1 </coefficients>
  </correlation>
</jastrow>
</wavefunction>
)";
  Libxml2Document doc;
  bool okay = doc.parseFromString(jasxml);
  REQUIRE(okay);
  xmlNodePtr jas1 = doc.getRoot();
  WaveFunctionFactory wf_factory(elec_, ptcl.getPool(), c);
  RuntimeOptions runtime_options;
  auto twf_ptr = wf_factory.buildTWF(jas1, runtime_options);
  auto& twf(*twf_ptr);
  twf.setMassTerm(elec_);
  auto& twf_component_list = twf.getOrbitals();
  auto cloned_j1spin       = twf_component_list[0]->makeClone(elec_);

  opt_variables_type active;
  twf.checkInVariables(active);
  active.removeInactive();
  int nparam = active.size_of_active();
  REQUIRE(nparam == 6);

  // check logs
  //evaluateLog += into G + L so reset
  elec_.G          = 0.0;
  elec_.L          = 0.0;
  LogValueType log = twf_component_list[0]->evaluateLog(elec_, elec_.G, elec_.L);
  LogValueType expected_log{-3.58983, 0.0};
  CHECK(log == LogComplexApprox(expected_log));
  //evaluateLog += into G + L so reset
  elec_.G                 = 0.0;
  elec_.L                 = 0.0;
  LogValueType cloned_log = cloned_j1spin->evaluateLog(elec_, elec_.G, elec_.L);
  CHECK(cloned_log == LogComplexApprox(expected_log));

  // check derivatives
  twf.evaluateLog(elec_);
  Vector<ValueType> dlogpsi(nparam);
  Vector<ValueType> dhpsioverpsi(nparam);
  Vector<ValueType> cloned_dlogpsi(nparam);
  Vector<ValueType> cloned_dhpsioverpsi(nparam);

  twf_component_list[0]->evaluateDerivatives(elec_, active, dlogpsi, dhpsioverpsi);
  cloned_j1spin->evaluateDerivatives(elec_, active, cloned_dlogpsi, cloned_dhpsioverpsi);
  // Numbers not validated
  std::vector<ValueType> expected_dlogpsi      = {-2.544, -4.70578, -2.544, -4.70578, -0.055314, -0.770138};
  std::vector<ValueType> expected_dhpsioverpsi = {-2.45001, 0.0794429, -2.45001, 0.0794429, 0.0462761, -0.330801};
  for (int i = 0; i < nparam; i++)
  {
    CHECK(dlogpsi[i] == ValueApprox(expected_dlogpsi[i]));
    CHECK(cloned_dlogpsi[i] == ValueApprox(expected_dlogpsi[i]));
    CHECK(dhpsioverpsi[i] == ValueApprox(expected_dhpsioverpsi[i]));
    CHECK(cloned_dhpsioverpsi[i] == ValueApprox(expected_dhpsioverpsi[i]));
  }
}
} // namespace qmcplusplus
