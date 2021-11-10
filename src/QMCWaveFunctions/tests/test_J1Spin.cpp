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

namespace qmcplusplus
{
using RealType     = WaveFunctionComponent::RealType;
using LogValueType = WaveFunctionComponent::LogValueType;
using ValueType    = QMCTraits::ValueType;
#if !defined(QMC_CUDA)
TEST_CASE("J1 spin evaluate derivatives Jastrow", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  auto ions_uptr = std::make_unique<ParticleSet>();
  auto elec_uptr = std::make_unique<ParticleSet>();
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion0");
  ions_.create(1);
  ions_.R[0]                 = {0.0, 0.0, 0.0};
  SpeciesSet& ispecies       = ions_.getSpeciesSet();
  int HIdx                   = ispecies.addSpecies("H");
  int ichargeIdx             = ispecies.addAttribute("charge");
  ispecies(ichargeIdx, HIdx) = 1.0;

  elec_.setName("e");
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

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.addParticleSet(std::move(elec_uptr));
  ptcl.addParticleSet(std::move(ions_uptr));

  ions_.update();
  elec_.addTable(elec_);
  elec_.addTable(ions_);
  elec_.update();

  const char* jasxml = "<wavefunction name=\"psi0\" target=\"e\"> \
<jastrow name=\"J1\" type=\"One-Body\" function=\"Bspline\" print=\"yes\" source=\"ion0\" spin=\"yes\"> \
  <correlation speciesA=\"H\" speciesB=\"u\" cusp=\"0.0\" size=\"2\" rcut=\"5.0\"> \
    <coefficients id=\"J1uH\" type=\"Array\"> 0.5 0.1 </coefficients> \
  </correlation> \
  <correlation speciesA=\"H\" speciesB=\"d\" cusp=\"0.0\" size=\"2\" rcut=\"5.0\"> \
    <coefficients id=\"J1dH\" type=\"Array\"> 0.5 0.1 </coefficients> \
  </correlation> \
</jastrow> \
</wavefunction> \
";
  Libxml2Document doc;
  bool okay = doc.parseFromString(jasxml);
  REQUIRE(okay);
  xmlNodePtr jas1 = doc.getRoot();
  WaveFunctionFactory wf_factory("psi0", elec_, ptcl.getPool(), c);
  wf_factory.put(jas1);
  auto& twf(*wf_factory.getTWF());
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
  std::vector<ValueType> dlogpsi(nparam);
  std::vector<ValueType> dhpsioverpsi(nparam);
  std::vector<ValueType> cloned_dlogpsi(nparam);
  std::vector<ValueType> cloned_dhpsioverpsi(nparam);

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
#endif
} // namespace qmcplusplus
