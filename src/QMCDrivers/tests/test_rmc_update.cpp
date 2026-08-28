//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>

#include "OhmmsData/Libxml2Doc.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCDrivers/RMC/RMCUpdateAll.h"
#include "QMCDrivers/RMC/RMCUpdatePbyP.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/RuntimeOptions.h"

namespace qmcplusplus
{
namespace testing
{
class RMCUpdateAllTests
{
public:
  static void check(const RMCUpdateAllWithDrift& update, bool scale_drift, int action_type, int equil_steps)
  {
    CHECK(update.scaleDrift == scale_drift);
    CHECK(update.actionType == action_type);
    CHECK(update.equilSteps == equil_steps);
  }
};
} // namespace testing

TEST_CASE("RMC particle-by-particle retired input controls", "[drivers][rmc]")
{
  const SimulationCell simulation_cell;
  MCWalkerConfiguration elec(simulation_cell);
  elec.setName("elec");
  elec.create({1});

  SpeciesSet& species         = elec.getSpeciesSet();
  int up_idx                  = species.addSpecies("u");
  int charge_idx              = species.addAttribute("charge");
  int mass_idx                = species.addAttribute("mass");
  species(charge_idx, up_idx) = -1;
  species(mass_idx, up_idx)   = 1.0;

  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);
  QMCHamiltonian h;
  FakeRandom<QMCTraits::FullPrecRealType> rng;
  RMCUpdatePbyPWithDrift update(elec, psi, h, rng, {}, {});

  Libxml2Document doc;
  REQUIRE(doc.parseFromString(R"(<qmc>
    <parameter name="useDrift">not-a-boolean</parameter>
    <parameter name="Action">DMC</parameter>
    <parameter name="equilsteps">not-an-integer</parameter>
    <parameter name="equilSteps">also-not-an-integer</parameter>
    <parameter name="debug_checks">checkGL_after_moves</parameter>
  </qmc>)"));

  // These names are intentionally ignored by the PbyP mover. Common live
  // QMCUpdateBase controls continue to be parsed through the inherited put().
  REQUIRE(update.put(doc.getRoot()));
  CHECK(update.debug_checks_ & DriverDebugChecks::CHECKGL_AFTER_MOVES);
}

TEST_CASE("RMC all-particle input controls remain live", "[drivers][rmc]")
{
  const SimulationCell simulation_cell;
  MCWalkerConfiguration elec(simulation_cell);
  elec.setName("elec");
  elec.create({1});

  SpeciesSet& species         = elec.getSpeciesSet();
  int up_idx                  = species.addSpecies("u");
  int charge_idx              = species.addAttribute("charge");
  int mass_idx                = species.addAttribute("mass");
  species(charge_idx, up_idx) = -1;
  species(mass_idx, up_idx)   = 1.0;

  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);
  QMCHamiltonian h;
  FakeRandom<QMCTraits::FullPrecRealType> rng;

  SECTION("lowercase equilsteps alias with scaled DMC action")
  {
    RMCUpdateAllWithDrift update(elec, psi, h, rng, {}, {});
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(R"(<qmc>
      <parameter name="useScaledDrift">yes</parameter>
      <parameter name="Action">DMC</parameter>
      <parameter name="equilsteps">7</parameter>
    </qmc>)"));
    REQUIRE(update.put(doc.getRoot()));
    testing::RMCUpdateAllTests::check(update, true, RMCUpdateAllWithDrift::DMC_ACTION, 7);
  }

  SECTION("camel-case equilSteps alias with symmetric action")
  {
    RMCUpdateAllWithDrift update(elec, psi, h, rng, {}, {});
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(R"(<qmc>
      <parameter name="useScaledDrift">no</parameter>
      <parameter name="Action">SLA</parameter>
      <parameter name="equilSteps">9</parameter>
    </qmc>)"));
    REQUIRE(update.put(doc.getRoot()));
    testing::RMCUpdateAllTests::check(update, false, RMCUpdateAllWithDrift::SYM_ACTION, 9);
  }
}
} // namespace qmcplusplus
