//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
// Copyright (c) 2026 QMCPACK developers.
//////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_test_macros.hpp>

#include "OhmmsData/Libxml2Doc.h"
#include "Particle/ParticleSet.h"
#include "QMCHamiltonians/ECPotentialBuilder.h"
#include "QMCHamiltonians/HamiltonianFactory.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Utilities/RuntimeOptions.h"

namespace qmcplusplus
{
namespace
{
struct ECPInput
{
  bool add_charge = true;
  double charge = 4;
  std::string atomic_attribute = "atomicnumber";
  int atomic_number = 6;
  std::string element = "C";
  std::string format = "xml";
};

void checkECPInput(const ECPInput& input, bool expect_abort)
{
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell);
  ParticleSet electrons(simulation_cell);
  ions.setName("ion0");
  ions.create({1});
  electrons.setName("e");
  electrons.create({1});
  auto& species = ions.getSpeciesSet();
  int carbon = species.addSpecies("C");
  if (input.add_charge)
    species(species.addAttribute("charge"), carbon) = input.charge;
  if (!input.atomic_attribute.empty())
    species(species.addAttribute(input.atomic_attribute), carbon) = input.atomic_number;

  QMCHamiltonian hamiltonian;
  ECPotentialBuilder builder(hamiltonian, ions, electrons, *OHMMS::Controller);
  Libxml2Document doc;
  REQUIRE(doc.parseFromString("<pairpot format=\"xml\"><pseudo elementType=\"" + input.element +
                             "\" href=\"C.BFD.xml\" format=\"" + input.format + "\"/></pairpot>"));
  REQUIRE(builder.put(doc.getRoot()));
  if (expect_abort)
    FAIL("Invalid ECP input must abort");
  REQUIRE(hamiltonian.size() == 2);
  CHECK(hamiltonian.getComponent(0)->getName() == "LocalECP");
  CHECK(hamiltonian.getComponent(1)->getName() == "NonLocalECP");
}
} // namespace

TEST_CASE("ECPotentialBuilder accepts matching ion species metadata", "[hamiltonian][ecp]")
{
  ECPInput input;
  SECTION("matching metadata") {}
  SECTION("atomic_number alias") { input.atomic_attribute = "atomic_number"; }
  SECTION("charge within positive tolerance") { input.charge += 5e-5; }
  SECTION("charge within negative tolerance") { input.charge -= 5e-5; }
  checkECPInput(input, false);
}

// Each section is selected explicitly by a separate CTest subprocess.
TEST_CASE("ECPotentialBuilder aborts on invalid ion species metadata", "[.ecp_abort]")
{
  ECPInput input;
  SECTION("missing_charge") { input.add_charge = false; }
  SECTION("charge_above_tolerance") { input.charge += 2e-4; }
  SECTION("charge_below_tolerance") { input.charge -= 2e-4; }
  SECTION("missing_atomic_number") { input.atomic_attribute.clear(); }
  SECTION("atomic_number_mismatch") { input.atomic_number = 5; }
  SECTION("atomic_number_alias_mismatch")
  {
    input.atomic_attribute = "atomic_number";
    input.atomic_number = 5;
  }
  SECTION("unknown_species") { input.element = "N"; }
  SECTION("unsupported_format") { input.format = "unsupported"; }
  checkECPInput(input, true);
}

// Run only in a separate process: invalid metadata must abort through the factory.
TEST_CASE("HamiltonianFactory aborts on invalid ECP metadata", "[.ecp_factory_abort]")
{
  const SimulationCell simulation_cell;
  auto ions = std::make_unique<ParticleSet>(simulation_cell);
  auto electrons = std::make_unique<ParticleSet>(simulation_cell);
  ions->setName("ion0");
  ions->create({1});
  auto& species = ions->getSpeciesSet();
  int carbon = species.addSpecies("C");
  species(species.addAttribute("charge"), carbon) = 3;
  species(species.addAttribute("atomicnumber"), carbon) = 6;
  electrons->setName("e");
  electrons->create({1});
  auto& electron_species = electrons->getSpeciesSet();
  int up = electron_species.addSpecies("u");
  electron_species(electron_species.addAttribute("mass"), up) = 1;
  auto& target = *electrons;
  HamiltonianFactory::PSetMap particles;
  particles.emplace("ion0", std::move(ions));
  particles.emplace("e", std::move(electrons));
  RuntimeOptions options;
  TrialWaveFunction psi(options, "psi0");
  HamiltonianFactory factory("h0", target, particles, psi, *OHMMS::Controller);
  Libxml2Document doc;
  REQUIRE(doc.parseFromString(R"(<hamiltonian name="h0" type="generic" target="e">
    <pairpot type="pseudo" name="PseudoPot" source="ion0" format="xml">
      <pseudo elementType="C" href="C.BFD.xml"/>
    </pairpot>
  </hamiltonian>)"));
  factory.put(doc.getRoot());
  FAIL("HamiltonianFactory must abort on invalid ECP metadata");
}
} // namespace qmcplusplus
