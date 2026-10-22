//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
// Copyright (c) 2026 QMCPACK developers.
//////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_test_macros.hpp>

#include "Message/UniformCommunicateError.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Particle/ParticleSet.h"
#include "QMCHamiltonians/ECPotentialBuilder.h"
#include "QMCHamiltonians/HamiltonianFactory.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Utilities/RuntimeOptions.h"

namespace qmcplusplus
{
TEST_CASE("ECPotentialBuilder validates ion species metadata", "[hamiltonian][ecp]")
{
  bool add_charge = true;
  double charge = 4;
  std::string atomic_attribute = "atomicnumber";
  int atomic_number = 6;
  std::string element = "C";
  std::string format = "xml";
  bool valid = false;

  SECTION("matching metadata") { valid = true; }
  SECTION("atomic_number alias")
  {
    atomic_attribute = "atomic_number";
    valid = true;
  }
  SECTION("charge within positive tolerance")
  {
    charge += 5e-5;
    valid = true;
  }
  SECTION("charge within negative tolerance")
  {
    charge -= 5e-5;
    valid = true;
  }
  SECTION("missing charge") { add_charge = false; }
  SECTION("charge above tolerance") { charge += 2e-4; }
  SECTION("charge below tolerance") { charge -= 2e-4; }
  SECTION("missing atomic number") { atomic_attribute.clear(); }
  SECTION("atomic number mismatch") { atomic_number = 5; }
  SECTION("atomic_number alias mismatch")
  {
    atomic_attribute = "atomic_number";
    atomic_number = 5;
  }
  SECTION("unknown ion species") { element = "N"; }
  SECTION("unsupported pseudo format") { format = "unsupported"; }

  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell);
  ParticleSet electrons(simulation_cell);
  ions.setName("ion0");
  ions.create({1});
  electrons.setName("e");
  electrons.create({1});
  auto& species = ions.getSpeciesSet();
  int carbon = species.addSpecies("C");
  if (add_charge)
    species(species.addAttribute("charge"), carbon) = charge;
  if (!atomic_attribute.empty())
    species(species.addAttribute(atomic_attribute), carbon) = atomic_number;

  QMCHamiltonian hamiltonian;
  ECPotentialBuilder builder(hamiltonian, ions, electrons, *OHMMS::Controller);
  Libxml2Document doc;
  REQUIRE(doc.parseFromString("<pairpot format=\"xml\"><pseudo elementType=\"" + element +
                             "\" href=\"C.BFD.xml\" format=\"" + format + "\"/></pairpot>"));
  if (valid)
  {
    REQUIRE(builder.put(doc.getRoot()));
    REQUIRE(hamiltonian.size() == 2);
    CHECK(hamiltonian.getComponent(0)->getName() == "LocalECP");
    CHECK(hamiltonian.getComponent(1)->getName() == "NonLocalECP");
  }
  else
  {
    try
    {
      builder.put(doc.getRoot());
      FAIL("Invalid ECP input must throw UniformCommunicateError");
    }
    catch (const UniformCommunicateError& error)
    {
      CHECK(std::string(error.what()) ==
            "ECPotentialBuilder::useXmlFormat failed: Failed to add pseudopotential for element " + element);
    }
    CHECK(hamiltonian.size() == 0);
  }
}

// Run only in a separate process: the factory must convert the builder exception into an abort.
TEST_CASE("HamiltonianFactory aborts on invalid ECP metadata", "[.ecp_abort]")
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
