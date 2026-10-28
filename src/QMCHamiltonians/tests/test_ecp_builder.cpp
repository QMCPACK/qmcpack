//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
// Copyright (c) 2026 QMCPACK developers.
//////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>

#include "Message/UniformCommunicateError.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include "QMCHamiltonians/ECPotentialBuilder.h"
#include "QMCHamiltonians/HamiltonianFactory.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Utilities/RuntimeOptions.h"

namespace qmcplusplus
{
TEST_CASE("ECPotentialBuilder validates ion species metadata", "[hamiltonian][ecp]")
{
  bool add_charge              = true;
  double charge                = 4;
  std::string atomic_attribute = "atomicnumber";
  int atomic_number            = 6;
  std::string element          = "C";
  std::string format           = "xml";
  bool valid                   = false;
  std::string diagnostic;

  SECTION("matching metadata") { valid = true; }
  SECTION("atomic_number alias")
  {
    atomic_attribute = "atomic_number";
    valid            = true;
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
  SECTION("missing charge")
  {
    add_charge = false;
    diagnostic = "Ion species C needs parameter 'charge'";
  }
  SECTION("charge above tolerance")
  {
    charge += 2e-4;
    diagnostic = "Ion species C charge 4.0002 pseudopotential charge 4 mismatch!";
  }
  SECTION("charge below tolerance")
  {
    charge -= 2e-4;
    diagnostic = "Ion species C charge 3.9998 pseudopotential charge 4 mismatch!";
  }
  SECTION("missing atomic number")
  {
    atomic_attribute.clear();
    diagnostic = "Ion species C needs parameter 'atomicnumber'";
  }
  SECTION("atomic number mismatch")
  {
    atomic_number = 5;
    diagnostic    = "Ion species C atomicnumber 5 pseudopotential atomic-number 6 mismatch!";
  }
  SECTION("atomic_number alias mismatch")
  {
    atomic_attribute = "atomic_number";
    atomic_number    = 5;
    diagnostic       = "Ion species C atomicnumber 5 pseudopotential atomic-number 6 mismatch!";
  }
  SECTION("unknown ion species")
  {
    element    = "N";
    diagnostic = "Ion species N is not found.";
  }
  SECTION("unsupported pseudo format") { format = "unsupported"; }

  SECTION("missing charge and atomic number")
  {
    add_charge = false;
    atomic_attribute.clear();
    diagnostic = "Ion species C needs parameter 'charge'";
  }

  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell);
  ParticleSet electrons(simulation_cell);
  ions.setName("ion0");
  ions.create({1});
  electrons.setName("e");
  electrons.create({1});
  auto& species = ions.getSpeciesSet();
  int carbon    = species.addSpecies("C");
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
    std::string expected_message = "ECPotentialBuilder::useXmlFormat failed:\n";
    if (!diagnostic.empty())
      expected_message += "  " + diagnostic + "\n";
    if (!add_charge && atomic_attribute.empty())
      expected_message += "  Ion species C needs parameter 'atomicnumber'\n";
    expected_message += "  Failed to add pseudopotential for element " + element + "\n";

    CHECK_THROWS_MATCHES(builder.put(doc.getRoot()), UniformCommunicateError,
                         Catch::Matchers::Message(expected_message));
    CHECK(hamiltonian.size() == 0);
  }
}

namespace
{
void checkInvalidFactoryInput()
{
  ParticleSetPool particle_pool(*OHMMS::Controller);
  const auto& simulation_cell = particle_pool.getSimulationCell();
  auto ions                   = std::make_unique<ParticleSet>(simulation_cell);
  auto electrons              = std::make_unique<ParticleSet>(simulation_cell);
  ions->setName("ion0");
  ions->create({1});
  auto& species                                         = ions->getSpeciesSet();
  int carbon                                            = species.addSpecies("C");
  species(species.addAttribute("charge"), carbon)       = 3;
  species(species.addAttribute("atomicnumber"), carbon) = 6;
  electrons->setName("e");
  electrons->create({1});
  auto& electron_species                                      = electrons->getSpeciesSet();
  int up                                                      = electron_species.addSpecies("u");
  electron_species(electron_species.addAttribute("mass"), up) = 1;
  auto& target                                                = *electrons;
  particle_pool.addParticleSet(std::move(ions));
  particle_pool.addParticleSet(std::move(electrons));
  RuntimeOptions options;
  Libxml2Document doc;
  REQUIRE(doc.parseFromString(R"(<hamiltonian name="h0" type="generic" target="e">
    <pairpot type="pseudo" name="PseudoPot" source="ion0" format="xml">
      <pseudo elementType="C" href="C.BFD.xml"/>
    </pairpot>
  </hamiltonian>)"));
  TrialWaveFunction psi(options, "psi0");
  HamiltonianFactory factory("h0", target, particle_pool.getPool(), psi, *OHMMS::Controller);
  REQUIRE_THROWS_MATCHES(factory.put(doc.getRoot()), UniformCommunicateError,
                         Catch::Matchers::Message("ECPotentialBuilder::useXmlFormat failed:\n"
                                                  "  Ion species C charge 3 pseudopotential charge 4 mismatch!\n"
                                                  "  Failed to add pseudopotential for element C\n"));
}
} // namespace

TEST_CASE("HamiltonianFactory throws on invalid ECP metadata", "[hamiltonian][ecp]") { checkInvalidFactoryInput(); }
} // namespace qmcplusplus
