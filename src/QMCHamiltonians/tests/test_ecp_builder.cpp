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
#include "QMCHamiltonians/ECPComponentBuilder.h"
#include "QMCHamiltonians/HamiltonianFactory.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Utilities/RuntimeOptions.h"
#include <filesystem>
#include <fstream>

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

TEST_CASE("ECPotentialBuilder unsupported format throws", "[hamiltonian][ecp]")
{
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell);
  ParticleSet electrons(simulation_cell);
  QMCHamiltonian hamiltonian;
  ECPotentialBuilder builder(hamiltonian, ions, electrons, *OHMMS::Controller);
  Libxml2Document doc;
  REQUIRE(doc.parseFromString(R"(<pairpot format="table"><pseudo elementType="C" href="C.BFD.xml"/></pairpot>)"));
  CHECK_THROWS_MATCHES(builder.put(doc.getRoot()), UniformCommunicateError,
                       Catch::Matchers::Message("Only XML format ECPs are currently supported."));
}

TEST_CASE("ECPComponentBuilder file reading exceptions", "[hamiltonian][ecp]")
{
  Communicate& c(*OHMMS::Controller);
  ECPComponentBuilder ecp("test_ecp", c);

  SECTION("missing file")
  {
    CHECK_THROWS_MATCHES(ecp.read_pp_file("nonexistent_pp_file.xml"), UniformCommunicateError,
                         Catch::Matchers::Message("ECPComponentBuilder::read_pp_file  Missing PP file nonexistent_pp_file.xml\n"));
  }

  SECTION("invalid xml")
  {
    const std::string bad_file = "invalid_test_pp.xml";
    if (c.rank() == 0)
    {
      std::ofstream fout(bad_file);
      fout << "<pseudo><unclosed>";
    }
    c.barrier();
    CHECK_THROWS_MATCHES(ecp.read_pp_file(bad_file), UniformCommunicateError,
                         Catch::Matchers::Message("ECPComponentBuilder::read_pp_file xml file " + bad_file + " is invalid"));
    c.barrier();
    if (c.rank() == 0)
      std::filesystem::remove(bad_file);
  }

  SECTION("empty document")
  {
    const std::string empty_file = "empty_test_pp.xml";
    if (c.rank() == 0)
    {
      std::ofstream fout(empty_file);
      fout << "<?xml version=\"1.0\"?>\n";
    }
    c.barrier();
    // Some XML parsers consider a document without root elements invalid; handle either expected UniformCommunicateError
    try
    {
      ecp.read_pp_file(empty_file);
      FAIL("Expected UniformCommunicateError");
    }
    catch (const UniformCommunicateError& uce)
    {
      std::string msg(uce.what());
      CHECK((msg == "Empty document" ||
             msg == "ECPComponentBuilder::read_pp_file xml file " + empty_file + " is invalid"));
    }
    c.barrier();
    if (c.rank() == 0)
      std::filesystem::remove(empty_file);
  }

  SECTION("parseCasino missing file")
  {
    Libxml2Document doc;
    REQUIRE(doc.parseFromString("<pseudo/>"));
    CHECK_THROWS_MATCHES(ecp.parseCasino("nonexistent_casino_file.dat", doc.getRoot()), UniformCommunicateError,
                         Catch::Matchers::Message("ECPComponentBuilder::parseCasino: Could not open file nonexistent_casino_file.dat"));
  }
}

TEST_CASE("ECPComponentBuilder L2 exceptions", "[hamiltonian][ecp]")
{
  Communicate& c(*OHMMS::Controller);
  ECPComponentBuilder ecp("test_ecp_l2", c);
  Libxml2Document doc;

  SECTION("no global grid")
  {
    REQUIRE(doc.parseFromString(R"(<L2 format="V" cutoff="1.0"/>)"));
    CHECK_THROWS_MATCHES(ecp.buildL2(doc.getRoot()), UniformCommunicateError,
                         Catch::Matchers::Message("ECPComponentBuilder::buildL2: Global grid needs to be defined."));
  }

  SECTION("unrecognized format")
  {
    ecp.grid_global = std::make_unique<LinearGrid<ECPComponentBuilder::mRealType>>();
    ecp.grid_global->set(0.0, 1.0, 10);
    REQUIRE(doc.parseFromString(R"(<L2 format="unknown_format" cutoff="1.0"/>)"));
    CHECK_THROWS_MATCHES(ecp.buildL2(doc.getRoot()), UniformCommunicateError,
                         Catch::Matchers::Message("ECPComponentBuilder::buildL2: Unrecognized format \"unknown_format\" in PP file."));
  }

  SECTION("missing or negative cutoff")
  {
    ecp.grid_global = std::make_unique<LinearGrid<ECPComponentBuilder::mRealType>>();
    ecp.grid_global->set(0.0, 1.0, 10);
    REQUIRE(doc.parseFromString(R"(<L2 format="V" cutoff="-0.5"/>)"));
    CHECK_THROWS_MATCHES(ecp.buildL2(doc.getRoot()), UniformCommunicateError,
                         Catch::Matchers::Message("ECPComponentBuilder::buildL2: L2 attribute cutoff is missing or negative. Cutoff must be a positive real number."));
  }
}

TEST_CASE("ECPComponentBuilder buildLocal exceptions", "[hamiltonian][ecp]")
{
  Communicate& c(*OHMMS::Controller);
  ECPComponentBuilder ecp("test_ecp_local", c);
  Libxml2Document doc;

  SECTION("missing grid information")
  {
    REQUIRE(doc.parseFromString(R"(<local/>)"));
    CHECK_THROWS_MATCHES(ecp.buildLocal(doc.getRoot()), UniformCommunicateError,
                         Catch::Matchers::Message("ECPComponentBuilder::buildLocal Missing grid information. "));
  }

  SECTION("illegal local pseudopotential")
  {
    ecp.grid_global = std::make_unique<LinearGrid<ECPComponentBuilder::mRealType>>();
    ecp.grid_global->set(0.0, 1.0, 10);
    ecp.Zeff = 4.0;
    REQUIRE(doc.parseFromString(R"(<local>
      <basisGroup type="Gaussian" n="0" l="0">
        <radfunc exponent="1.0" contraction="0.0"/>
      </basisGroup>
    </local>)"));
    CHECK_THROWS_MATCHES(ecp.buildLocal(doc.getRoot()), UniformCommunicateError,
                         Catch::Matchers::Message("ECPComponentBuilder::buildLocal. Illegal Local Pseudopotential"));
  }
}

TEST_CASE("ECPComponentBuilder buildSemiLocalAndLocal exceptions", "[hamiltonian][ecp]")
{
  Communicate& c(*OHMMS::Controller);
  Libxml2Document doc1, doc2;

  SECTION("no global grid")
  {
    ECPComponentBuilder ecp("test_ecp_semilocal", c, 4, 1);
    REQUIRE(doc1.parseFromString(R"(<semilocal/>)"));
    std::vector<xmlNodePtr> semiPtr = {doc1.getRoot()};
    CHECK_THROWS_MATCHES(ecp.buildSemiLocalAndLocal(semiPtr), UniformCommunicateError,
                         Catch::Matchers::Message("ECPComponentBuilder::buildSemiLocalAndLocal. Global grid needs to be defined.\n"));
  }

  SECTION("more than one semilocal tag")
  {
    ECPComponentBuilder ecp("test_ecp_semilocal", c, 4, 1);
    ecp.grid_global = std::make_unique<LinearGrid<ECPComponentBuilder::mRealType>>();
    ecp.grid_global->set(0.0, 1.0, 10);
    REQUIRE(doc1.parseFromString(R"(<semilocal/>)"));
    REQUIRE(doc2.parseFromString(R"(<semilocal/>)"));
    std::vector<xmlNodePtr> semiPtr = {doc1.getRoot(), doc2.getRoot()};
    CHECK_THROWS_MATCHES(ecp.buildSemiLocalAndLocal(semiPtr), UniformCommunicateError,
                         Catch::Matchers::Message("ECPComponentBuilder::buildSemiLocalAndLocal. We have more than one semilocal sections in the PP xml file"));
  }

  SECTION("default nrule not determined")
  {
    ECPComponentBuilder ecp("test_ecp_semilocal", c);
    ecp.grid_global = std::make_unique<LinearGrid<ECPComponentBuilder::mRealType>>();
    ecp.grid_global->set(0.0, 1.0, 10);
    ecp.pp_nonloc = std::make_unique<NonLocalECPComponent>();
    ecp.pp_nonloc->setLmax(12);
    REQUIRE(doc1.parseFromString(R"(<semilocal/>)"));
    std::vector<xmlNodePtr> semiPtr = {doc1.getRoot()};
    CHECK_THROWS_MATCHES(ecp.buildSemiLocalAndLocal(semiPtr), UniformCommunicateError,
                         Catch::Matchers::Message("Default value for pseudopotential nrule not determined."));
  }

  SECTION("unrecognized format")
  {
    ECPComponentBuilder ecp("test_ecp_semilocal", c, 4, 1);
    ecp.grid_global = std::make_unique<LinearGrid<ECPComponentBuilder::mRealType>>();
    ecp.grid_global->set(0.0, 1.0, 10);
    REQUIRE(doc1.parseFromString(R"(<semilocal format="bad_fmt"/>)"));
    std::vector<xmlNodePtr> semiPtr = {doc1.getRoot()};
    CHECK_THROWS_MATCHES(ecp.buildSemiLocalAndLocal(semiPtr), UniformCommunicateError,
                         Catch::Matchers::Message("ECPComponentBuilder::buildSemiLocalAndLocal.Unrecognized format \"bad_fmt\" in PP file.\n"));
  }

  SECTION("unrecognized angular momentum in vps")
  {
    ECPComponentBuilder ecp("test_ecp_semilocal", c, 4, 1);
    ecp.grid_global = std::make_unique<LinearGrid<ECPComponentBuilder::mRealType>>();
    ecp.grid_global->set(0.0, 1.0, 10);
    REQUIRE(doc1.parseFromString(R"(<semilocal format="V">
      <vps l="unknown"/>
    </semilocal>)"));
    std::vector<xmlNodePtr> semiPtr = {doc1.getRoot()};
    CHECK_THROWS_MATCHES(ecp.buildSemiLocalAndLocal(semiPtr), UniformCommunicateError,
                         Catch::Matchers::Message("ECPComponentBuilder::buildSemiLocalAndLocal. Requested angular momentum unknown not available.\n"));
  }

  SECTION("unrecognized angular momentum in vps_so")
  {
    ECPComponentBuilder ecp("test_ecp_semilocal", c, 4, 1);
    ecp.grid_global = std::make_unique<LinearGrid<ECPComponentBuilder::mRealType>>();
    ecp.grid_global->set(0.0, 1.0, 10);
    REQUIRE(doc1.parseFromString(R"(<semilocal format="V">
      <vps_so l="unknown"/>
    </semilocal>)"));
    std::vector<xmlNodePtr> semiPtr = {doc1.getRoot()};
    CHECK_THROWS_MATCHES(ecp.buildSemiLocalAndLocal(semiPtr), UniformCommunicateError,
                         Catch::Matchers::Message("ECPComponentBuilder::buildSemiLocalAndLocal. Requested angular momentum unknown not available for SO.\n"));
  }

  SECTION("local channel not specified when multiple vps")
  {
    ECPComponentBuilder ecp("test_ecp_semilocal", c, 4, -1);
    ecp.grid_global = std::make_unique<LinearGrid<ECPComponentBuilder::mRealType>>();
    ecp.grid_global->set(0.0, 1.0, 10);
    REQUIRE(doc1.parseFromString(R"(<semilocal format="V">
      <vps l="s"/>
      <vps l="p"/>
    </semilocal>)"));
    std::vector<xmlNodePtr> semiPtr = {doc1.getRoot()};
    CHECK_THROWS_MATCHES(ecp.buildSemiLocalAndLocal(semiPtr), UniformCommunicateError,
                         Catch::Matchers::Message("ECPComponentBuilder::doBreakUp: The local channel is specified in neither the pseudopotential file nor the "
                                                  "input file.\nPlease add 'l-local=\"n\"' attribute to either file.\n"));
  }

  SECTION("so channels count mismatch")
  {
    ECPComponentBuilder ecp("test_ecp_semilocal", c, 4, 0);
    ecp.grid_global = std::make_unique<LinearGrid<ECPComponentBuilder::mRealType>>();
    ecp.grid_global->set(0.0, 1.0, 10);
    REQUIRE(doc1.parseFromString(R"(<semilocal format="V" npots-so="2">
      <vps l="s"/>
    </semilocal>)"));
    std::vector<xmlNodePtr> semiPtr = {doc1.getRoot()};
    CHECK_THROWS_MATCHES(ecp.buildSemiLocalAndLocal(semiPtr), UniformCommunicateError,
                         Catch::Matchers::Message("Error. npots-so=0 while declared number of SO channels is 2\n"));
  }
}
} // namespace qmcplusplus
