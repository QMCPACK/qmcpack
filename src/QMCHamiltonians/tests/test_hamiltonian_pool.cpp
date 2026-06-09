//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, mdewing@anl.gov Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"


#include "Configuration.h"
#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCHamiltonians/HamiltonianPool.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/WaveFunctionPool.h"
#include "Utilities/RuntimeOptions.h"

#include <stdio.h>
#include <string>
#include <sstream>


namespace qmcplusplus
{
extern std::unique_ptr<ParticleSet> createElectronParticleSet(const SimulationCell& simulation_cell);

TEST_CASE("HamiltonianPool", "[qmcapp]")
{
  Communicate* c;
  c = OHMMS::Controller;

  // See src/QMCHamiltonians/tests/test_hamiltonian_factory for parsing tests
  const char* hamiltonian_xml = R"(<hamiltonian name="h0" type="generic" target="e">
         <pairpot type="coulomb" name="ElecElec" source="e" target="e"/>
</hamiltonian>)";

  Libxml2Document doc;
  REQUIRE(doc.parseFromString(hamiltonian_xml));

  xmlNodePtr root = doc.getRoot();

  ParticleSetPool pp(c);
  auto qp = createElectronParticleSet(pp.getSimulationCell());
  pp.addParticleSet(std::move(qp));

  RuntimeOptions runtime_options;
  WaveFunctionPool wfp(runtime_options, pp, c);
  wfp.add("psi0", std::make_unique<TrialWaveFunction>(runtime_options, "psi0"));

  HamiltonianPool hpool(pp, wfp, c);

  REQUIRE(hpool.empty());

  hpool.put(root);

  // test contains()
  REQUIRE(hpool.contains("h0"));
  REQUIRE(!hpool.contains("h1"));

  // test getHamiltonian()
  QMCHamiltonian& ham(hpool.getHamiltonian().value());
  // Bare kinetic energy is always added
  REQUIRE(ham.size() == 2);

  auto ham_noname_optional = hpool.getHamiltonian();
  REQUIRE(ham_noname_optional);
  QMCHamiltonian& ham_noname(*ham_noname_optional);
  REQUIRE(&ham == &ham_noname);
  auto ham_empty_optional = hpool.getHamiltonian("");
  REQUIRE(ham_empty_optional);
  QMCHamiltonian& ham_empty(*ham_empty_optional);
  REQUIRE(&ham == &ham_empty);
}

} // namespace qmcplusplus
