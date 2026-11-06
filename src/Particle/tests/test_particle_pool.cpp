//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by:  Mark Dewing, mdewing@anl.gov Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_test_macros.hpp>
#include "Utilities/for_testing/Catch2Approx.h"


#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Particle/ParticleSetPool.h"
#include "Particle/LongRange/StructFact.h"
#include "Particle/tests/MinimalParticlePool.h"


#include <stdio.h>
#include <string>
#include <sstream>
#include <cmath>
#include <limits>
#include <type_traits>


namespace qmcplusplus
{
TEST_CASE("ParticleSetPool", "[qmcapp]")
{
  Communicate& c(*OHMMS::Controller);

  ParticleSetPool pp(c);

  // See ParticleIO/tests/test_xml_io.cpp for particle parsing
  const char* particles = R"(
<particleset name="ion0" size="1">
  <group name="He">
    <parameter name="charge">2</parameter>
  </group>
  <attrib name="position" datatype="posArray">
    0.1 0.2 0.3
  </attrib>
</particleset>
)";
  Libxml2Document doc;
  REQUIRE(doc.parseFromString(particles));

  xmlNodePtr root = doc.getRoot();

  pp.put(root);

  ParticleSet* ions = pp.getParticleSet("ion0");
  REQUIRE(ions != NULL);

  ParticleSet* not_here = pp.getParticleSet("does_not_exist");
  REQUIRE(not_here == NULL);

  MCWalkerConfiguration* ws = pp.getWalkerSet("ion0");
  REQUIRE(ws != NULL);

  REQUIRE_THROWS_AS(pp.getWalkerSet("does_not_exist"), std::runtime_error);

  auto ps2 = std::make_unique<ParticleSet>(pp.getSimulationCell());
  ps2->setName("particle_set_2");
  pp.addParticleSet(std::move(ps2));

  // should do nothing, since no random particlesets were specified
  pp.randomize();

  std::stringstream out;
  pp.get(out);
  //app_log() << "ParticleSetPool::get returns  " << std::endl;
  //app_log() << out.str() << std::endl;
}

TEST_CASE("ParticleSetPool random", "[qmcapp]")
{
  Communicate& c(*OHMMS::Controller);

  ParticleSetPool pp(c);

  // See ParticleIO/tests/test_xml_io.cpp for particle parsing
  const char* particles = R"(
<tmp>
<particleset name="ion0" size="1">
  <group name="He">
    <parameter name="charge">2</parameter>
  </group>
  <attrib name="position" datatype="posArray">
    0.1 0.2 0.3
  </attrib>
</particleset>
<particleset name="elec" random="yes" randomsrc="ion0" spinor="yes">
   <group name="u" size="4">
      <parameter name="charge">-1</parameter>
   </group>
</particleset>
</tmp>
)";
  Libxml2Document doc;
  REQUIRE(doc.parseFromString(particles));

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr part_ion = xmlFirstElementChild(root);
  pp.put(part_ion);
  pp.put(xmlNextElementSibling(part_ion));

  ParticleSet* ions = pp.getParticleSet("ion0");
  REQUIRE(ions != NULL);
  REQUIRE(!ions->isSpinor());

  ParticleSet* elec = pp.getParticleSet("elec");
  REQUIRE(ions != NULL);
  REQUIRE(elec->isSpinor());
  REQUIRE(elec->R.size() == 4);
  REQUIRE(elec->spins.size() == 4);


  // should do something
  pp.randomize();

  REQUIRE(elec->R[0][0] != 0.0);
  REQUIRE(elec->spins[0] != 0.0);
}

TEST_CASE("MinimalParticlePool factory reset computes ionic structure factors", "[qmcapp][particle]")
{
  auto pool = MinimalParticlePool::make_diamondC_1x1x1(*OHMMS::Controller);
  const ParticleSet* ions_ptr = pool.getParticleSet("ion");
  REQUIRE(ions_ptr != nullptr);
  const ParticleSet& ions = *ions_ptr;
  REQUIRE(ions.hasSK());
  const auto& sk = ions.getSK();
  const auto& kpts = ions.getSimulationCell().getKLists().getKptsCart();

  REQUIRE(ions.getTotalNum() == 2);
  REQUIRE(sk.rhok_r.size(0) == ions.groups());
  REQUIRE(sk.rhok_r.size(1) == kpts.size());
  REQUIRE(sk.rhok_i.size(0) == ions.groups());
  REQUIRE(sk.rhok_i.size(1) == kpts.size());
  const double tolerance = 100 * std::numeric_limits<StructFact::RealType>::epsilon();

  // The first nonzero k vector avoids passing on an uninitialized zero-k density.
  bool checked_nonzero_k = false;
  for (std::size_t ik = 0; ik < kpts.size(); ++ik)
  {
    const auto& k = kpts[ik];
    if (dot(k, k) == 0.0)
      continue;

    for (int species = 0; species < sk.rhok_r.size(0); ++species)
    {
      double expected_r = 0.0;
      double expected_i = 0.0;
      for (int particle = 0; particle < ions.getTotalNum(); ++particle)
        if (ions.GroupID[particle] == species)
        {
          const double phase = dot(k, ions.R[particle]);
          expected_r += std::cos(phase);
          expected_i += std::sin(phase);
        }
      CHECK(sk.rhok_r(species, ik) == Approx(expected_r).margin(tolerance));
      CHECK(sk.rhok_i(species, ik) == Approx(expected_i).margin(tolerance));
    }
    checked_nonzero_k = true;
    break;
  }
  REQUIRE(checked_nonzero_k);
}

TEST_CASE("ParticleSetPool putLattice", "[qmcapp]")
{
  Communicate& c(*OHMMS::Controller);

  ParticleSetPool pp(c);

  const char* lattice = R"(<parameter name="lattice"> </parameter>)";

  Libxml2Document doc;
  REQUIRE(doc.parseFromString(lattice));

  xmlNodePtr root = doc.getRoot();
  pp.readSimulationCellXML(root);
}
} // namespace qmcplusplus
