//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"


#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/Tensor.h"
#include "Particle/ParticleSet.h"
#include "ParticleIO/XMLParticleIO.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("read_particleset_xml", "[particle_io][xml]")
{
  const char* particles = "<tmp> \
<particleset name=\"ion0\" size=\"1\"> \
  <group name=\"He\"> \
    <parameter name=\"charge\">2</parameter> \
  </group> \
  <attrib name=\"position\" datatype=\"posArray\"> \
    0.1 0.2 0.3 \
  </attrib> \
</particleset> \
<particleset name=\"e\"> \
  <group name=\"u\" size=\"1\"> \
    <parameter name=\"charge\">-1</parameter> \
    <attrib name=\"position\" datatype=\"posArray\"> \
      -0.28   0.0225     -2.709 \
    </attrib> \
  </group> \
  <group name=\"d\" size=\"1\"> \
    <parameter name=\"charge\">-1</parameter> \
    <attrib name=\"position\" datatype=\"posArray\"> \
      -1.08389   1.9679     -0.0128914 \
    </attrib> \
  </group> \
</particleset> \
</tmp> \
";
  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell), electrons(simulation_cell);

  XMLParticleParser parse_ions(ions);
  xmlNodePtr part1 = xmlFirstElementChild(root);
  parse_ions.put(part1);

  REQUIRE(ions.groups() == 1);
  REQUIRE(ions.R.size() == 1);
  REQUIRE(ions.R[0][0] == Approx(0.1));
  REQUIRE(ions.R[0][1] == Approx(0.2));
  REQUIRE(ions.R[0][2] == Approx(0.3));
  REQUIRE(ions.getName() == "ion0");

  XMLParticleParser parse_electrons(electrons);
  xmlNodePtr part2 = xmlNextElementSibling(part1);
  parse_electrons.put(part2);

  REQUIRE(electrons.groups() == 2);
  REQUIRE(electrons.R.size() == 2);
  REQUIRE(electrons.R[0][0] == Approx(-0.28));
  REQUIRE(electrons.R[0][1] == Approx(0.0225));
  REQUIRE(electrons.R[0][2] == Approx(-2.709));

  REQUIRE(electrons.R[1][0] == Approx(-1.08389));
  REQUIRE(electrons.R[1][1] == Approx(1.9679));
  REQUIRE(electrons.R[1][2] == Approx(-0.0128914));
  REQUIRE(electrons.getName() == "e");
}

TEST_CASE("read_dynamic_spin_eset_xml", "[particle_io][xml]")
{
  const char* particles = "<tmp> \
<particleset name=\"e\"> \
  <group name=\"e\" size=\"3\"> \
    <parameter name=\"charge\">-1</parameter> \
    <attrib name=\"position\" datatype=\"posArray\"> \
      -0.28   0.0225     -2.709 \
      -1.28   1.0225     -1.709 \
      -2.28   2.0225     -0.709 \
    </attrib> \
    <attrib name=\"spins\" datatype=\"scalarArray\"> \
      1.0 \
      0.2 \
      3.0 \
    </attrib> \
  </group> \
</particleset> \
</tmp> \
";
  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr part1 = xmlFirstElementChild(root);

  const SimulationCell simulation_cell;
  ParticleSet electrons(simulation_cell);

  XMLParticleParser parse_electrons(electrons);
  parse_electrons.put(part1);

  REQUIRE(electrons.groups() == 1);

  REQUIRE(electrons.R.size() == 3);

  REQUIRE(electrons.R[0][0] == Approx(-0.28));
  REQUIRE(electrons.R[0][1] == Approx(0.0225));
  REQUIRE(electrons.R[0][2] == Approx(-2.709));

  REQUIRE(electrons.R[1][0] == Approx(-1.28));
  REQUIRE(electrons.R[1][1] == Approx(1.0225));
  REQUIRE(electrons.R[1][2] == Approx(-1.709));

  REQUIRE(electrons.R[2][0] == Approx(-2.28));
  REQUIRE(electrons.R[2][1] == Approx(2.0225));
  REQUIRE(electrons.R[2][2] == Approx(-0.709));

  REQUIRE(electrons.spins[0] == Approx(1.0));
  REQUIRE(electrons.spins[1] == Approx(0.2));
  REQUIRE(electrons.spins[2] == Approx(3.0));

  REQUIRE(electrons.getName() == "e");
}
} // namespace qmcplusplus
