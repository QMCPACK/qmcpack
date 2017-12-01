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



#include "Message/catch_mpi_main.hpp"


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

  OHMMS::Controller->initialize(0, NULL);
  Communicate *c = OHMMS::Controller;


const char *particles = \
"<tmp> \
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

  Tensor<int, 3> tmat; // assuming OHMMSDIM==3
  tmat(0,0) = 1;
  tmat(1,1) = 1;
  tmat(2,2) = 1;

  ParticleSet ions, electrons;
  XMLParticleParser parse_ions(ions, tmat);
  xmlNodePtr part1 = xmlFirstElementChild(root);
  parse_ions.put(part1);

  REQUIRE(ions.groups() == 1);
  REQUIRE(ions.R.size() == 1);
  REQUIRE(ions.R[0][0] == Approx(0.1));
  REQUIRE(ions.R[0][1] == Approx(0.2));
  REQUIRE(ions.R[0][2] == Approx(0.3));
  REQUIRE(ions.getName() == "ion0");

  XMLParticleParser parse_electrons(electrons, tmat);
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
}


