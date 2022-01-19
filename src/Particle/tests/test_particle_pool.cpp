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


#include "catch.hpp"


#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Particle/ParticleSetPool.h"


#include <stdio.h>
#include <string>
#include <sstream>


namespace qmcplusplus
{
TEST_CASE("ParticleSetPool", "[qmcapp]")
{
  Communicate* c;
  c = OHMMS::Controller;

  ParticleSetPool pp(c);

  // See ParticleIO/tests/test_xml_io.cpp for particle parsing
  const char* particles = " \
<particleset name=\"ion0\" size=\"1\"> \
  <group name=\"He\"> \
    <parameter name=\"charge\">2</parameter> \
  </group> \
  <attrib name=\"position\" datatype=\"posArray\"> \
    0.1 0.2 0.3 \
  </attrib> \
</particleset> \
";
  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  pp.put(root);

  ParticleSet* ions = pp.getParticleSet("ion0");
  REQUIRE(ions != NULL);

  ParticleSet* not_here = pp.getParticleSet("does_not_exist");
  REQUIRE(not_here == NULL);

  ParticleSet* ws = pp.getWalkerSet("ion0");
  REQUIRE(ws != NULL);

  auto ps2 = std::make_unique<ParticleSet>(pp.getSimulationCell());
  ps2->setName("particle_set_2");
  pp.addParticleSet(std::move(ps2));

  // should do nothing, since no random particlesets were specified
  pp.randomize();

  std::stringstream out;
  pp.get(out);
  //std::cout << "ParticleSetPool::get returns  " << std::endl;
  //std::cout << out.str() << std::endl;
}

TEST_CASE("ParticleSetPool random", "[qmcapp]")
{
  Communicate* c;
  c = OHMMS::Controller;

  ParticleSetPool pp(c);

  // See ParticleIO/tests/test_xml_io.cpp for particle parsing
  const char* particles = " \
<tmp> \
<particleset name=\"ion0\" size=\"1\"> \
  <group name=\"He\"> \
    <parameter name=\"charge\">2</parameter> \
  </group> \
  <attrib name=\"position\" datatype=\"posArray\"> \
    0.1 0.2 0.3 \
  </attrib> \
</particleset> \
<particleset name=\"elec\" random=\"yes\" randomsrc=\"ion0\" spinor=\"yes\"> \
   <group name=\"u\" size=\"4\"> \
      <parameter name=\"charge\">-1</parameter> \
   </group> \
</particleset> \
</tmp> \
";
  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

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
#if !defined(QMC_CUDA)
  REQUIRE(elec->spins[0] != 0.0);
#endif
}

TEST_CASE("ParticleSetPool putLattice", "[qmcapp]")
{
  Communicate* c;
  c = OHMMS::Controller;

  ParticleSetPool pp(c);

  const char* lattice = "<parameter name='lattice'> </parameter>";

  Libxml2Document doc;
  bool okay = doc.parseFromString(lattice);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();
  pp.readSimulationCellXML(root);
}
} // namespace qmcplusplus
