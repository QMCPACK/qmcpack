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

TEST_CASE("read_particle_mass_same_xml", "[particle_io][xml]")
{
  // test that particle masses are properly read in

  OHMMS::Controller->initialize(0, NULL);
  Communicate *c = OHMMS::Controller;

const char *particles =
"<tmp> \
<particleset name=\"e\" random=\"yes\"> \
   <group name=\"u\" size=\"4\" mass=\"1.0\">\
      <parameter name=\"charge\"              >    -1                    </parameter>\
      <parameter name=\"mass\"                >    1.0                   </parameter>\
   </group>\
   <group name=\"d\" size=\"4\" mass=\"1.0\">\
      <parameter name=\"charge\"              >    -1                    </parameter>\
      <parameter name=\"mass\"                >    1.0                   </parameter>\
   </group>\
</particleset>\
<particleset name=\"ion0\">\
   <group name=\"H\" size=\"8\" mass=\"1836.15\">\
      <parameter name=\"charge\"              >    1                     </parameter>\
      <parameter name=\"valence\"             >    1                     </parameter>\
      <parameter name=\"atomicnumber\"        >    1                     </parameter>\
      <parameter name=\"mass\"                >    1836.15               </parameter>\
      <attrib name=\"position\" datatype=\"posArray\" condition=\"0\">\
               0.00000000        0.00000000        0.00000000\
               2.11170946        0.00000000        0.00000000\
               0.00000000        2.11170946        0.00000000\
               2.11170946        2.11170946        0.00000000\
               0.00000000        0.00000000        2.11170946\
               2.11170946        0.00000000        2.11170946\
               0.00000000        2.11170946        2.11170946\
               2.11170946        2.11170946        2.11170946\
      </attrib>\
   </group>\
</particleset>\
</tmp> \
"; // simple cubic lattice at rs=1.31

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root  = doc.getRoot();
  xmlNodePtr part1 = xmlFirstElementChild(root);
  xmlNodePtr part2 = xmlNextElementSibling(part1);

  Tensor<int, 3> tmat; // assuming OHMMSDIM==3
  tmat(0,0) = 1;
  tmat(1,1) = 1;
  tmat(2,2) = 1;

  ParticleSet ions, electrons;

  XMLParticleParser parse_electrons(electrons, tmat);
  parse_electrons.put(part1);
  REQUIRE(electrons.getName() == "e");

  XMLParticleParser parse_ions(ions, tmat);
  parse_ions.put(part2);
  REQUIRE(ions.getName() == "ion0");

  REQUIRE( ions.SameMass );
  REQUIRE( electrons.SameMass );

  // test electrons
  SpeciesSet& tspecies( electrons.getSpeciesSet() );
  int massind=tspecies.addAttribute("mass");
  char order[]="uuuudddd";
  for (int iat=0;iat<electrons.getTotalNum(); iat++){
    int species_id = electrons.GroupID[iat];
    std::string species_name = tspecies.speciesName[species_id];
    REQUIRE(*species_name.c_str() == order[iat]);
    REQUIRE(tspecies(massind,species_id) == Approx(1.0));
  }

  // test ions
  SpeciesSet& pspecies( ions.getSpeciesSet() );
  int pmassind=pspecies.addAttribute("mass");
  char porder[]="HHHHHHHH";
  for (int iat=0;iat<ions.getTotalNum(); iat++){
    int species_id = ions.GroupID[iat];
    std::string species_name = pspecies.speciesName[species_id];
    REQUIRE(*species_name.c_str() == porder[iat]);
    REQUIRE(pspecies(pmassind,species_id) == Approx(1836.15));
  }
} // TEST_CASE read_particle_mass_same_xml

} // namespace qmcplusplus
