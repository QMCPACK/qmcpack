//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"


#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/Tensor.h"
#include "Particle/ParticleSet.h"
#include "ParticleIO/XMLParticleIO.h"
#include "ParticleIO/ParticleLayoutIO.h"
#include "Particle/DistanceTableData.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{

TEST_CASE("distance_open_z", "[distance_table][xml]")
{
  // test that particle distances are properly calculated

  OHMMS::Controller->initialize(0, NULL);
  Communicate *c = OHMMS::Controller;

const char *particles =
"<tmp> \
<particleset name=\"e\" random=\"yes\"> \
   <group name=\"u\" size=\"1\" mass=\"1.0\">\
      <parameter name=\"charge\"              >    -1                    </parameter>\
      <parameter name=\"mass\"                >    1.0                   </parameter>\
      <attrib name=\"position\" datatype=\"posArray\" condition=\"0\">\
               0.00000000        0.00000000        0.20000000\
      </attrib>\
   </group>\
   <group name=\"d\" size=\"1\" mass=\"1.0\">\
      <parameter name=\"charge\"              >    -1                    </parameter>\
      <parameter name=\"mass\"                >    1.0                   </parameter>\
      <attrib name=\"position\" datatype=\"posArray\" condition=\"0\">\
               0.00000000        0.00000000       -0.20000000\
      </attrib>\
   </group>\
</particleset>\
<particleset name=\"ion0\">\
   <group name=\"H\" size=\"2\" mass=\"1836.15\">\
      <parameter name=\"charge\"              >    1                     </parameter>\
      <parameter name=\"valence\"             >    1                     </parameter>\
      <parameter name=\"atomicnumber\"        >    1                     </parameter>\
      <parameter name=\"mass\"                >    1836.15               </parameter>\
      <attrib name=\"position\" datatype=\"posArray\" condition=\"0\">\
               0.00000000        0.00000000        0.00000000\
               0.00000000        0.00000000        0.50000000\
      </attrib>\
   </group>\
</particleset>\
</tmp> \
";

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

  // read particle set
  ParticleSet ions, electrons;

  XMLParticleParser parse_electrons(electrons, tmat);
  parse_electrons.put(part1);

  XMLParticleParser parse_ions(ions, tmat);
  parse_ions.put(part2);

  REQUIRE( electrons.getName() == "e" );
  REQUIRE( ions.getName() == "ion0"   );
  REQUIRE( ions.SameMass      );
  REQUIRE( electrons.SameMass );

  // calculate particle distances
  electrons.addTable(ions,DT_AOS);
  electrons.update();

  // get target particle set's distance table data
  int tid = electrons.getTable(ions); //this is bad
  DistanceTableData* dtable = electrons.DistTables[tid];
  REQUIRE(dtable->getName() == "ion0_e");

  double expect[] = {0.2,0.2,0.3,0.7};
  int idx(0);
  for (int iat=0;iat<ions.getTotalNum();iat++){
    for (int jat=0;jat<electrons.getTotalNum();jat++,idx++){
      // note: target particle set is special (electrons in this case)
      // int tid = target_pset.getTable(source_pset)
      // DistanceTableData* dtable = target_pset.DistTables[tid] 
      // dtable->loc(source_ptcl_idx,target_ptcl_idx) !! source first target second !?
      double dist = dtable->r(dtable->loc(iat,jat));
      REQUIRE( dist == Approx(expect[idx]) );
    }
  }

  // get distance between target="e" group="u" iat=0 and source="ion0" group="H" jat=1

} // TEST_CASE distance_open_z

TEST_CASE("distance_open_xy", "[distance_table][xml]")
{

  OHMMS::Controller->initialize(0, NULL);
  Communicate *c = OHMMS::Controller;

const char *particles =
"<tmp> \
<particleset name=\"e\" random=\"yes\"> \
   <group name=\"u\" size=\"2\" mass=\"1.0\">\
      <parameter name=\"charge\"              >    -1                    </parameter>\
      <parameter name=\"mass\"                >    1.0                   </parameter>\
      <attrib name=\"position\" datatype=\"posArray\" condition=\"0\">\
               0.70000000        0.00000000        0.00000000\
               0.00000000        1.00000000        0.00000000\
      </attrib>\
   </group>\
   <group name=\"d\" size=\"1\" mass=\"1.0\">\
      <parameter name=\"charge\"              >    -1                    </parameter>\
      <parameter name=\"mass\"                >    1.0                   </parameter>\
      <attrib name=\"position\" datatype=\"posArray\" condition=\"0\">\
              -0.90000000        0.00000000        0.00000000\
      </attrib>\
   </group>\
</particleset>\
<particleset name=\"ion0\">\
   <group name=\"H\" size=\"2\" mass=\"1836.15\">\
      <parameter name=\"charge\"              >    1                     </parameter>\
      <parameter name=\"valence\"             >    1                     </parameter>\
      <parameter name=\"atomicnumber\"        >    1                     </parameter>\
      <parameter name=\"mass\"                >    1836.15               </parameter>\
      <attrib name=\"position\" datatype=\"posArray\" condition=\"0\">\
               0.00000000        0.00000000        0.00000000\
               1.00000000        0.00000000        0.00000000\
      </attrib>\
   </group>\
</particleset>\
</tmp> \
";

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

  // read particle set
  ParticleSet ions, electrons;

  XMLParticleParser parse_electrons(electrons, tmat);
  parse_electrons.put(part1);

  XMLParticleParser parse_ions(ions, tmat);
  parse_ions.put(part2);

  REQUIRE( electrons.getName() == "e" );
  REQUIRE( ions.getName() == "ion0"   );
  REQUIRE( ions.SameMass      );
  REQUIRE( electrons.SameMass );

  // calculate particle distances
  electrons.addTable(ions,DT_AOS);
  electrons.update();

  // get distance table attached to target particle set (electrons)
  int tid = electrons.getTable(ions);
  DistanceTableData* dtable = electrons.DistTables[tid];
  REQUIRE(dtable->getName() == "ion0_e");

  // calculate distance, one source particle at a time i.e.
  // H0 - e0: 0.7
  // H0 - e1: 1.0
  // H0 - e2: 0.9
  // H1 - e0: 0.3
  // etc.
  double expect[] = {
    0.7,1.0,0.9
   ,0.3,std::sqrt(2),1.9
  };
  int idx(0);
  for (int iat=0;iat<ions.getTotalNum();iat++){
    for (int jat=0;jat<electrons.getTotalNum();jat++,idx++){
      double dist = dtable->r(dtable->loc(iat,jat));
      REQUIRE( dist == Approx(expect[idx]) );
    }
  }

} // TEST_CASE distance_open_xy

TEST_CASE("distance_open_species_deviation", "[distance_table][xml]")
{
  // pull out distances between specific species

  OHMMS::Controller->initialize(0, NULL);
  Communicate *c = OHMMS::Controller;

const char *particles =
"<tmp> \
<particleset name=\"e\" random=\"yes\"> \
   <group name=\"u\" size=\"2\" mass=\"1.0\">\
      <parameter name=\"charge\"              >    -1                    </parameter>\
      <parameter name=\"mass\"                >    1.0                   </parameter>\
      <attrib name=\"position\" datatype=\"posArray\" condition=\"0\">\
               0.70000000        0.00000000        0.00000000\
               0.00000000        1.00000000        0.00000000\
      </attrib>\
   </group>\
   <group name=\"d\" size=\"1\" mass=\"1.0\">\
      <parameter name=\"charge\"              >    -1                    </parameter>\
      <parameter name=\"mass\"                >    1.0                   </parameter>\
      <attrib name=\"position\" datatype=\"posArray\" condition=\"0\">\
              -0.90000000        0.00000000        0.00000000\
      </attrib>\
   </group>\
</particleset>\
<particleset name=\"ion0\">\
   <group name=\"H\" size=\"2\" mass=\"1836.15\">\
      <parameter name=\"charge\"              >    1                     </parameter>\
      <parameter name=\"valence\"             >    1                     </parameter>\
      <parameter name=\"atomicnumber\"        >    1                     </parameter>\
      <parameter name=\"mass\"                >    1836.15               </parameter>\
      <attrib name=\"position\" datatype=\"posArray\" condition=\"0\">\
               0.00000000        0.00000000        0.00000000\
               1.00000000        0.00000000        0.00000000\
      </attrib>\
   </group>\
</particleset>\
</tmp> \
";

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

  // read particle set
  ParticleSet ions, electrons;

  XMLParticleParser parse_electrons(electrons, tmat);
  parse_electrons.put(part1);

  XMLParticleParser parse_ions(ions, tmat);
  parse_ions.put(part2);

  REQUIRE( electrons.getName() == "e" );
  REQUIRE( ions.getName() == "ion0"   );
  REQUIRE( ions.SameMass      );
  REQUIRE( electrons.SameMass );

  // calculate particle distances
  electrons.addTable(ions,DT_AOS);
  electrons.update();

  // get distance table attached to target particle set (electrons)
  int tid = electrons.getTable(ions);
  DistanceTableData* dtable = electrons.DistTables[tid];
  REQUIRE(dtable->getName() == "ion0_e");

  // get the electron species set
  SpeciesSet& especies(electrons.getSpeciesSet());

  // !! assume "u" and "H" groups have the same number of particles
  double latdev2 = 0.0; // mean-squared deviation from lattice
  int cur_jat(-1); // keep an index to the last found target particle
  double expect[] = {0.7,std::sqrt(2)};
  int idx(0);
  for (int iat=0;iat<ions.getTotalNum();iat++,idx++){
    for (int jat=cur_jat+1;jat<electrons.getTotalNum();jat++,idx++){
      // find next "u"
      int species_id = electrons.GroupID[jat];
      string species_name = especies.speciesName[species_id];
      if (species_name != "u") continue;

      // calculate distance from lattice site iat
      double dist = dtable->r(dtable->loc(iat,jat));
      latdev2 += std::pow(dist,2); // !? pow(x,2) does what?
      REQUIRE( dist == Approx(expect[idx]) );
      cur_jat = jat;
      break;
    }
  }
  REQUIRE( latdev2/ions.getTotalNum() == Approx(1.245) );

} // TEST_CASE distance_open_species_deviation

TEST_CASE("distance_pbc_z", "[distance_table][xml]")
{
  // test that particle distances are properly calculated under periodic boundary condition
  // There are many details in this example, but the main idea is simple: When a particle is moved by a full lattice vector, no distance should change.

  OHMMS::Controller->initialize(0, NULL);
  Communicate *c = OHMMS::Controller;

const char *particles =
"<tmp> \
  <simulationcell>\
     <parameter name=\"lattice\" units=\"bohr\">\
              6.00000000        0.00000000        0.00000000\
              0.00000000        6.00000000        0.00000000\
              0.00000000        0.00000000        6.00000000\
     </parameter>\
     <parameter name=\"bconds\">\
        p p p\
     </parameter>\
     <parameter name=\"LR_dim_cutoff\"       >    15                 </parameter>\
  </simulationcell>\
  <particleset name=\"e\">\
     <group name=\"u\" size=\"4\" mass=\"1.0\">\
        <parameter name=\"charge\"              >    -1                    </parameter>\
        <parameter name=\"mass\"                >    1.0                   </parameter>\
        <attrib name=\"position\" datatype=\"posArray\" condition=\"0\">\
                 0.00000000        0.00000000        0.00000000\
                 3.00000000        0.00000000        0.00000000\
                 0.00000000        3.00000000        0.00000000\
                 3.00000000        3.00000000        0.00000000\
        </attrib>\
     </group>\
     <group name=\"d\" size=\"4\" mass=\"1.0\">\
        <parameter name=\"charge\"              >    -1                    </parameter>\
        <parameter name=\"mass\"                >    1.0                   </parameter>\
        <attrib name=\"position\" datatype=\"posArray\" condition=\"0\">\
                 0.00000000        0.00000000        3.00000000\
                 3.00000000        0.00000000        3.00000000\
                 0.00000000        3.00000000        3.00000000\
                 3.00000000        3.00000000        3.00000000\
        </attrib>\
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
                 3.00000000        0.00000000        0.00000000\
                 0.00000000        3.00000000        0.00000000\
                 3.00000000        3.00000000        0.00000000\
                 0.00000000        0.00000000        3.00000000\
                 3.00000000        0.00000000        3.00000000\
                 0.00000000        3.00000000        3.00000000\
                 3.00000000        3.00000000        3.00000000\
        </attrib>\
     </group>\
  </particleset>\
</tmp> \
";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root  = doc.getRoot();
  xmlNodePtr part1 = xmlFirstElementChild(root);
  xmlNodePtr part2 = xmlNextElementSibling(part1);
  xmlNodePtr part3 = xmlNextElementSibling(part2);

  // read lattice
  ParticleSet::ParticleLayout_t* SimulationCell = new ParticleSet::ParticleLayout_t;
  LatticeParser lp(*SimulationCell);
  lp.put(part1);
  SimulationCell->print(app_log(), 0);

  // read particle set
  ParticleSet ions, electrons;
  Tensor<int, 3> tmat; // assuming OHMMSDIM==3
  tmat(0,0) = 1;
  tmat(1,1) = 1;
  tmat(2,2) = 1;

  XMLParticleParser parse_electrons(electrons, tmat);
  parse_electrons.put(part2);

  XMLParticleParser parse_ions(ions, tmat);
  parse_ions.put(part3);

  REQUIRE( electrons.getName() == "e" );
  REQUIRE( ions.getName() == "ion0"   );
  REQUIRE( ions.SameMass      );
  REQUIRE( electrons.SameMass );

  // calculate particle distances
  electrons.Lattice.copy(*SimulationCell);
  ions.Lattice.copy(*SimulationCell); // is this applied in qmcpack executable?
  // better be, electron-proton distances used in PairCorrelation estimator
  electrons.addTable(ions,DT_AOS);
  electrons.update();

  // get target particle set's distance table data
  int tid = electrons.getTable(ions);
  DistanceTableData* dtable = electrons.DistTables[tid];
  REQUIRE(dtable->getName() == "ion0_e");

  int num_src = ions.getTotalNum();
  int num_tar = electrons.getTotalNum();
  std::vector<double> expect;
  expect.resize(num_src*num_tar);
  int idx(0);
  for (int iat=0;iat<num_src;iat++){
    for (int jat=0;jat<num_tar;jat++,idx++){
      double dist = dtable->r(dtable->loc(iat,jat));
      expect[idx] = dist; 
    }
  }

  // move a particle by a lattice vector
  ParticleSet::SingleParticlePos_t disp(6.0,0,0);
  electrons.makeMove(0,disp); // temporary change written in distance table
  electrons.acceptMove(0);    // update distance table with temporary change
  
  // move some more
  disp[1] = -6.0;
  electrons.makeMove(1,disp);
  electrons.acceptMove(1);   

  // move some more more
  disp[2] = 12.0;
  electrons.makeMove(3,disp);
  electrons.acceptMove(3);   

  idx = 0;
  for (int iat=0;iat<num_src;iat++){
    for (int jat=0;jat<num_tar;jat++,idx++){
      double dist = dtable->r(dtable->loc(iat,jat));
      REQUIRE( expect[idx] == Approx(dist) );
    }
  }

} // TEST_CASE distance_pbc_z

} // namespace qmcplusplus
