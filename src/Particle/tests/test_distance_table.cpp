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

  const char* particles = "<tmp> \
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
  tmat(0, 0) = 1;
  tmat(1, 1) = 1;
  tmat(2, 2) = 1;

  // read particle set
  ParticleSet ions, electrons;

  XMLParticleParser parse_electrons(electrons, tmat);
  parse_electrons.put(part1);

  XMLParticleParser parse_ions(ions, tmat);
  parse_ions.put(part2);

  REQUIRE(electrons.getName() == "e");
  REQUIRE(ions.getName() == "ion0");
  REQUIRE(ions.SameMass);
  REQUIRE(electrons.SameMass);

  // calculate particle distances
#ifdef ENABLE_SOA
  const int tid = electrons.addTable(ions, DT_SOA);
#else
  const int tid = electrons.addTable(ions, DT_AOS);
#endif
  electrons.update();

  // get target particle set's distance table data
  const auto& dtable = electrons.getDistTable(tid);
  REQUIRE(dtable.getName() == "ion0_e");

  REQUIRE(dtable.sources() == ions.getTotalNum());
  REQUIRE(dtable.targets() == electrons.getTotalNum());

  double expect[] = {0.2, 0.2, 0.3, 0.7};
  int idx(0);
  for (int iat = 0; iat < dtable.sources(); iat++)
  {
    for (int jat = 0; jat < dtable.targets(); jat++, idx++)
    {
      // note: target particle set is special (electrons in this case)
      // int tid = target_pset.addTable(source_pset, DT_AOS);
      // const auto& dtable = target_pset.getDistTable(tid);
      // dtable.loc(source_ptcl_idx,target_ptcl_idx) !! source first target second !?
#ifdef ENABLE_SOA
      double dist = dtable.getDistRow(jat)[iat];
#else
      double dist = dtable.r(dtable.loc(iat, jat));
#endif
      REQUIRE(dist == Approx(expect[idx]));
    }
  }

#ifdef ENABLE_SOA
  TinyVector<double, 3> displ1 = dtable.getDisplacements()[0][0];
  REQUIRE(displ1[0] == Approx(0.0));
  REQUIRE(displ1[1] == Approx(0.0));
  REQUIRE(displ1[2] == Approx(-0.2));
#else
  TinyVector<double, 3> displ1 = dtable.displacement(0,0);
  REQUIRE(displ1[0] == Approx(0.0));
  REQUIRE(displ1[1] == Approx(0.0));
  REQUIRE(displ1[2] == Approx(0.2));
#endif

#ifdef ENABLE_SOA
  TinyVector<double, 3> displ2 = dtable.getDisplacements()[0][1];
  REQUIRE(displ2[0] == Approx(0.0));
  REQUIRE(displ2[1] == Approx(0.0));
  REQUIRE(displ2[2] == Approx(0.3));
#else
  TinyVector<double, 3> displ2 = dtable.displacement(1,0);
  REQUIRE(displ2[0] == Approx(0.0));
  REQUIRE(displ2[1] == Approx(0.0));
  REQUIRE(displ2[2] == Approx(-0.3));
#endif

  // get distance between target="e" group="u" iat=0 and source="ion0" group="H" jat=1

} // TEST_CASE distance_open_z

TEST_CASE("distance_open_xy", "[distance_table][xml]")
{
  OHMMS::Controller->initialize(0, NULL);

  const char* particles = "<tmp> \
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
  tmat(0, 0) = 1;
  tmat(1, 1) = 1;
  tmat(2, 2) = 1;

  // read particle set
  ParticleSet ions, electrons;

  XMLParticleParser parse_electrons(electrons, tmat);
  parse_electrons.put(part1);

  XMLParticleParser parse_ions(ions, tmat);
  parse_ions.put(part2);

  REQUIRE(electrons.getName() == "e");
  REQUIRE(ions.getName() == "ion0");
  REQUIRE(ions.SameMass);
  REQUIRE(electrons.SameMass);

  // calculate particle distances
#ifdef ENABLE_SOA
  const int tid = electrons.addTable(ions, DT_SOA);
#else
  const int tid = electrons.addTable(ions, DT_AOS);
#endif
  electrons.update();

  // get distance table attached to target particle set (electrons)
  const auto& dtable = electrons.getDistTable(tid);
  REQUIRE(dtable.getName() == "ion0_e");

  REQUIRE(dtable.sources() == ions.getTotalNum());
  REQUIRE(dtable.targets() == electrons.getTotalNum());

  // calculate distance, one source particle at a time i.e.
  // H0 - e0: 0.7
  // H0 - e1: 1.0
  // H0 - e2: 0.9
  // H1 - e0: 0.3
  // etc.
  double expect[] = {0.7, 1.0, 0.9, 0.3, std::sqrt(2), 1.9};
  int idx(0);
  for (int iat = 0; iat < dtable.sources(); iat++)
  {
    for (int jat = 0; jat < dtable.targets(); jat++, idx++)
    {
#ifdef ENABLE_SOA
      double dist = dtable.getDistRow(jat)[iat];
#else
      double dist = dtable.r(dtable.loc(iat, jat));
#endif
      REQUIRE(dist == Approx(expect[idx]));
    }
  }

} // TEST_CASE distance_open_xy

TEST_CASE("distance_open_species_deviation", "[distance_table][xml]")
{
  // pull out distances between specific species

  OHMMS::Controller->initialize(0, NULL);

  const char* particles = "<tmp> \
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
  tmat(0, 0) = 1;
  tmat(1, 1) = 1;
  tmat(2, 2) = 1;

  // read particle set
  ParticleSet ions, electrons;

  XMLParticleParser parse_electrons(electrons, tmat);
  parse_electrons.put(part1);

  XMLParticleParser parse_ions(ions, tmat);
  parse_ions.put(part2);

  REQUIRE(electrons.getName() == "e");
  REQUIRE(ions.getName() == "ion0");
  REQUIRE(ions.SameMass);
  REQUIRE(electrons.SameMass);

  // calculate particle distances
#ifdef ENABLE_SOA
  const int tid = electrons.addTable(ions, DT_SOA);
#else
  const int tid = electrons.addTable(ions, DT_AOS);
#endif
  electrons.update();

  // get distance table attached to target particle set (electrons)
  const auto& dtable = electrons.getDistTable(tid);
  REQUIRE(dtable.getName() == "ion0_e");

  // get the electron species set
  SpeciesSet& especies(electrons.getSpeciesSet());

  REQUIRE(dtable.sources() == ions.getTotalNum());
  REQUIRE(dtable.targets() == electrons.getTotalNum());

  // !! assume "u" and "H" groups have the same number of particles
  double latdev2 = 0.0; // mean-squared deviation from lattice
  int cur_jat(-1);      // keep an index to the last found target particle
  double expect[] = {0.7, std::sqrt(2)};
  int idx(0);
  for (int iat = 0; iat < dtable.sources(); iat++, idx++)
  {
    for (int jat = cur_jat + 1; jat < dtable.targets(); jat++, idx++)
    {
      // find next "u"
      int species_id      = electrons.GroupID[jat];
      string species_name = especies.speciesName[species_id];
      if (species_name != "u")
        continue;

      // calculate distance from lattice site iat
#ifdef ENABLE_SOA
      double dist = dtable.getDistRow(jat)[iat];
#else
      double dist = dtable.r(dtable.loc(iat, jat));
#endif
      latdev2 += std::pow(dist, 2); // !? pow(x,2) does what?
      REQUIRE(dist == Approx(expect[idx]));
      cur_jat = jat;
      break;
    }
  }
  REQUIRE(latdev2 / ions.getTotalNum() == Approx(1.245));

} // TEST_CASE distance_open_species_deviation

TEST_CASE("distance_pbc_z", "[distance_table][xml]")
{
  // test that particle distances are properly calculated under periodic boundary condition
  // There are many details in this example, but the main idea is simple: When a particle is moved by a full lattice vector, no distance should change.

  OHMMS::Controller->initialize(0, NULL);

  const char* particles = "<tmp> \
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
  ParticleSet::ParticleLayout_t SimulationCell;
  LatticeParser lp(SimulationCell);
  lp.put(part1);
  SimulationCell.print(app_log(), 0);

  // read particle set
  ParticleSet ions, electrons;
  Tensor<int, 3> tmat; // assuming OHMMSDIM==3
  tmat(0, 0) = 1;
  tmat(1, 1) = 1;
  tmat(2, 2) = 1;

  // enforce global Lattice on ions and electrons
  ions.Lattice = SimulationCell;
  electrons.Lattice = SimulationCell;

  XMLParticleParser parse_electrons(electrons, tmat);
  parse_electrons.put(part2);

  XMLParticleParser parse_ions(ions, tmat);
  parse_ions.put(part3);

  REQUIRE(electrons.getName() == "e");
  REQUIRE(ions.getName() == "ion0");
  REQUIRE(ions.SameMass);
  REQUIRE(electrons.SameMass);

  // calculate particle distances
#ifdef ENABLE_SOA
  const int ei_tid = electrons.addTable(ions, DT_SOA);
#else
  const int ei_tid = electrons.addTable(ions, DT_AOS);
#endif
  electrons.update();
  ions.update();

  // get target particle set's distance table data
  const auto& ei_dtable = electrons.getDistTable(ei_tid);
  REQUIRE(ei_dtable.getName() == "ion0_e");

  REQUIRE(ei_dtable.sources() == ions.getTotalNum());
  REQUIRE(ei_dtable.targets() == electrons.getTotalNum());

  int num_src = ions.getTotalNum();
  int num_tar = electrons.getTotalNum();
  std::vector<double> expect;
  expect.resize(num_src * num_tar);
  int idx(0);
  for (int iat = 0; iat < num_src; iat++)
  {
    for (int jat = 0; jat < num_tar; jat++, idx++)
    {
#ifdef ENABLE_SOA
      double dist = ei_dtable.getDistRow(jat)[iat];
#else
      double dist = ei_dtable.r(ei_dtable.loc(iat, jat));
#endif
      expect[idx] = dist;
    }
  }

  // move a particle by a lattice vector
  ParticleSet::SingleParticlePos_t disp(6.0, 0, 0);
  electrons.makeMove(0, disp); // temporary change written in distance table
  electrons.acceptMove(0);     // update distance table with temporary change

  // move some more
  disp[1] = -6.0;
  electrons.makeMove(1, disp);
  electrons.acceptMove(1);

  // move some more more
  disp[2] = 12.0;
  electrons.makeMove(3, disp);
  electrons.acceptMove(3);

  idx = 0;
  for (int iat = 0; iat < num_src; iat++)
  {
    for (int jat = 0; jat < num_tar; jat++, idx++)
    {
#ifdef ENABLE_SOA
      double dist = ei_dtable.getDistRow(jat)[iat];
#else
      double dist = ei_dtable.r(ei_dtable.loc(iat, jat));
#endif
      REQUIRE(expect[idx] == Approx(dist));
    }
  }

#ifdef ENABLE_SOA
  const int ee_tid = electrons.addTable(electrons, DT_SOA);
  // get target particle set's distance table data
  const auto& ee_dtable = electrons.getDistTable(ee_tid);
  REQUIRE(ee_dtable.getName() == "e_e");
  electrons.update();

  disp[0] = 0.2;
  disp[1] = 0.1;
  disp[2] = 0.3;

  electrons.makeMove(0, disp, false);
  REQUIRE(ee_dtable.getTempDists()[1] == Approx(2.8178005607));
  REQUIRE(ee_dtable.getTempDispls()[1][0] == Approx(2.8));
  REQUIRE(ee_dtable.getTempDispls()[1][1] == Approx(-0.1));
  REQUIRE(ee_dtable.getTempDispls()[1][2] == Approx(-0.3));
  REQUIRE(ee_dtable.getDistRow(1)[0] == Approx(3.0));
  REQUIRE(ee_dtable.getDisplRow(1)[0][0] == Approx(3.0));
  REQUIRE(ee_dtable.getDisplRow(1)[0][1] == Approx(0.0));
  REQUIRE(ee_dtable.getDisplRow(1)[0][2] == Approx(0.0));
  electrons.rejectMove(0);

  electrons.makeMove(0, disp);
  REQUIRE(ee_dtable.getTempDists()[1] == Approx(2.8178005607));
  REQUIRE(ee_dtable.getTempDispls()[1][0] == Approx(2.8));
  REQUIRE(ee_dtable.getTempDispls()[1][1] == Approx(-0.1));
  REQUIRE(ee_dtable.getTempDispls()[1][2] == Approx(-0.3));
  REQUIRE(ee_dtable.getOldDists()[1] == Approx(3.0));
  REQUIRE(ee_dtable.getOldDispls()[1][0] == Approx(-3.0));
  REQUIRE(ee_dtable.getOldDispls()[1][1] == Approx(0.0));
  REQUIRE(ee_dtable.getOldDispls()[1][2] == Approx(0.0));
  REQUIRE(ee_dtable.getDistRow(1)[0] == Approx(3.0));
  REQUIRE(ee_dtable.getDisplRow(1)[0][0] == Approx(3.0));
  REQUIRE(ee_dtable.getDisplRow(1)[0][1] == Approx(0.0));
  REQUIRE(ee_dtable.getDisplRow(1)[0][2] == Approx(0.0));
  electrons.acceptMove(0);

  REQUIRE(ee_dtable.getDistRow(1)[0] == Approx(2.8178005607));
  REQUIRE(ee_dtable.getDisplRow(1)[0][0] == Approx(-2.8));
  REQUIRE(ee_dtable.getDisplRow(1)[0][1] == Approx(0.1));
  REQUIRE(ee_dtable.getDisplRow(1)[0][2] == Approx(0.3));
#endif
} // TEST_CASE distance_pbc_z

} // namespace qmcplusplus
