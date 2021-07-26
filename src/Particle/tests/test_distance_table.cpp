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
  const int tid = electrons.addTable(ions);
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
      double dist = dtable.getDistRow(jat)[iat];
      CHECK(dist == Approx(expect[idx]));
    }
  }

  TinyVector<double, 3> displ1 = dtable.getDisplacements()[0][0];
  CHECK(displ1[0] == Approx(0.0));
  CHECK(displ1[1] == Approx(0.0));
  CHECK(displ1[2] == Approx(-0.2));

  TinyVector<double, 3> displ2 = dtable.getDisplacements()[0][1];
  CHECK(displ2[0] == Approx(0.0));
  CHECK(displ2[1] == Approx(0.0));
  CHECK(displ2[2] == Approx(0.3));

  // get distance between target="e" group="u" iat=0 and source="ion0" group="H" jat=1

} // TEST_CASE distance_open_z

TEST_CASE("distance_open_xy", "[distance_table][xml]")
{
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
  const int tid = electrons.addTable(ions);
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
      double dist = dtable.getDistRow(jat)[iat];
      CHECK(dist == Approx(expect[idx]));
    }
  }

} // TEST_CASE distance_open_xy

TEST_CASE("distance_open_species_deviation", "[distance_table][xml]")
{
  // pull out distances between specific species


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
  const int tid = electrons.addTable(ions);
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
      double dist = dtable.getDistRow(jat)[iat];
      latdev2 += std::pow(dist, 2); // !? pow(x,2) does what?
      CHECK(dist == Approx(expect[idx]));
      cur_jat = jat;
      break;
    }
  }
  REQUIRE(latdev2 / ions.getTotalNum() == Approx(1.245));

} // TEST_CASE distance_open_species_deviation

void parse_electron_ion_pbc_z(ParticleSet& ions, ParticleSet& electrons)
{
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
  Tensor<int, 3> tmat; // assuming OHMMSDIM==3
  tmat(0, 0) = 1;
  tmat(1, 1) = 1;
  tmat(2, 2) = 1;

  // enforce global Lattice on ions and electrons
  ions.Lattice      = SimulationCell;
  electrons.Lattice = SimulationCell;

  XMLParticleParser parse_electrons(electrons, tmat);
  parse_electrons.put(part2);

  XMLParticleParser parse_ions(ions, tmat);
  parse_ions.put(part3);

  REQUIRE(electrons.getName() == "e");
  REQUIRE(ions.getName() == "ion0");
  REQUIRE(ions.SameMass);
  REQUIRE(electrons.SameMass);

}

TEST_CASE("distance_pbc_z", "[distance_table][xml]")
{
  // test that particle distances are properly calculated under periodic boundary condition
  // There are many details in this example, but the main idea is simple: When a particle is moved by a full lattice vector, no distance should change.

  ParticleSet ions, electrons;

  parse_electron_ion_pbc_z(ions, electrons);

  // calculate particle distances
  const int ei_tid = electrons.addTable(ions);
  electrons.update();
  ions.update();

  // get target particle set's distance table data
  const auto& ei_dtable = electrons.getDistTable(ei_tid);
  CHECK(ei_dtable.getName() == "ion0_e");

  CHECK(ei_dtable.sources() == ions.getTotalNum());
  CHECK(ei_dtable.targets() == electrons.getTotalNum());

  int num_src = ions.getTotalNum();
  int num_tar = electrons.getTotalNum();
  std::vector<double> expect;
  expect.resize(num_src * num_tar);
  int idx(0);
  for (int iat = 0; iat < num_src; iat++)
  {
    for (int jat = 0; jat < num_tar; jat++, idx++)
    {
      double dist = ei_dtable.getDistRow(jat)[iat];
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
      double dist = ei_dtable.getDistRow(jat)[iat];
      CHECK(expect[idx] == Approx(dist));
    }
  }

  const int ee_tid = electrons.addTable(electrons);
  // get target particle set's distance table data
  const auto& ee_dtable = electrons.getDistTable(ee_tid);
  CHECK(ee_dtable.getName() == "e_e");
  electrons.update();

  // shift electron 0 a bit to avoid box edges.
  ParticleSet::SingleParticlePos_t shift(0.1, 0.2, -0.1);
  electrons.makeMove(0, shift);
  electrons.acceptMove(0);

  disp[0] = 0.2;
  disp[1] = 0.1;
  disp[2] = 0.3;

  electrons.makeMove(0, disp, false);
  CHECK(ee_dtable.getTempDists()[1] == Approx(2.7239676944));
  CHECK(ee_dtable.getTempDispls()[1][0] == Approx(2.7));
  CHECK(ee_dtable.getTempDispls()[1][1] == Approx(-0.3));
  CHECK(ee_dtable.getTempDispls()[1][2] == Approx(-0.2));
  CHECK(ee_dtable.getDistRow(1)[0] == Approx(2.908607914));
  CHECK(ee_dtable.getDisplRow(1)[0][0] == Approx(-2.9));
  CHECK(ee_dtable.getDisplRow(1)[0][1] == Approx(0.2));
  CHECK(ee_dtable.getDisplRow(1)[0][2] == Approx(-0.1));
  electrons.rejectMove(0);

  electrons.makeMove(0, disp);
  CHECK(ee_dtable.getTempDists()[1] == Approx(2.7239676944));
  CHECK(ee_dtable.getTempDispls()[1][0] == Approx(2.7));
  CHECK(ee_dtable.getTempDispls()[1][1] == Approx(-0.3));
  CHECK(ee_dtable.getTempDispls()[1][2] == Approx(-0.2));
  CHECK(ee_dtable.getOldDists()[1] == Approx(2.908607914));
  CHECK(ee_dtable.getOldDispls()[1][0] == Approx(2.9));
  CHECK(ee_dtable.getOldDispls()[1][1] == Approx(-0.2));
  CHECK(ee_dtable.getOldDispls()[1][2] == Approx(0.1));
  CHECK(ee_dtable.getDistRow(1)[0] == Approx(2.908607914));
  CHECK(ee_dtable.getDisplRow(1)[0][0] == Approx(-2.9));
  CHECK(ee_dtable.getDisplRow(1)[0][1] == Approx(0.2));
  CHECK(ee_dtable.getDisplRow(1)[0][2] == Approx(-0.1));
  electrons.acceptMove(0);

  CHECK(ee_dtable.getDistRow(1)[0] == Approx(2.7239676944));
  CHECK(ee_dtable.getDisplRow(1)[0][0] == Approx(-2.7));
  CHECK(ee_dtable.getDisplRow(1)[0][1] == Approx(0.3));
  CHECK(ee_dtable.getDisplRow(1)[0][2] == Approx(0.2));

  // revert previous move
  electrons.makeMove(0, -disp);
  electrons.acceptMove(0);

  // test forward mode
  electrons.makeMove(0, disp);
  electrons.accept_rejectMove(0, true, true);

  CHECK(ee_dtable.getDistRow(1)[0] == Approx(2.908607914));
  CHECK(ee_dtable.getDisplRow(1)[0][0] == Approx(-2.9));
  CHECK(ee_dtable.getDisplRow(1)[0][1] == Approx(0.2));
  CHECK(ee_dtable.getDisplRow(1)[0][2] == Approx(-0.1));

  electrons.makeMove(1, disp);
  electrons.accept_rejectMove(1, false, true);
  CHECK(ee_dtable.getDistRow(1)[0] == Approx(2.7239676944));
  CHECK(ee_dtable.getDisplRow(1)[0][0] == Approx(-2.7));
  CHECK(ee_dtable.getDisplRow(1)[0][1] == Approx(0.3));
  CHECK(ee_dtable.getDisplRow(1)[0][2] == Approx(0.2));
} // TEST_CASE distance_pbc_z

void test_distance_pbc_z_batched_APIs(DynamicCoordinateKind test_kind)
{
  // test that particle distances are properly calculated under periodic boundary condition
  // There are many details in this example, but the main idea is simple: When a particle is moved by a full lattice vector, no distance should change.

  ParticleSet ions, electrons(test_kind);
  parse_electron_ion_pbc_z(ions, electrons);

  // calculate particle distances
  ions.update();
  const int ee_tid = electrons.addTable(electrons);
  // get target particle set's distance table data
  const auto& ee_dtable = electrons.getDistTable(ee_tid);
  CHECK(ee_dtable.getName() == "e_e");
  electrons.update();

  // shift electron 0 a bit to avoid box edges.
  ParticleSet::SingleParticlePos_t shift(0.1, 0.2, -0.1);
  electrons.makeMove(0, shift);
  electrons.accept_rejectMove(0, true, false);
  electrons.donePbyP();

  ParticleSet electrons_clone(electrons);
  RefVectorWithLeader<ParticleSet> p_list(electrons);
  p_list.push_back(electrons);
  p_list.push_back(electrons_clone);

  std::vector<ParticleSet::SingleParticlePos_t> disp{{0.2, 0.1, 0.3}, {0.2, 0.1, 0.3}};

  ParticleSet::mw_makeMove(p_list, 0, disp);
  ParticleSet::mw_accept_rejectMove(p_list, 0, {true, true}, true);
  ParticleSet::mw_makeMove(p_list, 1, disp);
  ParticleSet::mw_accept_rejectMove(p_list, 1, {false, false}, true);

  ParticleSet::mw_donePbyP(p_list);
  CHECK(ee_dtable.getDistRow(1)[0] == Approx(2.7239676944));
  CHECK(ee_dtable.getDisplRow(1)[0][0] == Approx(-2.7));
  CHECK(ee_dtable.getDisplRow(1)[0][1] == Approx(0.3));
  CHECK(ee_dtable.getDisplRow(1)[0][2] == Approx(0.2));
} // test_distance_pbc_z_batched_APIs

TEST_CASE("distance_pbc_z batched APIs", "[distance_table][xml]")
{
  test_distance_pbc_z_batched_APIs(DynamicCoordinateKind::DC_POS);
#if defined(ENABLE_OFFLOAD)
  test_distance_pbc_z_batched_APIs(DynamicCoordinateKind::DC_POS_OFFLOAD);
#endif
}

void test_distance_pbc_z_batched_APIs_ee_NEED_TEMP_DATA_ON_HOST(DynamicCoordinateKind test_kind)
{
  // test that particle distances are properly calculated under periodic boundary condition
  // There are many details in this example, but the main idea is simple: When a particle is moved by a full lattice vector, no distance should change.

  ParticleSet ions, electrons(test_kind);
  parse_electron_ion_pbc_z(ions, electrons);

  // calculate particle distances
  ions.update();
  const int ee_tid = electrons.addTable(electrons, DTModes::NEED_TEMP_DATA_ON_HOST);
  // get target particle set's distance table data
  const auto& ee_dtable = electrons.getDistTable(ee_tid);
  CHECK(ee_dtable.getName() == "e_e");
  electrons.update();

  // shift electron 0 a bit to avoid box edges.
  ParticleSet::SingleParticlePos_t shift(0.1, 0.2, -0.1);
  electrons.makeMove(0, shift);
  electrons.accept_rejectMove(0, true, false);
  electrons.donePbyP();

  ParticleSet electrons_clone(electrons);
  RefVectorWithLeader<ParticleSet> p_list(electrons);
  p_list.push_back(electrons);
  p_list.push_back(electrons_clone);

  std::vector<ParticleSet::SingleParticlePos_t> disp{{0.2, 0.1, 0.3}, {0.2, 0.1, 0.3}};

  ParticleSet::mw_makeMove(p_list, 0, disp);
  CHECK(ee_dtable.getTempDists()[1] == Approx(2.7239676944));
  CHECK(ee_dtable.getTempDispls()[1][0] == Approx(2.7));
  CHECK(ee_dtable.getTempDispls()[1][1] == Approx(-0.3));
  CHECK(ee_dtable.getTempDispls()[1][2] == Approx(-0.2));
  CHECK(ee_dtable.getOldDists()[1] == Approx(2.908607914));
  CHECK(ee_dtable.getOldDispls()[1][0] == Approx(2.9));
  CHECK(ee_dtable.getOldDispls()[1][1] == Approx(-0.2));
  CHECK(ee_dtable.getOldDispls()[1][2] == Approx(0.1));
  CHECK(ee_dtable.getDistRow(1)[0] == Approx(2.908607914));
  CHECK(ee_dtable.getDisplRow(1)[0][0] == Approx(-2.9));
  CHECK(ee_dtable.getDisplRow(1)[0][1] == Approx(0.2));
  CHECK(ee_dtable.getDisplRow(1)[0][2] == Approx(-0.1));
  ParticleSet::mw_accept_rejectMove(p_list, 0, {true, true}, true);

  CHECK(ee_dtable.getDistRow(1)[0] == Approx(2.908607914));
  CHECK(ee_dtable.getDisplRow(1)[0][0] == Approx(-2.9));
  CHECK(ee_dtable.getDisplRow(1)[0][1] == Approx(0.2));
  CHECK(ee_dtable.getDisplRow(1)[0][2] == Approx(-0.1));

  ParticleSet::mw_makeMove(p_list, 1, disp);
  ParticleSet::mw_accept_rejectMove(p_list, 1, {false, false}, true);
  CHECK(ee_dtable.getDistRow(1)[0] == Approx(2.7239676944));
  CHECK(ee_dtable.getDisplRow(1)[0][0] == Approx(-2.7));
  CHECK(ee_dtable.getDisplRow(1)[0][1] == Approx(0.3));
  CHECK(ee_dtable.getDisplRow(1)[0][2] == Approx(0.2));
} // test_distance_pbc_z_batched_APIs_ee_NEED_TEMP_DATA_ON_HOST

TEST_CASE("distance_pbc_z batched APIs ee NEED_TEMP_DATA_ON_HOST", "[distance_table][xml]")
{
  test_distance_pbc_z_batched_APIs_ee_NEED_TEMP_DATA_ON_HOST(DynamicCoordinateKind::DC_POS);
#if defined(ENABLE_OFFLOAD)
  test_distance_pbc_z_batched_APIs_ee_NEED_TEMP_DATA_ON_HOST(DynamicCoordinateKind::DC_POS_OFFLOAD);
#endif
}
} // namespace qmcplusplus
