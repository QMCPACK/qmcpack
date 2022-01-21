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
#include "Lattice/CrystalLattice.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/PairCorrEstimator.h"
#include "Particle/ParticleSetPool.h"

#include <stdio.h>
#include <string>

using std::string;


/*
  Minimal input block (based off HEG problem).
  This is a dumb way to construct a test, but I can't figure out
  the right way to initialize a ParticleSet, so instead hijack 
  the ParticleSetPool class to make a ParticleSet for me.

  NB: Hardcoded 8 particles (4 each of u/d) on a B1 (NaCl) lattice.
*/

// Lattice block
const char* lat_xml = "<simulationcell> \
                        <parameter name=\"lattice\" units=\"bohr\"> \
                          2.0 0.0 0.0 \
                          0.0 2.0 0.0 \
                          0.0 0.0 2.0 \
                        </parameter> \
                        <parameter name=\"bconds\"> \
                           p p p \
                        </parameter> \
                        <parameter name=\"LR_dim_cutoff\" >     6  </parameter> \
                       </simulationcell>";

// Particleset block
const char* pset_xml = "<particleset name=\"e\" random=\"yes\"> \
                          <group name=\"u\" size=\"4\" mass=\"1.0\"> \
                            <parameter name=\"charge\" >   -1  </parameter> \
                            <parameter name=\"mass\"   >  1.0  </parameter> \
                          </group> \
                          <group name=\"d\" size=\"4\" mass=\"1.0\"> \
                            <parameter name=\"charge\" >   -1  </parameter> \
                            <parameter name=\"mass\"   >  1.0  </parameter> \
                          </group> \
                        </particleset>";

// PairCorrEstimator block
const char* gofr_xml = "<estimator type=\"gofr\" name=\"gofr\" rmax=\"2.0\" num_bin=\"99\" />";


namespace qmcplusplus
{
TEST_CASE("Pair Correlation", "[hamiltonian]")
{
  std::cout << std::fixed;
  std::cout << std::setprecision(8);
  typedef QMCTraits::RealType RealType;

  Communicate* c = OHMMS::Controller;

  // XML parser
  Libxml2Document doc;

  std::cout << "\n\n\ntest_paircorr:: START\n";

  // TEST new idea: ParticlesetPool to make a ParticleSet
  bool lat_okay = doc.parseFromString(lat_xml);
  REQUIRE(lat_okay);
  xmlNodePtr lat_xml_root = doc.getRoot();

  ParticleSetPool pset_builder(c, "pset_builder");
  pset_builder.readSimulationCellXML(lat_xml_root); // Builds lattice

  bool pset_okay = doc.parseFromString(pset_xml);
  REQUIRE(pset_okay);
  xmlNodePtr pset_xml_root = doc.getRoot();
  pset_builder.put(pset_xml_root); // Builds ParticleSet

  // Get the (now assembled) ParticleSet, do simple sanity checks, then print info
  ParticleSet* elec = pset_builder.getParticleSet("e");

  std::cout << "cheeeee " << elec->getLattice().R << std::endl;
  REQUIRE(elec->isSameMass());
  REQUIRE(elec->getName() == "e");

  // Move the particles manually onto B1 lattice
  // NB: Spins are grouped contiguously
  // Up spins
  elec->R[0][0] = 0.0;
  elec->R[0][1] = 0.0;
  elec->R[0][2] = 0.0;
  elec->R[1][0] = 1.0;
  elec->R[1][1] = 1.0;
  elec->R[1][2] = 0.0;
  elec->R[2][0] = 1.0;
  elec->R[2][1] = 0.0;
  elec->R[2][2] = 1.0;
  elec->R[3][0] = 0.0;
  elec->R[3][1] = 1.0;
  elec->R[3][2] = 1.0;

  // Down spins
  elec->R[4][0] = 1.0;
  elec->R[4][1] = 0.0;
  elec->R[4][2] = 0.0;
  elec->R[5][0] = 0.0;
  elec->R[5][1] = 1.0;
  elec->R[5][2] = 0.0;
  elec->R[6][0] = 0.0;
  elec->R[6][1] = 0.0;
  elec->R[6][2] = 1.0;
  elec->R[7][0] = 1.0;
  elec->R[7][1] = 1.0;
  elec->R[7][2] = 1.0;

  elec->get(std::cout); // print particleset info to stdout

  // Set up the distance table, match expected layout
  const int ee_table_id = elec->addTable(*elec);

  const auto& dii(elec->getDistTable(ee_table_id));
  elec->update(); // distance table evaluation here

  // Make a PairCorrEstimator, call put() to set up internals
  std::string name = elec->getName();
  PairCorrEstimator paircorr(*elec, name);
  bool gofr_okay = doc.parseFromString(gofr_xml);
  REQUIRE(gofr_okay);
  xmlNodePtr gofr_xml_root = doc.getRoot();
  paircorr.put(gofr_xml_root);
  paircorr.addObservables(elec->PropertyList, elec->Collectables);

  // Compute g(r) then print it to stdout
  // NB: Hardcoded to match hardcoded xml above. ***Fragile!!!***
  paircorr.evaluate(*elec);

  // Hardcoded gofr parameters - MUST match xml input for PairCorrEstimator!!!
  // This might cause a segfault if it does not match xml input!
  const RealType Rmax   = 2.0;
  const int Nbins       = 99;
  const RealType deltaR = Rmax / static_cast<RealType>(Nbins);

  auto gofr = elec->Collectables;
  std::cout << "\n";
  std::cout << "gofr:\n";
  std::cout << std::fixed;
  std::cout << std::setprecision(6);
  std::cout << std::setw(4) << "i"
            << "  " << std::setw(12) << "r"
            << "  " << std::setw(12) << "uu"
            << "  " << std::setw(12) << "ud"
            << "  " << std::setw(12) << "dd"
            << "\n";
  std::cout << "============================================================\n";

  for (int i = 0; i < Nbins; i++)
  {
    std::cout << std::setw(4) << i << "  " << std::setw(12) << i * deltaR << "  " << std::setw(12) << gofr[i] << "  "
              << std::setw(12) << gofr[i + Nbins] << "  " << std::setw(12) << gofr[i + 2 * Nbins] << "\n";
  }

  // Verify against analytic result.
  // NB: At the moment the tolerance is fairly loose because there
  //     is about 0.08-0.01 disagreement between QMCPACK and analytic
  //     result, depending on the bin. I'm not sure about the source
  //     of the disagreement because norm_factor is now exact.
  //     The only thing left is the distance table (I think)...
  const RealType eps = 1E-02; // tolerance

  // Nearest neighbor peak (ud) | Distance = 1
  const int bin_nn = 49;
  REQUIRE(std::fabs(gofr[49] - 0.00000000) < eps);
  REQUIRE(std::fabs(gofr[148] - 23.6361163) < eps);
  REQUIRE(std::fabs(gofr[247] - 0.00000000) < eps);

  // 2nd-nearest neighbor peak (uu/dd) | Distance = sqrt(2)
  const int bin_2n = 70;
  REQUIRE(std::fabs(gofr[70] - 15.536547) < eps);
  REQUIRE(std::fabs(gofr[169] - 0.0000000) < eps);
  REQUIRE(std::fabs(gofr[268] - 15.536547) < eps);

  // 3rd-nearest neighbor peak (ud) | Distance = sqrt(3)
  REQUIRE(std::fabs(gofr[85] - 0.00000000) < eps);
  REQUIRE(std::fabs(gofr[184] - 2.6408410) < eps);
  REQUIRE(std::fabs(gofr[283] - 0.0000000) < eps);

  std::cout << "test_paircorr:: STOP\n";
}

TEST_CASE("Pair Correlation Pair Index", "[hamiltonian]")
{
  // Check generation of pair id for 2 species groups
  REQUIRE(PairCorrEstimator::gen_pair_id(0, 0, 2) == 0);
  REQUIRE(PairCorrEstimator::gen_pair_id(0, 1, 2) == 1);
  REQUIRE(PairCorrEstimator::gen_pair_id(1, 0, 2) == 1);
  REQUIRE(PairCorrEstimator::gen_pair_id(1, 1, 2) == 2);

  // Check generation of pair id for 3 species groups
  REQUIRE(PairCorrEstimator::gen_pair_id(0, 0, 3) == 0);
  REQUIRE(PairCorrEstimator::gen_pair_id(0, 1, 3) == 1);
  REQUIRE(PairCorrEstimator::gen_pair_id(0, 2, 3) == 2);
  REQUIRE(PairCorrEstimator::gen_pair_id(1, 0, 3) == 1);
  REQUIRE(PairCorrEstimator::gen_pair_id(1, 1, 3) == 3);
  REQUIRE(PairCorrEstimator::gen_pair_id(1, 2, 3) == 4);
  REQUIRE(PairCorrEstimator::gen_pair_id(2, 0, 3) == 2);
  REQUIRE(PairCorrEstimator::gen_pair_id(2, 1, 3) == 4);
  REQUIRE(PairCorrEstimator::gen_pair_id(2, 2, 3) == 5);
}
} // namespace qmcplusplus
