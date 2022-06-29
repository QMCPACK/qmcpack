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
#include "LongRange/StructFact.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/SkAllEstimator.h"
#include "Particle/ParticleSetPool.h"
#include <stdio.h>
#include <string>
using std::string;


/*
  -- jpt 25/03/2020 --
  Deterministic unit test for the SkAllEstimator class.
  Despite the plumbing, it seems that _for the moment_ (March 2020)
  that only the e-e structure factor is computed and stored.

  This is the important piece, as we use it for finite size corrections.
  However, it looks like one might be able to generalize it for e-i and i-i
  at a later date. If that happens, this test should be updated!
 */


namespace qmcplusplus
{
TEST_CASE("SkAll", "[hamiltonian]")
{
  // Boiler plate setup
  std::cout << std::fixed;
  std::cout << std::setprecision(8);
  using RealType = QMCTraits::RealType;

  Communicate* c;
  c = OHMMS::Controller;

  // XML parser
  Libxml2Document doc;


  /*
    Minimal xml input blocks
  */
  // Lattice block
  const char* lat_xml = "<simulationcell>                                         \
                        <parameter name=\"lattice\" units=\"bohr\"> \
                          2.0 0.0 0.0 \
                          0.0 2.0 0.0 \
                          0.0 0.0 2.0 \
                        </parameter> \
                         <parameter name=\"bconds\">                            \
                           p p p                                                \
                         </parameter>                                           \
                         <parameter name=\"LR_dim_cutoff\" >     6  </parameter>\
                         <parameter name=\"rs\"            >   1.0  </parameter>\
                         <parameter name=\"nparticles\"    >     8  </parameter>\
                       </simulationcell>";

  // Particleset block for the electrons
  const char* elec_pset_xml = "<particleset name=\"e\" random=\"yes\">            \
                               <group name=\"u\" size=\"4\">	                \
                                 <parameter name=\"charge\" >   -1  </parameter>\
                                 <parameter name=\"mass\"   >  1.0  </parameter>\
                               </group>                                         \
                               <group name=\"d\" size=\"4\">                    \
                                 <parameter name=\"charge\" >   -1  </parameter>\
                                 <parameter name=\"mass\"   >  1.0  </parameter>\
                               </group>                                         \
                             </particleset>";

  // Particleset block for the ions
  const char* ion_pset_xml = "<particleset name=\"i\" size=\"4\">                 \
                              <group name=\"He\">                               \
                                <parameter name=\"charge\"      > 2 </parameter>\
                                <parameter name=\"valence\"     > 2 </parameter>\
                                <parameter name=\"atomicnumber\"> 2 </parameter>\
                              </group>                                          \
                              <attrib name=\"position\" datatype=\"posArray\"   \
                                                        condition=\"1\">        \
                                0.00 0.00 0.00                                  \
                                0.25 0.25 0.25                                  \
                                0.50 0.50 0.50                                  \
                                0.75 0.75 0.75                                  \
                              </attrib>                                         \
                              <attrib name=\"ionid\" datatype=\"stringArray\">  \
                                He He He He                                     \
                              </attrib>                                         \
                            </particleset>";

  // SKAllEstimator block, seems that writeionion does nothing?
  const char* skall_xml = "<estimator                                           \
                             name=\"Skall\" type=\"skall\"                      \
                             source=\"i\" target=\"e\" hdf5=\"yes\"             \
                           />";

  // Read in xml, add to pset_builder below
  bool lat_okay = doc.parseFromString(lat_xml);
  REQUIRE(lat_okay);
  xmlNodePtr lat_xml_root = doc.getRoot();

  std::cout << "\n\n\ntest_SkAllEstimator: START\n";

  // Build a ParticleSetPool - makes ParticleSets
  ParticleSetPool pset_builder(c, "pset_builder");

  // First attach the Lattice defined above
  pset_builder.readSimulationCellXML(lat_xml_root);

  // Now build the elec ParticleSet
  bool elec_pset_okay = doc.parseFromString(elec_pset_xml);
  REQUIRE(elec_pset_okay);
  xmlNodePtr elec_pset_xml_root = doc.getRoot();
  pset_builder.put(elec_pset_xml_root);

  // Get the (now assembled) elec ParticleSet, sanity check, report
  ParticleSet* elec = pset_builder.getParticleSet("e");
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


  // Get the (now assembled) ion ParticleSet, sanity check, report
  bool ion_pset_okay = doc.parseFromString(ion_pset_xml);
  REQUIRE(ion_pset_okay);
  xmlNodePtr ion_pset_xml_root = doc.getRoot();
  pset_builder.put(ion_pset_xml_root);

  // Finally, make the ion ParticleSet, sanity check, report
  // It seems that we need this only to construct skall, but it
  // is never used to evaluate the estimator in this test.
  ParticleSet* ion = pset_builder.getParticleSet("i");
  REQUIRE(ion->isSameMass());
  REQUIRE(ion->getName() == "i");
  ion->get(std::cout); // print particleset info to stdout


  // Set up the distance table, match expected layout
  const int ee_table_id = elec->addTable(*elec);

  const auto& dii(elec->getDistTable(ee_table_id));
  elec->update(); // distance table evaluation here

  // Check that the skall xml block is valid
  bool skall_okay = doc.parseFromString(skall_xml);
  REQUIRE(skall_okay);
  xmlNodePtr skall_xml_root = doc.getRoot();

  // Make a SkAllEstimator, call put() to set up internals
  SkAllEstimator skall(*ion, *elec);
  skall.put(skall_xml_root);
  skall.addObservables(elec->PropertyList, elec->Collectables);
  skall.get(app_log()); // pretty print settings

  // Hack to make a walker so that t_walker_ points to something
  // Only used to set t_walker_->Weight = 1 so that skall->evaluate()
  // doesn't segfault.
  // NB: setHistories(dummy) attaches dummy to t_walker_
  ParticleSet::Walker_t dummy = ParticleSet::Walker_t(1);
  skall.setHistories(dummy);
  skall.evaluate(*elec);

  // s(k) is computed by & held by the ParticleSet. SkAll just
  // takes that pre-computed s(k) and pulls out rho(k)..
  // In order to compare to analytic result, need the list
  // of k-vectors in cartesian coordinates.
  // Luckily, ParticleSet stores that in SK->getKLists().kpts_cart
  int nkpts = elec->getSimulationCell().getKLists().numk;
  std::cout << "\n";
  std::cout << "SkAll results:\n";
  std::cout << std::fixed;
  std::cout << std::setprecision(6);
  std::cout << std::setw(4) << "i"
            << "  " << std::setw(8) << "kx"
            << "  " << std::setw(8) << "ky"
            << "  " << std::setw(8) << "kz"
            << "  " << std::setw(8) << "rhok_r"
            << "  " << std::setw(8) << "rhok_i"
            << "  " << std::setw(8) << "c.c."
            << "\n";
  std::cout << "================================================================\n";

  // Extract rhok out of Collectables, print values
  auto rhok = elec->Collectables;
  std::cout << std::fixed;
  std::cout << std::setprecision(5);
  for (int k = 0; k < nkpts; k++)
  {
    auto kvec      = elec->getSimulationCell().getKLists().kpts_cart[k];
    RealType kx    = kvec[0];
    RealType ky    = kvec[1];
    RealType kz    = kvec[2];
    RealType rk_r  = rhok[nkpts + k];
    RealType rk_i  = rhok[2 * nkpts + k];
    RealType rk_cc = rhok[k];
    std::cout << std::setw(4) << k << "  " << std::setw(8) << kx << "  " << std::setw(8) << ky << "  " << std::setw(8)
              << kz << "  " << std::setw(8) << rk_r << "  " << std::setw(8) << rk_i << "  " << std::setw(8) << rk_cc
              << "\n";
  }

  /*    
    Verify against analytic result:
    NB: MUST match xml input!!!
    rho(-1.94889,  0.00000,  0.00000)=  2.52341 + -3.71748i
    rho( 0.00000, -1.94889,  0.00000)=  2.52341 + -3.71748i
    rho( 0.00000,  0.00000, -1.94889)=  2.52341 + -3.71748i
    rho( 0.00000,  0.00000,  1.94889)=  2.52341 +  3.71748i
    rho( 0.00000,  1.94889,  0.00000)=  2.52341 +  3.71748i
    rho( 1.94889,  0.00000,  0.00000)=  2.52341 +  3.71748i
    rho(-1.94889, -1.94889,  0.00000)= -0.93151 + -2.34518i
    rho(-1.94889,  0.00000, -1.94889)= -0.93151 + -2.34518i
    rho(-1.94889,  0.00000,  1.94889)=  2.52341 +  0.00000i
    rho(-1.94889,  1.94889,  0.00000)=  2.52341 +  0.00000i
    rho( 0.00000, -1.94889, -1.94889)= -0.93151 + -2.34518i
    rho( 0.00000, -1.94889,  1.94889)=  2.52341 +  0.00000i
    rho( 0.00000,  1.94889, -1.94889)=  2.52341 +  0.00000i
    rho( 0.00000,  1.94889,  1.94889)= -0.93151 +  2.34518i
    rho( 1.94889, -1.94889,  0.00000)=  2.52341 +  0.00000i
    rho( 1.94889,  0.00000, -1.94889)=  2.52341 +  0.00000i
    rho( 1.94889,  0.00000,  1.94889)= -0.93151 +  2.34518i
    rho( 1.94889,  1.94889,  0.00000)= -0.93151 +  2.34518i
    rho(-1.94889, -1.94889, -1.94889)= -1.38359 + -0.30687i
    rho(-1.94889, -1.94889,  1.94889)=  0.79595 + -1.17259i
    rho(-1.94889,  1.94889, -1.94889)=  0.79595 + -1.17259i
    rho(-1.94889,  1.94889,  1.94889)=  0.79595 +  1.17259i
    rho( 1.94889, -1.94889, -1.94889)=  0.79595 + -1.17259i
    rho( 1.94889, -1.94889,  1.94889)=  0.79595 +  1.17259i
    rho( 1.94889,  1.94889, -1.94889)=  0.79595 +  1.17259i
    rho( 1.94889,  1.94889,  1.94889)= -1.38359 +  0.30687i
  */

  const RealType eps = 1E-04; // tolerance
  REQUIRE(std::fabs(rhok[26] - 2.52341) < eps);
  REQUIRE(std::fabs(rhok[52] + 3.71748) < eps);
  REQUIRE(std::fabs(rhok[26 + 1] - 2.52341) < eps);
  REQUIRE(std::fabs(rhok[52 + 1] + 3.71748) < eps);
  REQUIRE(std::fabs(rhok[26 + 2] - 2.52341) < eps);
  REQUIRE(std::fabs(rhok[52 + 2] + 3.71748) < eps);

  REQUIRE(std::fabs(rhok[26 + 3] - 2.52341) < eps);
  REQUIRE(std::fabs(rhok[52 + 3] - 3.71748) < eps);
  REQUIRE(std::fabs(rhok[26 + 4] - 2.52341) < eps);
  REQUIRE(std::fabs(rhok[52 + 4] - 3.71748) < eps);
  REQUIRE(std::fabs(rhok[26 + 5] - 2.52341) < eps);
  REQUIRE(std::fabs(rhok[52 + 5] - 3.71748) < eps);

  REQUIRE(std::fabs(rhok[26 + 6] + 0.93151) < eps);
  REQUIRE(std::fabs(rhok[52 + 6] + 2.34518) < eps);
  REQUIRE(std::fabs(rhok[26 + 7] + 0.93151) < eps);
  REQUIRE(std::fabs(rhok[52 + 7] + 2.34518) < eps);
  REQUIRE(std::fabs(rhok[26 + 8] - 2.52341) < eps);
  REQUIRE(std::fabs(rhok[52 + 8] - 0.00000) < eps);

  REQUIRE(std::fabs(rhok[26 + 9] - 2.52341) < eps);
  REQUIRE(std::fabs(rhok[52 + 9] - 0.00000) < eps);
  REQUIRE(std::fabs(rhok[26 + 10] + 0.93151) < eps);
  REQUIRE(std::fabs(rhok[52 + 10] + 2.34518) < eps);
  REQUIRE(std::fabs(rhok[26 + 11] - 2.52341) < eps);
  REQUIRE(std::fabs(rhok[52 + 11] - 0.00000) < eps);

  REQUIRE(std::fabs(rhok[26 + 12] - 2.52341) < eps);
  REQUIRE(std::fabs(rhok[52 + 12] - 0.00000) < eps);
  REQUIRE(std::fabs(rhok[26 + 13] + 0.93151) < eps);
  REQUIRE(std::fabs(rhok[52 + 13] - 2.34518) < eps);
  REQUIRE(std::fabs(rhok[26 + 14] - 2.52341) < eps);
  REQUIRE(std::fabs(rhok[52 + 14] - 0.00000) < eps);

  REQUIRE(std::fabs(rhok[26 + 15] - 2.52341) < eps);
  REQUIRE(std::fabs(rhok[52 + 15] - 0.00000) < eps);
  REQUIRE(std::fabs(rhok[26 + 16] + 0.93151) < eps);
  REQUIRE(std::fabs(rhok[52 + 16] - 2.34518) < eps);
  REQUIRE(std::fabs(rhok[26 + 17] + 0.93151) < eps);
  REQUIRE(std::fabs(rhok[52 + 17] - 2.34518) < eps);

  REQUIRE(std::fabs(rhok[26 + 18] + 1.38359) < eps);
  REQUIRE(std::fabs(rhok[52 + 18] + 0.30687) < eps);
  REQUIRE(std::fabs(rhok[26 + 19] - 0.79595) < eps);
  REQUIRE(std::fabs(rhok[52 + 19] + 1.17259) < eps);
  REQUIRE(std::fabs(rhok[26 + 20] - 0.79595) < eps);
  REQUIRE(std::fabs(rhok[52 + 20] + 1.17259) < eps);

  REQUIRE(std::fabs(rhok[26 + 21] - 0.79595) < eps);
  REQUIRE(std::fabs(rhok[52 + 21] - 1.17259) < eps);
  REQUIRE(std::fabs(rhok[26 + 22] - 0.79595) < eps);
  REQUIRE(std::fabs(rhok[52 + 22] + 1.17259) < eps);
  REQUIRE(std::fabs(rhok[26 + 23] - 0.79595) < eps);
  REQUIRE(std::fabs(rhok[52 + 23] - 1.17259) < eps);

  REQUIRE(std::fabs(rhok[26 + 24] - 0.79595) < eps);
  REQUIRE(std::fabs(rhok[52 + 24] - 1.17259) < eps);
  REQUIRE(std::fabs(rhok[26 + 25] + 1.38359) < eps);
  REQUIRE(std::fabs(rhok[52 + 25] - 0.30687) < eps);

  std::cout << "test_SkAllEstimator:: STOP\n";
}
} // namespace qmcplusplus
