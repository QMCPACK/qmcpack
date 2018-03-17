//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Configuration.h"
#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCWaveFunctions/WaveFunctionFactory.h"

namespace qmcplusplus
{

TEST_CASE("WaveFunctionFactory","[wavefunction]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate *c = OHMMS::Controller;

  ParticleSet *qp = new ParticleSet;
  std::vector<int> agroup(1);
  agroup[0] = 2;
  qp->setName("e");
  qp->create(agroup);
  qp->R[0][0] = 1.0;
  qp->R[0][1] = 2.0;
  qp->R[0][2] = 3.0;
  qp->R[1][0] = 0.0;
  qp->R[1][1] = 1.1;
  qp->R[1][2] = 2.2;

  SpeciesSet &tspecies = qp->getSpeciesSet();
  int upIdx = tspecies.addSpecies("u");
  int downIdx = tspecies.addSpecies("d");
  int massIdx = tspecies.addAttribute("mass");
  tspecies(massIdx, upIdx) = 1.0;
  tspecies(massIdx, downIdx) = 1.0;

  qp->update();

  WaveFunctionFactory::PtclPoolType particle_set_map;
  particle_set_map["e"] = qp;


  WaveFunctionFactory wff(qp, particle_set_map, c);

  const char* wavefunction_xml = \
"<wavefunction> \
         <jastrow type=\"Two-Body\" name=\"J2\" function=\"bspline\" print=\"yes\"> \
            <correlation speciesA=\"u\" speciesB=\"d\" size=\"8\" cutoff=\"10.0\"> \
               <coefficients id=\"ud\" type=\"Array\"> \
0.5954603818 0.5062051797 0.3746940461 0.2521010502 0.1440163317 0.07796688253 \
0.03804420551 0.01449320872 \
               </coefficients> \
            </correlation> \
         </jastrow> \
</wavefunction>";
  Libxml2Document doc;
  bool okay = doc.parseFromString(wavefunction_xml);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();
  wff.put(root);

  REQUIRE(wff.targetPsi != NULL);
  REQUIRE(wff.targetPsi->size() == 1);

  OrbitalBase *j2_base = wff.targetPsi->getOrbitals()[0];
  REQUIRE(j2_base != NULL);
}
}

