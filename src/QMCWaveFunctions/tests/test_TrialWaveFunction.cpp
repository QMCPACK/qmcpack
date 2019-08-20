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
#include "Particle/ParticleSet.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminant.h"
#include "QMCWaveFunctions/Jastrow/RadialJastrowBuilder.h"


namespace qmcplusplus
{

TEST_CASE("TrialWaveFunction", "[wavefunction]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate* c = OHMMS::Controller;

  ParticleSet ions_;
  ParticleSet elec_;

  ions_.setName("ion");
  ions_.create(2);
  ions_.R[0][0] = 0.0;
  ions_.R[0][1] = 0.0;
  ions_.R[0][2] = 0.0;
  ions_.R[1][0] = 1.68658058;
  ions_.R[1][1] = 1.68658058;
  ions_.R[1][2] = 1.68658058;
  ions_.update();


  elec_.setName("elec");
  elec_.create(4);
  elec_.R[0] = {0.0, 0.0, 0.0};
  elec_.R[1] = {0.0, 1.0, 0.0};
  elec_.R[2] = {1.0, 0.0, 0.0};
  elec_.R[3] = {0.0, 0.0, 1.0};

  // diamondC_1x1x1
  elec_.Lattice.R(0, 0) = 3.37316115;
  elec_.Lattice.R(0, 1) = 3.37316115;
  elec_.Lattice.R(0, 2) = 0.0;
  elec_.Lattice.R(1, 0) = 0.0;
  elec_.Lattice.R(1, 1) = 3.37316115;
  elec_.Lattice.R(1, 2) = 3.37316115;
  elec_.Lattice.R(2, 0) = 3.37316115;
  elec_.Lattice.R(2, 1) = 0.0;
  elec_.Lattice.R(2, 2) = 3.37316115;

  SpeciesSet& tspecies         = elec_.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;

#ifdef ENABLE_SOA
  elec_.addTable(ions_, DT_SOA);
#else
  elec_.addTable(ions_, DT_AOS);
#endif
  elec_.resetGroups();
  elec_.update();

  ParticleSetPool ptcl{c};
  ptcl.addParticleSet(&elec_);
  ptcl.addParticleSet(&ions_);

  //diamondC_1x1x1
  const char* particles = "<tmp> \
<determinantset type=\"einspline\" href=\"pwscf.pwscf.h5\" tilematrix=\"1 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\" meshfactor=\"1.0\" precision=\"float\" size=\"4\"/> \
</tmp> \
";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr spo_root = doc.getRoot();
  xmlNodePtr ein1 = xmlFirstElementChild(spo_root);

  EinsplineSetBuilder einSet(elec_, ptcl.getPool(), c, ein1);
  SPOSet* spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo != nullptr);

  auto* det = new DiracDeterminant<>(spo);
  det->set(0, 4);

  TrialWaveFunction psi{c};
  psi.addComponent(det, "SingleDet");
  const char* jas_input = "<tmp> \
<jastrow name=\"J2\" type=\"Two-Body\" function=\"Bspline\" print=\"yes\"> \
   <correlation size=\"10\" speciesA=\"u\" speciesB=\"d\"> \
      <coefficients id=\"ud\" type=\"Array\"> 0.02904699284 -0.1004179 -0.1752703883 -0.2232576505 -0.2728029201 -0.3253286875 -0.3624525145 -0.3958223107 -0.4268582166 -0.4394531176</coefficients> \
    </correlation> \
</jastrow> \
</tmp> \
";
  Libxml2Document doc_jas;
  okay = doc.parseFromString(jas_input);
  REQUIRE(okay);

  xmlNodePtr jas_root = doc.getRoot();
  xmlNodePtr jas1 = xmlFirstElementChild(jas_root);

  RadialJastrowBuilder jb(elec_, psi);
  bool build_okay = jb.put(jas1);
/*
  REQUIRE(build_okay);
*/
#if !defined(QMC_CUDA)
  double logpsi = psi.evaluateLog(elec_);
  std::cout << std::setprecision(20) << logpsi << "debug YYYY" << std::endl;
  REQUIRE(logpsi == Approx(0.9023637518659));
#endif

}

} // namespace qmcplusplus
