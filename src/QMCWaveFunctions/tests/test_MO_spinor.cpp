//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Configuration.h"
#include "Message/Communicate.h"
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include "Particle/DistanceTableData.h"
#include "QMCWaveFunctions/LCAO/LCAOSpinorBuilder.h"

namespace qmcplusplus
{
void test_lcao_spinor()
{
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  app_log() << "!!!!!!!   LCAO SpinorSet from HDF   !!!!!!!\n";
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

  using ValueType = SPOSet::ValueType;
  using RealType  = SPOSet::RealType;
  Communicate* c;
  c = OHMMS::Controller;

  ParticleSet ions_;
  ParticleSet elec_;

  ions_.setName("ion");
  ions_.create(1);

  ions_.R[0][0]        = 0.00000000;
  ions_.R[0][1]        = 0.00000000;
  ions_.R[0][2]        = 0.00000000;
  SpeciesSet& ispecies = ions_.getSpeciesSet();
  int hIdx             = ispecies.addSpecies("H");
  ions_.update();

  elec_.setName("elec");
  elec_.create(1);
  elec_.R[0][0]  = 0.1;
  elec_.R[0][1]  = -0.3;
  elec_.R[0][2]  = 1.7;
  elec_.spins[0] = 0.6;

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;


  elec_.addTable(ions_, DT_SOA);
  elec_.update();

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.addParticleSet(&elec_);
  ptcl.addParticleSet(&ions_);

  const char* particles = "<tmp> \
   <sposet_builder name=\"spinorbuilder\" type=\"molecularspinor\" href=\"lcao_spinor.h5\" source=\"ion\" precision=\"float\"> \
     <sposet name=\"myspo\" size=\"1\"/> \
   </sposet_builder> \
   </tmp> \
";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr bnode   = xmlFirstElementChild(root);
  xmlNodePtr sponode = xmlFirstElementChild(bnode);

  LCAOSpinorBuilder bb(elec_, ions_, c, bnode);
  std::unique_ptr<SPOSet> spo(bb.createSPOSetFromXML(sponode));
  REQUIRE(spo);

  SPOSet::ValueMatrix_t psiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix_t dpsiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t dspsiM(elec_.R.size(), spo->getOrbitalSetSize()); //spin gradient
  SPOSet::ValueMatrix_t d2psiM(elec_.R.size(), spo->getOrbitalSetSize());

  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, d2psiM);

  ValueType val(5.584596565578567e-05, 0.0012321335483993093);
  ValueType vdx(-2.7922982853738495e-05, 0.0006160667747698895);
  ValueType vdy(8.376894847079449e-05, 0.0018482003223147029);
  ValueType vdz(-0.00047469070804637896, -0.010473135160780798);
  ValueType vlp(0.0033369958606517744, 0.07362437917398082);

  for (unsigned int iat = 0; iat < 1; iat++)
  {
    REQUIRE(psiM[iat][0] == ComplexApprox(val).epsilon(6));
    REQUIRE(dpsiM[iat][0][0] == ComplexApprox(vdx).epsilon(6));
    REQUIRE(dpsiM[iat][0][1] == ComplexApprox(vdy).epsilon(6));
    REQUIRE(dpsiM[iat][0][2] == ComplexApprox(vdz).epsilon(6));
    REQUIRE(d2psiM[iat][0] == ComplexApprox(vlp).epsilon(6));
  }
}

TEST_CASE("ReadMolecularOrbital GTO spinor", "[wavefunction]") { test_lcao_spinor(); }

} // namespace qmcplusplus
