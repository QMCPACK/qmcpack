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
  std::cout << "orb size" << spo->getOrbitalSetSize() << std::endl;

  SPOSet::ValueMatrix_t psiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix_t dpsiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t dspsiM(elec_.R.size(), spo->getOrbitalSetSize()); //spin gradient
  SPOSet::ValueMatrix_t d2psiM(elec_.R.size(), spo->getOrbitalSetSize());

  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, d2psiM);

  ValueType val(5.584596565578567e-05,0.0012321335483993093);
  ValueType vdx(-2.7922982853738495e-05,0.0006160667747698895);
  ValueType vdy(8.376894847079449e-05,0.0018482003223147029);
  ValueType vdz(-0.00047469070804637896,-0.010473135160780798);
  ValueType vlp(0.0033369958606517744,0.07362437917398082);

  for (unsigned int iat = 0; iat < 1; iat++)
  {
    REQUIRE(psiM[iat][0] == ComplexApprox(val).epsilon(6));
    REQUIRE(dpsiM[iat][0][0] == ComplexApprox(vdx).epsilon(6));
    REQUIRE(dpsiM[iat][0][1] == ComplexApprox(vdy).epsilon(6));
    REQUIRE(dpsiM[iat][0][2] == ComplexApprox(vdz).epsilon(6));
    REQUIRE(d2psiM[iat][0] == ComplexApprox(vlp).epsilon(6));
  }
  /*
    REQUIRE(psiM[iat][0] == ComplexApprox(psiM_ref[iat][0]));
    REQUIRE(psiM[iat][1] == ComplexApprox(psiM_ref[iat][1]));
    REQUIRE(psiM[iat][2] == ComplexApprox(psiM_ref[iat][2]));

    REQUIRE(dpsiM[iat][0][0] == ComplexApprox(dpsiM_ref[iat][0][0]));
    REQUIRE(dpsiM[iat][0][1] == ComplexApprox(dpsiM_ref[iat][0][1]));
    REQUIRE(dpsiM[iat][0][2] == ComplexApprox(dpsiM_ref[iat][0][2]));

    REQUIRE(dpsiM[iat][1][0] == ComplexApprox(dpsiM_ref[iat][1][0]));
    REQUIRE(dpsiM[iat][1][1] == ComplexApprox(dpsiM_ref[iat][1][1]));
    REQUIRE(dpsiM[iat][1][2] == ComplexApprox(dpsiM_ref[iat][1][2]));

    REQUIRE(dpsiM[iat][2][0] == ComplexApprox(dpsiM_ref[iat][2][0]));
    REQUIRE(dpsiM[iat][2][1] == ComplexApprox(dpsiM_ref[iat][2][1]));
    REQUIRE(dpsiM[iat][2][2] == ComplexApprox(dpsiM_ref[iat][2][2]));


    REQUIRE(d2psiM[iat][0] == ComplexApprox(d2psiM_ref[iat][0]));
    REQUIRE(d2psiM[iat][1] == ComplexApprox(d2psiM_ref[iat][1]));
    REQUIRE(d2psiM[iat][2] == ComplexApprox(d2psiM_ref[iat][2]));
  }

  //Now we're going to test evaluateValue and evaluateVGL:

  int OrbitalSetSize = spo->getOrbitalSetSize();
  //temporary arrays for holding the values of the up and down channels respectively.
  SPOSet::ValueVector_t psi_work;

  //temporary arrays for holding the gradients of the up and down channels respectively.
  SPOSet::GradVector_t dpsi_work;

  //temporary arrays for holding the laplacians of the up and down channels respectively.
  SPOSet::ValueVector_t d2psi_work;


  psi_work.resize(OrbitalSetSize);

  dpsi_work.resize(OrbitalSetSize);

  d2psi_work.resize(OrbitalSetSize);

  //We worked hard to generate nice reference data above.  Let's generate a test for evaluateV
  //and evaluateVGL by perturbing the electronic configuration by dR, and then make
  //single particle moves that bring it back to our original R reference values.

  //Our perturbation vector.
  ParticleSet::ParticlePos_t dR;
  dR.resize(3);
  //Values chosen based on divine inspiration.  Not important.
  dR[0][0] = 0.1;
  dR[0][1] = 0.2;
  dR[0][2] = 0.1;
  dR[1][0] = -0.1;
  dR[1][1] = -0.05;
  dR[1][2] = 0.05;
  dR[2][0] = -0.1;
  dR[2][1] = 0.0;
  dR[2][2] = 0.0;

  //The new R of our perturbed particle set.
  ParticleSet::ParticlePos_t Rnew;
  Rnew.resize(3);
  Rnew    = elec_.R + dR;
  elec_.R = Rnew;
  elec_.update();

  //Now we test evaluateValue()
  for (unsigned int iat = 0; iat < 3; iat++)
  {
    psi_work = 0.0;
    elec_.makeMove(iat, -dR[iat], false);
    spo->evaluateValue(elec_, iat, psi_work);

    REQUIRE(psi_work[0] == ComplexApprox(psiM_ref[iat][0]));
    REQUIRE(psi_work[1] == ComplexApprox(psiM_ref[iat][1]));
    REQUIRE(psi_work[2] == ComplexApprox(psiM_ref[iat][2]));
    elec_.rejectMove(iat);
  }

  //Now we test evaluateVGL()
  for (unsigned int iat = 0; iat < 3; iat++)
  {
    psi_work   = 0.0;
    dpsi_work  = 0.0;
    d2psi_work = 0.0;

    elec_.makeMove(iat, -dR[iat], false);
    spo->evaluateVGL(elec_, iat, psi_work, dpsi_work, d2psi_work);

    REQUIRE(psi_work[0] == ComplexApprox(psiM_ref[iat][0]));
    REQUIRE(psi_work[1] == ComplexApprox(psiM_ref[iat][1]));
    REQUIRE(psi_work[2] == ComplexApprox(psiM_ref[iat][2]));

    REQUIRE(dpsi_work[0][0] == ComplexApprox(dpsiM_ref[iat][0][0]));
    REQUIRE(dpsi_work[0][1] == ComplexApprox(dpsiM_ref[iat][0][1]));
    REQUIRE(dpsi_work[0][2] == ComplexApprox(dpsiM_ref[iat][0][2]));

    REQUIRE(dpsi_work[1][0] == ComplexApprox(dpsiM_ref[iat][1][0]));
    REQUIRE(dpsi_work[1][1] == ComplexApprox(dpsiM_ref[iat][1][1]));
    REQUIRE(dpsi_work[1][2] == ComplexApprox(dpsiM_ref[iat][1][2]));

    REQUIRE(dpsi_work[2][0] == ComplexApprox(dpsiM_ref[iat][2][0]));
    REQUIRE(dpsi_work[2][1] == ComplexApprox(dpsiM_ref[iat][2][1]));
    REQUIRE(dpsi_work[2][2] == ComplexApprox(dpsiM_ref[iat][2][2]));

    REQUIRE(d2psi_work[0] == ComplexApprox(d2psiM_ref[iat][0]));
    REQUIRE(d2psi_work[1] == ComplexApprox(d2psiM_ref[iat][1]));
    REQUIRE(d2psi_work[2] == ComplexApprox(d2psiM_ref[iat][2]));
    elec_.rejectMove(iat);
  }

  //Now we test evaluateSpin:
  SPOSet::ValueVector_t dspsi_work;
  dspsi_work.resize(OrbitalSetSize);

  for (unsigned int iat = 0; iat < 3; iat++)
  {
    psi_work   = 0.0;
    dspsi_work = 0.0;

    elec_.makeMove(iat, -dR[iat], false);
    spo->evaluate_spin(elec_, iat, psi_work, dspsi_work);

    REQUIRE(psi_work[0] == ComplexApprox(psiM_ref[iat][0]));
    REQUIRE(psi_work[1] == ComplexApprox(psiM_ref[iat][1]));
    REQUIRE(psi_work[2] == ComplexApprox(psiM_ref[iat][2]));

    REQUIRE(dspsi_work[0] == ComplexApprox(dspsiM_ref[iat][0]));
    REQUIRE(dspsi_work[1] == ComplexApprox(dspsiM_ref[iat][1]));
    REQUIRE(dspsi_work[2] == ComplexApprox(dspsiM_ref[iat][2]));

    elec_.rejectMove(iat);
  }
    */
}

TEST_CASE("ReadMolecularOrbital GTO spinor", "[wavefunction]") { test_lcao_spinor(); }

} // namespace qmcplusplus
