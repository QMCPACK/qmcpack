//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
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
#include "Particle/DistanceTable.h"
#include "QMCWaveFunctions/SPOSetBuilderFactory.h"
#include "Utilities/ResourceCollection.h"

namespace qmcplusplus
{
void test_lcao_spinor()
{
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  app_log() << "!!!!!!!   LCAO SpinorSet from HDF   !!!!!!!\n";
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

  using ValueType = SPOSet::ValueType;
  using RealType  = SPOSet::RealType;
  Communicate* c = OHMMS::Controller;

  ParticleSetPool ptcl = ParticleSetPool(c);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create(1);

  ions_.R[0][0]        = 0.00000000;
  ions_.R[0][1]        = 0.00000000;
  ions_.R[0][2]        = 0.00000000;
  SpeciesSet& ispecies = ions_.getSpeciesSet();
  int hIdx             = ispecies.addSpecies("H");
  ions_.update();

  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create(1);
  elec_.R[0][0]  = 0.1;
  elec_.R[0][1]  = -0.3;
  elec_.R[0][2]  = 1.7;
  elec_.spins[0] = 0.6;
  elec_.setSpinor(true);

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;


  elec_.addTable(ions_);
  elec_.update();

  const char* particles = "<tmp> \
   <sposet_builder name=\"spinorbuilder\" type=\"molecularorbital\" href=\"lcao_spinor.h5\" source=\"ion\" precision=\"float\"> \
     <basisset transform=\"yes\"/> \
     <sposet name=\"myspo\" size=\"1\"/> \
   </sposet_builder> \
   </tmp> \
";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr bnode = xmlFirstElementChild(root);
  SPOSetBuilderFactory fac(c, elec_, ptcl.getPool());
  auto& bb = fac.createSPOSetBuilder(bnode);

  // only pick up the last sposet
  SPOSet* spo = nullptr;
  processChildren(bnode, [&](const std::string& cname, const xmlNodePtr element) {
    if (cname == "sposet")
      spo = bb.createSPOSet(element);
  });
  REQUIRE(spo);

  SPOSet::ValueMatrix_t psiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix_t dpsiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t dspsiM(elec_.R.size(), spo->getOrbitalSetSize()); //spin gradient
  SPOSet::ValueMatrix_t d2psiM(elec_.R.size(), spo->getOrbitalSetSize());

  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, d2psiM);

  ValueType val(5.584596565578567e-05, 0.0012321335483993093);
  ValueType vdx(-2.7922982853738495e-05, -0.0006160667747698895);
  ValueType vdy(8.376894847079449e-05, 0.0018482003223147029);
  ValueType vdz(-0.00047469070804637896, -0.010473135160780798);
  ValueType vlp(0.0033369958606517744, 0.07362437917398082);
  ValueType vds(1.0474021389417806e-05, -0.00040482519442528657);
  const RealType eps = 1e-4;

  for (unsigned int iat = 0; iat < 1; iat++)
  {
    CHECK(psiM[iat][0] == ComplexApprox(val).epsilon(eps));
    CHECK(dpsiM[iat][0][0] == ComplexApprox(vdx).epsilon(eps));
    CHECK(dpsiM[iat][0][1] == ComplexApprox(vdy).epsilon(eps));
    CHECK(dpsiM[iat][0][2] == ComplexApprox(vdz).epsilon(eps));
    CHECK(d2psiM[iat][0] == ComplexApprox(vlp).epsilon(eps));
  }

  int OrbitalSetSize = spo->getOrbitalSetSize();
  //temporary arrays for holding the values of the up and down channels respectively.
  SPOSet::ValueVector_t psi_work;

  //temporary arrays for holding the gradients of the up and down channels respectively.
  SPOSet::GradVector_t dpsi_work;

  //temporary arrays for holding the laplacians of the up and down channels respectively.
  SPOSet::ValueVector_t d2psi_work;

  //temporary arrays for holding spin gradient
  SPOSet::ValueVector_t dspsi_work;

  psi_work.resize(OrbitalSetSize);
  dpsi_work.resize(OrbitalSetSize);
  d2psi_work.resize(OrbitalSetSize);
  dspsi_work.resize(OrbitalSetSize);

  //We worked hard to generate nice reference data above.  Let's generate a test for evaluateV
  //and evaluateVGL by perturbing the electronic configuration by dR, and then make
  //single particle moves that bring it back to our original R reference values.

  //Our perturbation vector.
  ParticleSet::ParticlePos_t dR;
  dR.resize(1);
  dR[0][0] = 0.1;
  dR[0][1] = 0.2;
  dR[0][2] = 0.1;

  //The new R of our perturbed particle set. Ma
  ParticleSet::ParticlePos_t Rnew;
  Rnew.resize(1);
  Rnew    = elec_.R + dR;
  elec_.R = Rnew;
  elec_.update();

  //Now we test evaluateValue()
  for (unsigned int iat = 0; iat < 1; iat++)
  {
    psi_work = 0.0;
    elec_.makeMove(iat, -dR[iat], false);
    spo->evaluateValue(elec_, iat, psi_work);

    CHECK(psi_work[0] == ComplexApprox(val));
    elec_.rejectMove(iat);
  }

  //Now we test evaluateVGL()
  for (unsigned int iat = 0; iat < 1; iat++)
  {
    psi_work   = 0.0;
    dpsi_work  = 0.0;
    d2psi_work = 0.0;
    dspsi_work = 0.0;

    elec_.makeMove(iat, -dR[iat], false);
    spo->evaluateVGL_spin(elec_, iat, psi_work, dpsi_work, d2psi_work, dspsi_work);

    CHECK(psi_work[0] == ComplexApprox(val).epsilon(eps));
    CHECK(dpsi_work[0][0] == ComplexApprox(vdx).epsilon(eps));
    CHECK(dpsi_work[0][1] == ComplexApprox(vdy).epsilon(eps));
    CHECK(dpsi_work[0][2] == ComplexApprox(vdz).epsilon(eps));
    CHECK(d2psi_work[0] == ComplexApprox(vlp).epsilon(eps));
    CHECK(dspsi_work[0] == ComplexApprox(vds).epsilon(eps));
    elec_.rejectMove(iat);
  }

  //Now we test evaluateSpin:

  for (unsigned int iat = 0; iat < 1; iat++)
  {
    psi_work   = 0.0;
    dspsi_work = 0.0;

    elec_.makeMove(iat, -dR[iat], false);
    spo->evaluate_spin(elec_, iat, psi_work, dspsi_work);

    CHECK(psi_work[0] == ComplexApprox(val).epsilon(eps));
    CHECK(dspsi_work[0] == ComplexApprox(vds).epsilon(eps));

    elec_.rejectMove(iat);
  }

  // test batched interface
  // first move elec_ back to original positions for reference
  Rnew    = elec_.R - dR;
  elec_.R = Rnew;
  elec_.update();

  //now create second walker
  ParticleSet elec_2(elec_);
  elec_2.R[0][0]  = -0.4;
  elec_2.R[0][1]  = 1.5;
  elec_2.R[0][2]  = -0.2;
  elec_2.spins[0] = -1.3;

  //create new reference values for new positions
  ValueType val2(-0.00010787670610075059, -5.882498404872149e-05);
  ValueType vdx2(-0.0002157534121903495, -0.00011764996809136147);
  ValueType vdy2(0.0008090752956289096, 0.0004411873802963092);
  ValueType vdz2(-0.00010787670612158852, -5.8824984060083926e-05);
  ValueType vlp2(-0.004989237947754119, -0.0027206229528162103);
  ValueType vds2(0.0001907917398151183, 0.005002478563410625);

  ResourceCollection pset_res("test_pset_res");
  elec_.createResource(pset_res);
  RefVectorWithLeader<ParticleSet> p_list(elec_);
  p_list.push_back(elec_);
  p_list.push_back(elec_2);

  ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(pset_res, p_list);

  elec_.mw_update(p_list);
  std::unique_ptr<SPOSet> spo_2(spo->makeClone());
  RefVectorWithLeader<SPOSet> spo_list(*spo);
  spo_list.push_back(*spo);
  spo_list.push_back(*spo_2);

  SPOSet::ValueMatrix_t psiM_2(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix_t dpsiM_2(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t d2psiM_2(elec_.R.size(), spo->getOrbitalSetSize());

  RefVector<SPOSet::ValueMatrix_t> logdet_list;
  RefVector<SPOSet::GradMatrix_t> dlogdet_list;
  RefVector<SPOSet::ValueMatrix_t> d2logdet_list;

  logdet_list.push_back(psiM);
  logdet_list.push_back(psiM_2);
  dlogdet_list.push_back(dpsiM);
  dlogdet_list.push_back(dpsiM_2);
  d2logdet_list.push_back(d2psiM);
  d2logdet_list.push_back(d2psiM_2);

  spo->mw_evaluate_notranspose(spo_list, p_list, 0, 1, logdet_list, dlogdet_list, d2logdet_list);
  for (unsigned int iat = 0; iat < 1; iat++)
  {
    //walker 0
    CHECK(logdet_list[0].get()[iat][0] == ComplexApprox(val).epsilon(eps));
    CHECK(dlogdet_list[0].get()[iat][0][0] == ComplexApprox(vdx).epsilon(eps));
    CHECK(dlogdet_list[0].get()[iat][0][1] == ComplexApprox(vdy).epsilon(eps));
    CHECK(dlogdet_list[0].get()[iat][0][2] == ComplexApprox(vdz).epsilon(eps));
    CHECK(d2logdet_list[0].get()[iat][0] == ComplexApprox(vlp).epsilon(eps));
    //walker 1
    CHECK(logdet_list[1].get()[iat][0] == ComplexApprox(val2).epsilon(eps));
    CHECK(dlogdet_list[1].get()[iat][0][0] == ComplexApprox(vdx2).epsilon(eps));
    CHECK(dlogdet_list[1].get()[iat][0][1] == ComplexApprox(vdy2).epsilon(eps));
    CHECK(dlogdet_list[1].get()[iat][0][2] == ComplexApprox(vdz2).epsilon(eps));
    CHECK(d2logdet_list[1].get()[iat][0] == ComplexApprox(vlp2).epsilon(eps));
  }

  //first, lets displace all the elec in each walker
  for (int iat = 0; iat < 1; iat++)
  {
    std::vector<ParticleSet::SingleParticlePos_t> displs = {dR[iat], dR[iat]};
    elec_.mw_makeMove(p_list, iat, displs);
    std::vector<bool> accept = {true, true};
    elec_.mw_accept_rejectMove(p_list, iat, accept);
  }
  elec_.mw_update(p_list);

  SPOSet::ValueVector_t psi_work_2(OrbitalSetSize);
  SPOSet::GradVector_t dpsi_work_2(OrbitalSetSize);
  SPOSet::ValueVector_t d2psi_work_2(OrbitalSetSize);
  SPOSet::ValueVector_t dspsi_work_2(OrbitalSetSize);

  RefVector<SPOSet::ValueVector_t> psi_v_list   = {psi_work, psi_work_2};
  RefVector<SPOSet::GradVector_t> dpsi_v_list   = {dpsi_work, dpsi_work_2};
  RefVector<SPOSet::ValueVector_t> d2psi_v_list = {d2psi_work, d2psi_work_2};
  RefVector<SPOSet::ValueVector_t> dspsi_v_list = {dspsi_work, dspsi_work_2};
  //check mw_evaluateVGLWithSpin
  for (int iat = 0; iat < 1; iat++)
  {
    //reset values to zero, updates the ref vectors to zero as well
    psi_work     = 0.0;
    dpsi_work    = 0.0;
    d2psi_work   = 0.0;
    dspsi_work   = 0.0;
    psi_work_2   = 0.0;
    dpsi_work_2  = 0.0;
    d2psi_work_2 = 0.0;
    dspsi_work_2 = 0.0;

    std::vector<ParticleSet::SingleParticlePos_t> displs = {-dR[iat], -dR[iat]};
    elec_.mw_makeMove(p_list, iat, displs);
    spo->mw_evaluateVGLWithSpin(spo_list, p_list, iat, psi_v_list, dpsi_v_list, d2psi_v_list, dspsi_v_list);
    //walker 0
    CHECK(psi_v_list[0].get()[0] == ComplexApprox(val).epsilon(eps));
    CHECK(dpsi_v_list[0].get()[0][0] == ComplexApprox(vdx).epsilon(eps));
    CHECK(dpsi_v_list[0].get()[0][1] == ComplexApprox(vdy).epsilon(eps));
    CHECK(dpsi_v_list[0].get()[0][2] == ComplexApprox(vdz).epsilon(eps));
    CHECK(d2psi_v_list[0].get()[0] == ComplexApprox(vlp).epsilon(eps));
    CHECK(dspsi_v_list[0].get()[0] == ComplexApprox(vds).epsilon(eps));
    //walker 1
    CHECK(psi_v_list[1].get()[0] == ComplexApprox(val2).epsilon(eps));
    CHECK(dpsi_v_list[1].get()[0][0] == ComplexApprox(vdx2).epsilon(eps));
    CHECK(dpsi_v_list[1].get()[0][1] == ComplexApprox(vdy2).epsilon(eps));
    CHECK(dpsi_v_list[1].get()[0][2] == ComplexApprox(vdz2).epsilon(eps));
    CHECK(d2psi_v_list[1].get()[0] == ComplexApprox(vlp2).epsilon(eps));
    CHECK(dspsi_v_list[1].get()[0] == ComplexApprox(vds2).epsilon(eps));

    std::vector<bool> accept = {false, false};
    elec_.mw_accept_rejectMove(p_list, iat, accept);
  }
}

void test_lcao_spinor_excited()
{
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  app_log() << "!! LCAO SpinorSet from HDF with excited  !!\n";
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

  using ValueType = SPOSet::ValueType;
  using RealType  = SPOSet::RealType;
  Communicate* c = OHMMS::Controller;

  ParticleSetPool ptcl = ParticleSetPool(c);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create(1);

  ions_.R[0][0]        = 0.00000000;
  ions_.R[0][1]        = 0.00000000;
  ions_.R[0][2]        = 0.00000000;
  SpeciesSet& ispecies = ions_.getSpeciesSet();
  int hIdx             = ispecies.addSpecies("H");
  ions_.update();

  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create(1);
  elec_.R[0][0]  = 0.1;
  elec_.R[0][1]  = -0.3;
  elec_.R[0][2]  = 1.7;
  elec_.spins[0] = 0.6;
  elec_.setSpinor(true);

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;


  elec_.addTable(ions_);
  elec_.update();

  const char* particles = "<tmp> \
   <sposet_builder name=\"spinorbuilder\" type=\"molecularorbital\" href=\"lcao_spinor.h5\" source=\"ion\" precision=\"float\"> \
     <basisset name=\"myset\" transform=\"yes\"/> \
     <sposet name=\"myspo\" basisset=\"myset\" size=\"1\"> \
       <occupation mode=\"excited\"> \
         -1 2 \
       </occupation> \
     </sposet> \
   </sposet_builder> \
   </tmp> \
";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr bnode = xmlFirstElementChild(root);
  SPOSetBuilderFactory fac(c, elec_, ptcl.getPool());
  auto& bb = fac.createSPOSetBuilder(bnode);

  // only pick up the last sposet
  SPOSet* spo = nullptr;
  processChildren(bnode, [&](const std::string& cname, const xmlNodePtr element) {
    if (cname == "sposet")
      spo = bb.createSPOSet(element);
  });
  REQUIRE(spo);

  const int OrbitalSetSize = spo->getOrbitalSetSize();
  CHECK(OrbitalSetSize == 1);

  SPOSet::ValueMatrix_t psiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix_t dpsiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t dspsiM(elec_.R.size(), spo->getOrbitalSetSize()); //spin gradient
  SPOSet::ValueMatrix_t d2psiM(elec_.R.size(), spo->getOrbitalSetSize());

  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, d2psiM);

  ValueType val(0.0008237860500019983, 1.0474021389417806e-05);
  ValueType vdx(-0.00041189302538224967, -5.237010699556317e-06);
  ValueType vdy(0.0012356790748129446, 1.5711032081710294e-05);
  ValueType vdz(-0.007002181424606922, -8.90291818048377e-05);
  ValueType vlp(0.04922415803252472, 0.0006258601782677606);
  ValueType vds(-0.0010017050778321178, -5.584596565578559e-05);
  const RealType eps = 1e-4;

  for (unsigned int iat = 0; iat < 1; iat++)
  {
    CHECK(psiM[iat][0] == ComplexApprox(val).epsilon(eps));
    CHECK(dpsiM[iat][0][0] == ComplexApprox(vdx).epsilon(eps));
    CHECK(dpsiM[iat][0][1] == ComplexApprox(vdy).epsilon(eps));
    CHECK(dpsiM[iat][0][2] == ComplexApprox(vdz).epsilon(eps));
    CHECK(d2psiM[iat][0] == ComplexApprox(vlp).epsilon(eps));
  }

  //temporary arrays for holding the values of the up and down channels respectively.
  SPOSet::ValueVector_t psi_work;

  //temporary arrays for holding the gradients of the up and down channels respectively.
  SPOSet::GradVector_t dpsi_work;

  //temporary arrays for holding the laplacians of the up and down channels respectively.
  SPOSet::ValueVector_t d2psi_work;

  //temporary arrays for holding spin gradient
  SPOSet::ValueVector_t dspsi_work;

  psi_work.resize(OrbitalSetSize);
  dpsi_work.resize(OrbitalSetSize);
  d2psi_work.resize(OrbitalSetSize);
  dspsi_work.resize(OrbitalSetSize);

  //We worked hard to generate nice reference data above.  Let's generate a test for evaluateV
  //and evaluateVGL by perturbing the electronic configuration by dR, and then make
  //single particle moves that bring it back to our original R reference values.

  //Our perturbation vector.
  ParticleSet::ParticlePos_t dR;
  dR.resize(1);
  dR[0][0] = 0.1;
  dR[0][1] = 0.2;
  dR[0][2] = 0.1;

  //The new R of our perturbed particle set. Ma
  ParticleSet::ParticlePos_t Rnew;
  Rnew.resize(1);
  Rnew    = elec_.R + dR;
  elec_.R = Rnew;
  elec_.update();

  //Now we test evaluateValue()
  for (unsigned int iat = 0; iat < 1; iat++)
  {
    psi_work = 0.0;
    elec_.makeMove(iat, -dR[iat], false);
    spo->evaluateValue(elec_, iat, psi_work);

    CHECK(psi_work[0] == ComplexApprox(val));
    elec_.rejectMove(iat);
  }

  //Now we test evaluateVGL()
  for (unsigned int iat = 0; iat < 1; iat++)
  {
    psi_work   = 0.0;
    dpsi_work  = 0.0;
    d2psi_work = 0.0;
    dspsi_work = 0.0;

    elec_.makeMove(iat, -dR[iat], false);
    spo->evaluateVGL_spin(elec_, iat, psi_work, dpsi_work, d2psi_work, dspsi_work);

    CHECK(psi_work[0] == ComplexApprox(val).epsilon(eps));
    CHECK(dpsi_work[0][0] == ComplexApprox(vdx).epsilon(eps));
    CHECK(dpsi_work[0][1] == ComplexApprox(vdy).epsilon(eps));
    CHECK(dpsi_work[0][2] == ComplexApprox(vdz).epsilon(eps));
    CHECK(d2psi_work[0] == ComplexApprox(vlp).epsilon(eps));
    CHECK(dspsi_work[0] == ComplexApprox(vds).epsilon(eps));
    elec_.rejectMove(iat);
  }

  //Now we test evaluateSpin:
  for (unsigned int iat = 0; iat < 1; iat++)
  {
    psi_work   = 0.0;
    dspsi_work = 0.0;

    elec_.makeMove(iat, -dR[iat], false);
    spo->evaluate_spin(elec_, iat, psi_work, dspsi_work);

    CHECK(psi_work[0] == ComplexApprox(val).epsilon(eps));
    CHECK(dspsi_work[0] == ComplexApprox(vds).epsilon(eps));

    elec_.rejectMove(iat);
  }

  //test batched interface
  // first move elec_ back to orginal positions for reference
  Rnew    = elec_.R - dR;
  elec_.R = Rnew;
  elec_.update();

  //now create second walker
  ParticleSet elec_2(elec_);
  elec_2.R[0][0]  = -0.4;
  elec_2.R[0][1]  = 1.5;
  elec_2.R[0][2]  = -0.2;
  elec_2.spins[0] = -1.3;

  //create new reference values for new positions
  ValueType val2(0.0026291910291941015, 0.00019079173981511807);
  ValueType vdx2(0.005258382058116388, 0.00038158347961051147);
  ValueType vdy2(-0.019718932715867252, -0.0014309380483892627);
  ValueType vdz2(0.002629191029701947, 0.00019079173985197097);
  ValueType vlp2(0.1215986298515522, 0.008824012363842379);
  ValueType vds2(0.004256243259981321, 0.00010787670610075102);

  ResourceCollection pset_res("test_pset_res");
  elec_.createResource(pset_res);
  RefVectorWithLeader<ParticleSet> p_list(elec_);
  p_list.push_back(elec_);
  p_list.push_back(elec_2);

  ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(pset_res, p_list);

  elec_.mw_update(p_list);
  std::unique_ptr<SPOSet> spo_2(spo->makeClone());
  RefVectorWithLeader<SPOSet> spo_list(*spo);
  spo_list.push_back(*spo);
  spo_list.push_back(*spo_2);

  SPOSet::ValueMatrix_t psiM_2(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix_t dpsiM_2(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t d2psiM_2(elec_.R.size(), spo->getOrbitalSetSize());

  RefVector<SPOSet::ValueMatrix_t> logdet_list;
  RefVector<SPOSet::GradMatrix_t> dlogdet_list;
  RefVector<SPOSet::ValueMatrix_t> d2logdet_list;

  logdet_list.push_back(psiM);
  logdet_list.push_back(psiM_2);
  dlogdet_list.push_back(dpsiM);
  dlogdet_list.push_back(dpsiM_2);
  d2logdet_list.push_back(d2psiM);
  d2logdet_list.push_back(d2psiM_2);

  spo->mw_evaluate_notranspose(spo_list, p_list, 0, 1, logdet_list, dlogdet_list, d2logdet_list);
  for (unsigned int iat = 0; iat < 1; iat++)
  {
    //walker 0
    CHECK(logdet_list[0].get()[iat][0] == ComplexApprox(val).epsilon(eps));
    CHECK(dlogdet_list[0].get()[iat][0][0] == ComplexApprox(vdx).epsilon(eps));
    CHECK(dlogdet_list[0].get()[iat][0][1] == ComplexApprox(vdy).epsilon(eps));
    CHECK(dlogdet_list[0].get()[iat][0][2] == ComplexApprox(vdz).epsilon(eps));
    CHECK(d2logdet_list[0].get()[iat][0] == ComplexApprox(vlp).epsilon(eps));
    //walker 1
    CHECK(logdet_list[1].get()[iat][0] == ComplexApprox(val2).epsilon(eps));
    CHECK(dlogdet_list[1].get()[iat][0][0] == ComplexApprox(vdx2).epsilon(eps));
    CHECK(dlogdet_list[1].get()[iat][0][1] == ComplexApprox(vdy2).epsilon(eps));
    CHECK(dlogdet_list[1].get()[iat][0][2] == ComplexApprox(vdz2).epsilon(eps));
    CHECK(d2logdet_list[1].get()[iat][0] == ComplexApprox(vlp2).epsilon(eps));
  }

  //first, lets displace all the elec in each walker
  for (int iat = 0; iat < 1; iat++)
  {
    std::vector<ParticleSet::SingleParticlePos_t> displs = {dR[iat], dR[iat]};
    elec_.mw_makeMove(p_list, iat, displs);
    std::vector<bool> accept = {true, true};
    elec_.mw_accept_rejectMove(p_list, iat, accept);
  }
  elec_.mw_update(p_list);

  SPOSet::ValueVector_t psi_work_2(OrbitalSetSize);
  SPOSet::GradVector_t dpsi_work_2(OrbitalSetSize);
  SPOSet::ValueVector_t d2psi_work_2(OrbitalSetSize);
  SPOSet::ValueVector_t dspsi_work_2(OrbitalSetSize);

  RefVector<SPOSet::ValueVector_t> psi_v_list   = {psi_work, psi_work_2};
  RefVector<SPOSet::GradVector_t> dpsi_v_list   = {dpsi_work, dpsi_work_2};
  RefVector<SPOSet::ValueVector_t> d2psi_v_list = {d2psi_work, d2psi_work_2};
  RefVector<SPOSet::ValueVector_t> dspsi_v_list = {dspsi_work, dspsi_work_2};
  //check mw_evaluateVGLWithSpin
  for (int iat = 0; iat < 1; iat++)
  {
    //reset values to zero, updates the ref vectors to zero as well
    psi_work     = 0.0;
    dpsi_work    = 0.0;
    d2psi_work   = 0.0;
    dspsi_work   = 0.0;
    psi_work_2   = 0.0;
    dpsi_work_2  = 0.0;
    d2psi_work_2 = 0.0;
    dspsi_work_2 = 0.0;

    std::vector<ParticleSet::SingleParticlePos_t> displs = {-dR[iat], -dR[iat]};
    elec_.mw_makeMove(p_list, iat, displs);
    spo->mw_evaluateVGLWithSpin(spo_list, p_list, iat, psi_v_list, dpsi_v_list, d2psi_v_list, dspsi_v_list);
    //walker 0
    CHECK(psi_v_list[0].get()[0] == ComplexApprox(val).epsilon(eps));
    CHECK(dpsi_v_list[0].get()[0][0] == ComplexApprox(vdx).epsilon(eps));
    CHECK(dpsi_v_list[0].get()[0][1] == ComplexApprox(vdy).epsilon(eps));
    CHECK(dpsi_v_list[0].get()[0][2] == ComplexApprox(vdz).epsilon(eps));
    CHECK(d2psi_v_list[0].get()[0] == ComplexApprox(vlp).epsilon(eps));
    CHECK(dspsi_v_list[0].get()[0] == ComplexApprox(vds).epsilon(eps));
    //walker 1
    CHECK(psi_v_list[1].get()[0] == ComplexApprox(val2).epsilon(eps));
    CHECK(dpsi_v_list[1].get()[0][0] == ComplexApprox(vdx2).epsilon(eps));
    CHECK(dpsi_v_list[1].get()[0][1] == ComplexApprox(vdy2).epsilon(eps));
    CHECK(dpsi_v_list[1].get()[0][2] == ComplexApprox(vdz2).epsilon(eps));
    CHECK(d2psi_v_list[1].get()[0] == ComplexApprox(vlp2).epsilon(eps));
    CHECK(dspsi_v_list[1].get()[0] == ComplexApprox(vds2).epsilon(eps));

    std::vector<bool> accept = {false, false};
    elec_.mw_accept_rejectMove(p_list, iat, accept);
  }
}


TEST_CASE("ReadMolecularOrbital GTO spinor", "[wavefunction]") { test_lcao_spinor(); }
TEST_CASE("ReadMolecularOrbital GTO spinor with excited", "[wavefunction]") { test_lcao_spinor_excited(); }

} // namespace qmcplusplus
