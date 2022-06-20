//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminant.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "QMCWaveFunctions/Jastrow/RadialJastrowBuilder.h"
#include "TWFGrads.hpp"
#include <ResourceCollection.h>

namespace qmcplusplus
{
#if defined(ENABLE_CUDA) && !defined(QMC_CUDA2HIP)
using DiracDet = DiracDeterminant<DelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
#else
using DiracDet = DiracDeterminant<DelayedUpdate<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
#endif

using LogValueType = TrialWaveFunction::LogValueType;
using PsiValueType = TrialWaveFunction::PsiValueType;
using GradType     = TrialWaveFunction::GradType;

TEST_CASE("TrialWaveFunction_diamondC_1x1x1", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

#if defined(ENABLE_OFFLOAD)
  const DynamicCoordinateKind kind_selected = DynamicCoordinateKind::DC_POS_OFFLOAD;
#else
  const DynamicCoordinateKind kind_selected = DynamicCoordinateKind::DC_POS;
#endif
  // diamondC_1x1x1
  ParticleSet::ParticleLayout lattice;
  lattice.R(0, 0)   = 3.37316115;
  lattice.R(0, 1)   = 3.37316115;
  lattice.R(0, 2)   = 0.0;
  lattice.R(1, 0)   = 0.0;
  lattice.R(1, 1)   = 3.37316115;
  lattice.R(1, 2)   = 3.37316115;
  lattice.R(2, 0)   = 3.37316115;
  lattice.R(2, 1)   = 0.0;
  lattice.R(2, 2)   = 3.37316115;
  lattice.BoxBConds = {1, 1, 1};
  lattice.reset();

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.setSimulationCell(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell(), kind_selected);
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell(), kind_selected);
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({2});
  ions_.R[0] = {0.0, 0.0, 0.0};
  ions_.R[1] = {1.68658058, 1.68658058, 1.68658058};
  ions_.update();


  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({2, 2});
  elec_.R[0] = {0.0, 0.0, 0.0};
  elec_.R[1] = {0.0, 1.0, 1.0};
  elec_.R[2] = {1.0, 1.0, 0.0};
  elec_.R[3] = {1.0, 0.0, 1.0};

  SpeciesSet& tspecies         = elec_.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;

  elec_.addTable(ions_);
  elec_.resetGroups();
  elec_.createSK(); // needed by AoS J2 for ChiesaKEcorrection

  // make a ParticleSet Clone
  ParticleSet elec_clone(elec_);

  //diamondC_1x1x1
  const char* spo_xml = "<tmp> \
<determinantset type=\"einspline\" href=\"diamondC_1x1x1.pwscf.h5\" tilematrix=\"1 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\" meshfactor=\"1.0\" precision=\"float\" size=\"2\"/> \
</tmp> \
";

  Libxml2Document doc;
  bool okay = doc.parseFromString(spo_xml);
  REQUIRE(okay);

  xmlNodePtr spo_root = doc.getRoot();
  xmlNodePtr ein1     = xmlFirstElementChild(spo_root);

  EinsplineSetBuilder einSet(elec_, ptcl.getPool(), c, ein1);
  auto spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo != nullptr);

  std::vector<std::unique_ptr<DiracDeterminantBase>> dets;
  dets.push_back(std::make_unique<DiracDet>(spo->makeClone(), 0, 2));
  dets.push_back(std::make_unique<DiracDet>(spo->makeClone(), 2, 4));

  auto slater_det = std::make_unique<SlaterDet>(elec_, std::move(dets));

  TrialWaveFunction psi;
  psi.addComponent(std::move(slater_det));

  const char* jas_input = "<tmp> \
<jastrow name=\"J2\" type=\"Two-Body\" function=\"Bspline\" print=\"yes\"> \
   <correlation size=\"10\" speciesA=\"u\" speciesB=\"u\"> \
      <coefficients id=\"uu\" type=\"Array\"> 0.02904699284 -0.1004179 -0.1752703883 -0.2232576505 -0.2728029201 -0.3253286875 -0.3624525145 -0.3958223107 -0.4268582166 -0.4394531176</coefficients> \
   </correlation> \
</jastrow> \
</tmp> \
";
  Libxml2Document doc_jas;
  okay = doc.parseFromString(jas_input);
  REQUIRE(okay);

  xmlNodePtr jas_root = doc.getRoot();
  xmlNodePtr jas1     = xmlFirstElementChild(jas_root);

  RadialJastrowBuilder jb(c, elec_);
  psi.addComponent(jb.buildComponent(jas1));

#if !defined(QMC_CUDA)
  // initialize distance tables.
  elec_.update();
  double logpsi = psi.evaluateLog(elec_);

  //std::cout << "debug before YYY " << std::setprecision(16) << psi.getLogPsi() << " " << psi.getPhase()<< std::endl;
#if defined(QMC_COMPLEX)
  REQUIRE(logpsi == Approx(-0.1201465271523596));
#else
  REQUIRE(logpsi == Approx(-1.471840358291562));
#endif

  // make a TrialWaveFunction Clone
  std::unique_ptr<TrialWaveFunction> psi_clone(psi.makeClone(elec_clone));

  elec_clone.update();
  double logpsi_clone = psi_clone->evaluateLog(elec_clone);
#if defined(QMC_COMPLEX)
  REQUIRE(logpsi_clone == Approx(-0.1201465271523596));
#else
  REQUIRE(logpsi_clone == Approx(-1.471840358291562));
#endif

  const int moved_elec_id = 0;

  using PosType   = QMCTraits::PosType;
  using RealType  = QMCTraits::RealType;
  using ValueType = QMCTraits::ValueType;
  PosType delta(0.1, 0.1, 0.2);

  elec_.makeMove(moved_elec_id, delta);

  ValueType r_all_val       = psi.calcRatio(elec_, moved_elec_id);
  ValueType r_fermionic_val = psi.calcRatio(elec_, moved_elec_id, TrialWaveFunction::ComputeType::FERMIONIC);
  ValueType r_bosonic_val   = psi.calcRatio(elec_, moved_elec_id, TrialWaveFunction::ComputeType::NONFERMIONIC);

  //std::cout << "debug YYY " << std::setprecision(16) << r_all_val << std::endl;
  //std::cout << "debug YYY " << std::setprecision(16) << r_fermionic_val << std::endl;
  //std::cout << "debug YYY " << std::setprecision(16) << r_bosonic_val << std::endl;
#if defined(QMC_COMPLEX)
  CHECK(r_all_val == ComplexApprox(ValueType(1.653821746120792, 0.5484992491019633)));
  CHECK(r_fermionic_val == ComplexApprox(ValueType(1.804065087219802, 0.598328295048828)));
#else
  CHECK(r_all_val == Approx(2.305591774210242));
  CHECK(r_fermionic_val == ValueApprox(2.515045914101833));
#endif
  CHECK(r_bosonic_val == ValueApprox(0.9167195562048454));

  psi.acceptMove(elec_, moved_elec_id);
  elec_.acceptMove(moved_elec_id);
#if defined(QMC_COMPLEX)
  REQUIRE(psi.getLogPsi() == Approx(0.4351202455204972));
#else
  REQUIRE(psi.getLogPsi() == Approx(-0.63650297977845492));
#endif

  elec_.update(true);
  psi.evaluateLog(elec_);
#if defined(QMC_COMPLEX)
  REQUIRE(psi.getLogPsi() == Approx(0.4351202455204972));
#else
  REQUIRE(psi.getLogPsi() == Approx(-0.63650297977845492));
#endif

  // testing batched interfaces
  ResourceCollection pset_res("test_pset_res");
  ResourceCollection twf_res("test_twf_res");

  elec_.createResource(pset_res);
  psi.createResource(twf_res);

  //Temporary as switch to std::reference_wrapper proceeds
  // testing batched interfaces
  RefVectorWithLeader<ParticleSet> p_ref_list(elec_, {elec_, elec_clone});
  RefVectorWithLeader<TrialWaveFunction> wf_ref_list(psi, {psi, *psi_clone});

  ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(pset_res, p_ref_list);
  ResourceCollectionTeamLock<TrialWaveFunction> mw_twf_lock(twf_res, wf_ref_list);

  ParticleSet::mw_update(p_ref_list);
  TrialWaveFunction::mw_evaluateLog(wf_ref_list, p_ref_list);
#if defined(QMC_COMPLEX)
  REQUIRE(std::complex<RealType>(wf_ref_list[0].getLogPsi(), wf_ref_list[0].getPhase()) ==
          LogComplexApprox(std::complex<RealType>(0.4351202455204972, 6.665972664860828)));
  REQUIRE(std::complex<RealType>(wf_ref_list[1].getLogPsi(), wf_ref_list[1].getPhase()) ==
          LogComplexApprox(std::complex<RealType>(-0.1201465271523596, 6.345732826640545)));
#else
  REQUIRE(std::complex<RealType>(wf_ref_list[0].getLogPsi(), wf_ref_list[0].getPhase()) ==
          LogComplexApprox(std::complex<RealType>(-0.6365029797784554, 3.141592653589793)));
  REQUIRE(std::complex<RealType>(wf_ref_list[1].getLogPsi(), wf_ref_list[1].getPhase()) ==
          LogComplexApprox(std::complex<RealType>(-1.471840358291562, 3.141592653589793)));
#endif

  TWFGrads<CoordsType::POS> grad_old(2);

  grad_old.grads_positions[0] = wf_ref_list[0].evalGrad(p_ref_list[0], moved_elec_id);
  grad_old.grads_positions[1] = wf_ref_list[1].evalGrad(p_ref_list[1], moved_elec_id);

  std::cout << "evalGrad " << std::setprecision(14) << grad_old.grads_positions[0][0] << " "
            << grad_old.grads_positions[0][1] << " " << grad_old.grads_positions[0][2] << " "
            << grad_old.grads_positions[1][0] << " " << grad_old.grads_positions[1][1] << " "
            << grad_old.grads_positions[1][2] << std::endl;

  TrialWaveFunction::mw_evalGrad(wf_ref_list, p_ref_list, moved_elec_id, grad_old);
#if defined(QMC_COMPLEX)
  CHECK(grad_old.grads_positions[0][0] == ComplexApprox(ValueType(18.817970466022, -6.5837500306076)));
  CHECK(grad_old.grads_positions[0][1] == ComplexApprox(ValueType(-22.840838391977, 3.9963373883645)));
  CHECK(grad_old.grads_positions[0][2] == ComplexApprox(ValueType(3.8805320617146, 1.5825508129169)));
  CHECK(grad_old.grads_positions[1][0] == ComplexApprox(ValueType(47.387717528888, -8.7703065253151e-06)));
  CHECK(grad_old.grads_positions[1][1] == ComplexApprox(ValueType(-54.671696901113, -7.3126138879524)));
  CHECK(grad_old.grads_positions[1][2] == ComplexApprox(ValueType(6.6288917088321, 7.3126230586018)));
#else
  CHECK(grad_old.grads_positions[0][0] == Approx(14.77249702264));
  CHECK(grad_old.grads_positions[0][1] == Approx(-20.385235323777));
  CHECK(grad_old.grads_positions[0][2] == Approx(4.8529516184558));
  CHECK(grad_old.grads_positions[1][0] == Approx(47.38770710732));
  CHECK(grad_old.grads_positions[1][1] == Approx(-63.361119579044));
  CHECK(grad_old.grads_positions[1][2] == Approx(15.318325284049));
#endif

  PosType delta_sign_changed(0.1, 0.1, -0.2);

  std::vector<PosType> displs{delta_sign_changed, delta_sign_changed};
  ParticleSet::mw_makeMove(p_ref_list, moved_elec_id, displs);

  if (kind_selected != DynamicCoordinateKind::DC_POS_OFFLOAD)
  {
    ValueType r_0 = wf_ref_list[0].calcRatio(p_ref_list[0], moved_elec_id);
    GradType grad_temp;
    ValueType r_1 = wf_ref_list[1].calcRatioGrad(p_ref_list[1], moved_elec_id, grad_temp);
#if defined(QMC_COMPLEX)
    CHECK(r_0 == ComplexApprox(ValueType(-0.045474407700114, -0.59956233350555)));
    CHECK(r_1 == ComplexApprox(ValueType(-0.44602867091608, -1.8105588403509)));
    CHECK(grad_temp[0] == ComplexApprox(ValueType(-6.6139971152489, 22.82304260002)));
    CHECK(grad_temp[1] == ComplexApprox(ValueType(8.3367501707711, -23.362154838104)));
    CHECK(grad_temp[2] == ComplexApprox(ValueType(-2.6347597529645, 0.67383144279783)));
#else
    CHECK(r_0 == Approx(-0.4138835449));
    CHECK(r_1 == Approx(-2.5974770159));
    CHECK(grad_temp[0] == Approx(-17.865723259764));
    CHECK(grad_temp[1] == Approx(19.854257889369));
    CHECK(grad_temp[2] == Approx(-2.9669578650441));
#endif
  }

  PosType delta_zero(0, 0, 0);
  displs = {delta_zero, delta};
  ParticleSet::mw_makeMove(p_ref_list, moved_elec_id, displs);

  std::vector<PsiValueType> ratios(2);
  TrialWaveFunction::mw_calcRatio(wf_ref_list, p_ref_list, moved_elec_id, ratios);
  std::cout << "calcRatio " << std::setprecision(14) << ratios[0] << " " << ratios[1] << std::endl;
#if defined(QMC_COMPLEX)
  CHECK(ratios[0] == ComplexApprox(PsiValueType(1, 0)));
  CHECK(ratios[1] == ComplexApprox(PsiValueType(1.6538214581548, 0.54849918598717)));
#else
  CHECK(ratios[0] == Approx(1));
  CHECK(ratios[1] == Approx(2.3055913093424));
#endif

  std::fill(ratios.begin(), ratios.end(), 0);
  TWFGrads<CoordsType::POS> grad_new(2);

  if (kind_selected != DynamicCoordinateKind::DC_POS_OFFLOAD)
  {
    ratios[0] = wf_ref_list[0].calcRatioGrad(p_ref_list[0], moved_elec_id, grad_new.grads_positions[0]);
    ratios[1] = wf_ref_list[1].calcRatioGrad(p_ref_list[1], moved_elec_id, grad_new.grads_positions[1]);

    std::cout << "calcRatioGrad " << std::setprecision(14) << ratios[0] << " " << ratios[1] << std::endl
              << grad_new.grads_positions[0][0] << " " << grad_new.grads_positions[0][1] << " "
              << grad_new.grads_positions[0][2] << " " << grad_new.grads_positions[1][0] << " "
              << grad_new.grads_positions[1][1] << " " << grad_new.grads_positions[1][2] << std::endl;
  }
  //Temporary as switch to std::reference_wrapper proceeds
  // testing batched interfaces

  TrialWaveFunction::mw_calcRatioGrad(wf_ref_list, p_ref_list, moved_elec_id, ratios, grad_new);
#if defined(QMC_COMPLEX)
  CHECK(ratios[0] == ComplexApprox(ValueType(1, 0)));
  CHECK(grad_new.grads_positions[0][0] == ComplexApprox(ValueType(18.817970466022, -6.5837500306076)));
  CHECK(grad_new.grads_positions[0][1] == ComplexApprox(ValueType(-22.840838391977, 3.9963373883645)));
  CHECK(grad_new.grads_positions[0][2] == ComplexApprox(ValueType(3.8805320617146, 1.5825508129169)));
  CHECK(ratios[1] == ComplexApprox(ValueType(1.6538214581548, 0.54849918598717)));
  CHECK(grad_new.grads_positions[1][0] == ComplexApprox(ValueType(18.817970466022, -6.5837500306076)));
  CHECK(grad_new.grads_positions[1][1] == ComplexApprox(ValueType(-22.840838391977, 3.9963373883645)));
  CHECK(grad_new.grads_positions[1][2] == ComplexApprox(ValueType(3.8805320617146, 1.5825508129169)));
#else
  CHECK(ratios[0] == Approx(1));
  CHECK(grad_new.grads_positions[0][0] == Approx(14.77249702264));
  CHECK(grad_new.grads_positions[0][1] == Approx(-20.385235323777));
  CHECK(grad_new.grads_positions[0][2] == Approx(4.8529516184558));
  CHECK(ratios[1] == Approx(2.3055913093424));
  CHECK(grad_new.grads_positions[1][0] == Approx(14.77249702264));
  CHECK(grad_new.grads_positions[1][1] == Approx(-20.385235323777));
  CHECK(grad_new.grads_positions[1][2] == Approx(4.8529516184558));
#endif

  std::vector<bool> isAccepted(2, true);
  TrialWaveFunction::mw_accept_rejectMove(wf_ref_list, p_ref_list, moved_elec_id, isAccepted);
#if defined(QMC_COMPLEX)
  REQUIRE(std::complex<RealType>(wf_ref_list[0].getLogPsi(), wf_ref_list[0].getPhase()) ==
          LogComplexApprox(std::complex<RealType>(0.4351202455204972, 6.665972664860828)));
  REQUIRE(std::complex<RealType>(wf_ref_list[1].getLogPsi(), wf_ref_list[1].getPhase()) ==
          LogComplexApprox(std::complex<RealType>(0.4351202455204972, 6.665972664860828)));
#else
  REQUIRE(std::complex<RealType>(wf_ref_list[0].getLogPsi(), wf_ref_list[0].getPhase()) ==
          LogComplexApprox(std::complex<RealType>(-0.6365029797784554, 3.141592653589793)));
  REQUIRE(std::complex<RealType>(wf_ref_list[1].getLogPsi(), wf_ref_list[1].getPhase()) ==
          LogComplexApprox(std::complex<RealType>(-0.6365029797784554, 3.141592653589793)));
#endif

  TrialWaveFunction::mw_evalGrad(wf_ref_list, p_ref_list, moved_elec_id, grad_old);
#if defined(QMC_COMPLEX)
  CHECK(grad_old.grads_positions[0][0] == ComplexApprox(ValueType(18.817970466022, -6.5837500306076)));
  CHECK(grad_old.grads_positions[0][1] == ComplexApprox(ValueType(-22.840838391977, 3.9963373883645)));
  CHECK(grad_old.grads_positions[0][2] == ComplexApprox(ValueType(3.8805320617146, 1.5825508129169)));
  CHECK(grad_old.grads_positions[1][0] == ComplexApprox(ValueType(18.817970466022, -6.5837500306076)));
  CHECK(grad_old.grads_positions[1][1] == ComplexApprox(ValueType(-22.840838391977, 3.9963373883645)));
  CHECK(grad_old.grads_positions[1][2] == ComplexApprox(ValueType(3.8805320617146, 1.5825508129169)));
#else
  CHECK(grad_old.grads_positions[0][0] == Approx(14.77249702264));
  CHECK(grad_old.grads_positions[0][1] == Approx(-20.385235323777));
  CHECK(grad_old.grads_positions[0][2] == Approx(4.8529516184558));
  CHECK(grad_old.grads_positions[1][0] == Approx(14.77249702264));
  CHECK(grad_old.grads_positions[1][1] == Approx(-20.385235323777));
  CHECK(grad_old.grads_positions[1][2] == Approx(4.8529516184558));
#endif

#endif
}


} // namespace qmcplusplus
