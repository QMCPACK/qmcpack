//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include "WaveFunctionFactory.h"
#include "LCAO/LCAOrbitalSet.h"
#include "TWFGrads.hpp"
#include "Utilities/RuntimeOptions.h"
#include <ResourceCollection.h>

#include <stdio.h>
#include <string>
#include <limits>

using std::string;

namespace qmcplusplus
{
using PosType      = ParticleSet::PosType;
using RealType     = ParticleSet::RealType;
using ValueType    = ParticleSet::ValueType;
using GradType     = ParticleSet::GradType;
using LogValueType = WaveFunctionComponent::LogValueType;
using PsiValueType = WaveFunctionComponent::PsiValueType;

void test_LiH_msd(const std::string& spo_xml_string,
                  const std::string& check_sponame,
                  int check_spo_size,
                  int check_basisset_size,
                  int test_nlpp_algorithm_batched,
                  int test_batched_api)
{
  Communicate* c = OHMMS::Controller;

  ParticleSetPool ptcl = ParticleSetPool(c);
  auto ions_uptr       = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr       = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion0");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({1, 1});
  ions_.R[0]           = {0.0, 0.0, 0.0};
  ions_.R[1]           = {0.0, 0.0, 3.0139239693};
  SpeciesSet& ispecies = ions_.getSpeciesSet();
  int LiIdx            = ispecies.addSpecies("Li");
  int HIdx             = ispecies.addSpecies("H");

  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({2, 2});
  elec_.R[0] = {0.5, 0.5, 0.5};
  elec_.R[1] = {0.1, 0.1, 1.1};
  elec_.R[2] = {-0.5, -0.5, -0.5};
  elec_.R[3] = {-0.1, -0.1, 1.5};

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int downIdx                = tspecies.addSpecies("d");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(massIdx, upIdx)   = 1.0;
  tspecies(massIdx, downIdx) = 1.0;
  // Necessary to set mass
  elec_.resetGroups();

  Libxml2Document doc;
  bool okay = doc.parseFromString(spo_xml_string);
  REQUIRE(okay);

  xmlNodePtr ein_xml = doc.getRoot();

  WaveFunctionFactory wf_factory(elec_, ptcl.getPool(), c);
  RuntimeOptions runtime_options;
  auto twf_ptr = wf_factory.buildTWF(ein_xml, runtime_options);

  auto& spo = dynamic_cast<const LCAOrbitalSet&>(twf_ptr->getSPOSet(check_sponame));
  CHECK(spo.getOrbitalSetSize() == check_spo_size);
  CHECK(spo.getBasisSetSize() == check_basisset_size);

  ions_.update();
  elec_.update();

  auto& twf(*twf_ptr);
  twf.setMassTerm(elec_);
  twf.evaluateLog(elec_);

  app_log() << "twf.evaluateLog logpsi " << std::setprecision(16) << twf.getLogPsi() << " " << twf.getPhase()
            << std::endl;
  CHECK(std::complex<double>(twf.getLogPsi(), twf.getPhase()) ==
        LogComplexApprox(std::complex<double>(-7.646027846242066, 3.141592653589793)));
  CHECK(elec_.G[0][0] == ValueApprox(-2.181896934));
  CHECK(elec_.G[1][1] == ValueApprox(0.120821033));
  CHECK(elec_.G[2][2] == ValueApprox(1.2765987657));
  CHECK(elec_.L[0] == ValueApprox(-15.460736911));
  CHECK(elec_.L[3] == ValueApprox(-0.328013327566));

  twf.prepareGroup(elec_, 0);
  auto grad_old = twf.evalGrad(elec_, 1);
  app_log() << "twf.evalGrad grad_old " << std::setprecision(16) << grad_old << std::endl;
  CHECK(grad_old[0] == ValueApprox(0.1204183219));
  CHECK(grad_old[1] == ValueApprox(0.120821033));
  CHECK(grad_old[2] == ValueApprox(2.05904174));

  PosType delta(0.1, 0.1, 0.2);
  elec_.makeMove(1, delta);

  ParticleSet::GradType grad_new;
  auto ratio = twf.calcRatioGrad(elec_, 1, grad_new);
  app_log() << "twf.calcRatioGrad ratio " << ratio << " grad_new " << grad_new << std::endl;
  CHECK(ratio == ValueApprox(1.374307585));
  CHECK(grad_new[0] == ValueApprox(0.05732804333));
  CHECK(grad_new[1] == ValueApprox(0.05747775029));
  CHECK(grad_new[2] == ValueApprox(1.126889742));

  ratio = twf.calcRatio(elec_, 1);
  app_log() << "twf.calcRatio ratio " << ratio << std::endl;
  CHECK(ratio == ValueApprox(1.374307585));


  opt_variables_type active;
  twf.checkInVariables(active);

  const int nparam = active.size_of_active();
  REQUIRE(nparam == 1486);

  using ValueType = QMCTraits::ValueType;
  Vector<ValueType> dlogpsi(nparam);
  Vector<ValueType> dhpsioverpsi(nparam);
  twf.evaluateDerivatives(elec_, active, dlogpsi, dhpsioverpsi);

  // Numbers not validated
  CHECK(dlogpsi[0] == ValueApprox(0.006449058893092842));
  CHECK(dlogpsi[1] == ValueApprox(-0.01365690177395768));
  CHECK(dlogpsi[nparam - 1] == ValueApprox(0.1641910574099575));

  CHECK(dhpsioverpsi[0] == ValueApprox(0.2207480131794138));
  CHECK(dhpsioverpsi[1] == ValueApprox(0.009316665149067847));
  CHECK(dhpsioverpsi[nparam - 1] == ValueApprox(0.982665984797896));

  if (test_nlpp_algorithm_batched)
  {
    // set virtutal particle position
    PosType newpos(0.3, 0.2, 0.5);

    //elec_.makeVirtualMoves(newpos);
    //std::vector<ValueType> ratios(elec_.getTotalNum());
    //twf.evaluateRatiosAlltoOne(elec_, ratios);

    //CHECK(std::real(ratios[0]) == Approx());
    //CHECK(std::real(ratios[1]) == Approx());
    //CHECK(std::real(ratios[2]) == Approx());

    elec_.makeMove(0, newpos - elec_.R[0]);
    ValueType ratio_0 = twf.calcRatio(elec_, 0);
    elec_.rejectMove(0);

    CHECK(std::real(ratio_0) == Approx(2.350046921));

    VirtualParticleSet VP(elec_, 2);
    std::vector<PosType> newpos2(2);
    std::vector<ValueType> ratios2(2);
    newpos2[0] = newpos - elec_.R[1];
    newpos2[1] = PosType(0.2, 0.5, 0.3) - elec_.R[1];
    VP.makeMoves(elec_, 1, newpos2);
    twf.evaluateRatios(VP, ratios2);

    CHECK(std::real(ratios2[0]) == Approx(-0.8544310407));
    CHECK(std::real(ratios2[1]) == Approx(-1.0830708458));

    std::fill(ratios2.begin(), ratios2.end(), 0);
    Matrix<ValueType> dratio(2, nparam);
    twf.evaluateDerivRatios(VP, active, ratios2, dratio);

    CHECK(std::real(ratios2[0]) == Approx(-0.8544310407));
    CHECK(std::real(ratios2[1]) == Approx(-1.0830708458));

    CHECK(std::real(dratio[0][0]) == Approx(0.248887465));
    CHECK(std::real(dratio[0][1]) == Approx(0.135021218));
  }

  //test acceptMove
  {
    PosType newpos(0.3, 0.2, 0.5);
    elec_.makeMove(1, newpos - elec_.R[1]);
    ValueType ratio_1 = twf.calcRatio(elec_, 1);
    twf.acceptMove(elec_, 1);
    elec_.acceptMove(1);

    CHECK(std::real(ratio_1) == Approx(-0.8544310407));
    CHECK(twf.getLogPsi() == Approx(-7.8033473273));

    twf.evaluateLog(elec_);

    app_log() << "twf.evaluateLog logpsi " << std::setprecision(16) << twf.getLogPsi() << " " << twf.getPhase()
              << std::endl;
    CHECK(std::complex<double>(twf.getLogPsi(), twf.getPhase()) ==
          LogComplexApprox(std::complex<double>(-7.803347327300154, 0.0)));
    CHECK(elec_.G[0][0] == ValueApprox(1.63020975849953));
    CHECK(elec_.G[1][1] == ValueApprox(-1.795375999646262));
    CHECK(elec_.G[2][2] == ValueApprox(1.215768958589418));
    CHECK(elec_.L[0] == ValueApprox(-21.84021387509693));
    CHECK(elec_.L[3] == ValueApprox(-1.332448295858972));
  }

  // testing batched interfaces
  if (test_batched_api)
  {
    ParticleSet elec_clone(elec_);
    elec_clone.update();

    std::unique_ptr<TrialWaveFunction> twf_clone(twf.makeClone(elec_clone));

    std::vector<PosType> displ(2);
    displ[0] = {0.1, 0.2, 0.3};
    displ[1] = {-0.2, -0.3, 0.0};

    // testing batched interfaces
    ResourceCollection pset_res("test_pset_res");
    ResourceCollection twf_res("test_twf_res");

    elec_.createResource(pset_res);
    twf.createResource(twf_res);

    // testing batched interfaces
    RefVectorWithLeader<ParticleSet> p_ref_list(elec_, {elec_, elec_clone});
    RefVectorWithLeader<TrialWaveFunction> wf_ref_list(twf, {twf, *twf_clone});

    ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(pset_res, p_ref_list);
    ResourceCollectionTeamLock<TrialWaveFunction> mw_twf_lock(twf_res, wf_ref_list);

    ParticleSet::mw_update(p_ref_list);
    TrialWaveFunction::mw_evaluateLog(wf_ref_list, p_ref_list);
    app_log() << "before YYY [0] getLogPsi getPhase " << std::setprecision(16) << wf_ref_list[0].getLogPsi() << " "
              << wf_ref_list[0].getPhase() << std::endl;
    app_log() << "before YYY [1] getLogPsi getPhase " << std::setprecision(16) << wf_ref_list[1].getLogPsi() << " "
              << wf_ref_list[1].getPhase() << std::endl;
    CHECK(std::complex<RealType>(wf_ref_list[0].getLogPsi(), wf_ref_list[0].getPhase()) ==
          LogComplexApprox(std::complex<RealType>(-7.803347327300153, 0.0)));
    CHECK(std::complex<RealType>(wf_ref_list[1].getLogPsi(), wf_ref_list[1].getPhase()) ==
          LogComplexApprox(std::complex<RealType>(-7.803347327300153, 0.0)));

    TrialWaveFunction::mw_prepareGroup(wf_ref_list, p_ref_list, 0);

    TWFGrads<CoordsType::POS> grad_old(2);

    const int moved_elec_id = 1;
    TrialWaveFunction::mw_evalGrad(wf_ref_list, p_ref_list, moved_elec_id, grad_old);

    CHECK(grad_old.grads_positions[0][0] == ValueApprox(-2.6785305398));
    CHECK(grad_old.grads_positions[0][1] == ValueApprox(-1.7953759996));
    CHECK(grad_old.grads_positions[0][2] == ValueApprox(-5.8209379274));
    CHECK(grad_old.grads_positions[1][0] == ValueApprox(-2.6785305398));
    CHECK(grad_old.grads_positions[1][1] == ValueApprox(-1.7953759996));
    CHECK(grad_old.grads_positions[1][2] == ValueApprox(-5.8209379274));

    TWFGrads<CoordsType::POS> grad_new(2);
    std::vector<PsiValueType> ratios(2);

    ParticleSet::mw_makeMove(p_ref_list, moved_elec_id, displ);
    TrialWaveFunction::mw_calcRatio(wf_ref_list, p_ref_list, moved_elec_id, ratios);

    CHECK(ratios[0] == ValueApprox(PsiValueType(-0.6181619459)));
    CHECK(ratios[1] == ValueApprox(PsiValueType(1.6186330488)));

    TrialWaveFunction::mw_calcRatioGrad(wf_ref_list, p_ref_list, moved_elec_id, ratios, grad_new);

    CHECK(ratios[0] == ValueApprox(PsiValueType(-0.6181619459)));
    CHECK(ratios[1] == ValueApprox(PsiValueType(1.6186330488)));

    CHECK(grad_new.grads_positions[0][0] == ValueApprox(1.2418467899));
    CHECK(grad_new.grads_positions[0][1] == ValueApprox(1.2425653495));
    CHECK(grad_new.grads_positions[0][2] == ValueApprox(4.4273237873));
    CHECK(grad_new.grads_positions[1][0] == ValueApprox(-0.8633778143));
    CHECK(grad_new.grads_positions[1][1] == ValueApprox(0.8245347691));
    CHECK(grad_new.grads_positions[1][2] == ValueApprox(-5.1513380151));

    std::vector<bool> isAccepted{false, true};
    TrialWaveFunction::mw_accept_rejectMove(wf_ref_list, p_ref_list, moved_elec_id, isAccepted);
    ParticleSet::mw_accept_rejectMove(p_ref_list, moved_elec_id, isAccepted);

    CHECK(wf_ref_list[0].getLogPsi() == Approx(-7.803347327300152));
    CHECK(wf_ref_list[1].getLogPsi() == Approx(-7.321765331299484));

    // move the next electron
    TrialWaveFunction::mw_prepareGroup(wf_ref_list, p_ref_list, 1);
    const int moved_elec_id_next = 2;
    TrialWaveFunction::mw_evalGrad(wf_ref_list, p_ref_list, moved_elec_id_next, grad_old);

    CHECK(grad_old.grads_positions[0][0] == ValueApprox(1.3325558736));
    CHECK(grad_old.grads_positions[0][1] == ValueApprox(1.3327966725));
    CHECK(grad_old.grads_positions[0][2] == ValueApprox(1.2157689586));
    CHECK(grad_old.grads_positions[1][0] == ValueApprox(1.3222514142));
    CHECK(grad_old.grads_positions[1][1] == ValueApprox(1.3230108868));
    CHECK(grad_old.grads_positions[1][2] == ValueApprox(1.2035047435));

    ParticleSet::mw_makeMove(p_ref_list, moved_elec_id_next, displ);
    TrialWaveFunction::mw_calcRatio(wf_ref_list, p_ref_list, moved_elec_id_next, ratios);

    CHECK(ratios[0] == ValueApprox(PsiValueType(2.1080036144)));
    CHECK(ratios[1] == ValueApprox(PsiValueType(0.4947158435)));

    TrialWaveFunction::mw_calcRatioGrad(wf_ref_list, p_ref_list, moved_elec_id_next, ratios, grad_new);

    CHECK(ratios[0] == ValueApprox(PsiValueType(2.1080036144)));
    CHECK(ratios[1] == ValueApprox(PsiValueType(0.4947158435)));

    CHECK(grad_new.grads_positions[0][0] == ValueApprox(1.8412365668));
    CHECK(grad_new.grads_positions[0][1] == ValueApprox(1.3736370007));
    CHECK(grad_new.grads_positions[0][2] == ValueApprox(0.8043818454));
    CHECK(grad_new.grads_positions[1][0] == ValueApprox(1.3553132105));
    CHECK(grad_new.grads_positions[1][1] == ValueApprox(1.5552132255));
    CHECK(grad_new.grads_positions[1][2] == ValueApprox(0.804301246));
  }
}

TEST_CASE("LiH multi Slater dets table_method", "[wavefunction]")
{
  app_log() << "-----------------------------------------------------------------" << std::endl;
  app_log() << "LiH_msd using the table method no precomputation" << std::endl;
  app_log() << "-----------------------------------------------------------------" << std::endl;
  const char* spo_xml_string1 = R"(<wavefunction name="psi0" target="e">
    <sposet_collection type="MolecularOrbital" name="LCAOBSet" source="ion0" cuspCorrection="no" href="LiH.orbs.h5">
      <basisset name="LCAOBSet" key="GTO" transform="yes">
        <grid type="log" ri="1.e-6" rf="1.e2" npts="1001"/>
      </basisset>
      <sposet basisset="LCAOBSet" name="spo-up" size="85">
        <occupation mode="ground"/>
        <coefficient size="85" spindataset="0"/>
      </sposet>
      <sposet basisset="LCAOBSet" name="spo-dn" size="85">
        <occupation mode="ground"/>
        <coefficient size="85" spindataset="0"/>
      </sposet>
    </sposet_collection>
    <determinantset>
      <multideterminant optimize="yes" spo_up="spo-up" spo_dn="spo-dn" algorithm="table_method">
        <detlist size="1487" type="DETS" cutoff="1e-20" href="LiH.orbs.h5"/>
      </multideterminant>
    </determinantset>
</wavefunction>
)";
  test_LiH_msd(spo_xml_string1, "spo-up", 85, 105, true, true);
}

TEST_CASE("LiH multi Slater dets precomputed_table_method", "[wavefunction]")
{
  app_log() << "-----------------------------------------------------------------" << std::endl;
  app_log() << "LiH_msd using the table method with new optimization" << std::endl;
  app_log() << "-----------------------------------------------------------------" << std::endl;
  const char* spo_xml_string1_new = R"(<wavefunction name="psi0" target="e">
    <sposet_collection type="MolecularOrbital" name="LCAOBSet" source="ion0" cuspCorrection="no" href="LiH.orbs.h5">
      <basisset name="LCAOBSet" key="GTO" transform="yes">
        <grid type="log" ri="1.e-6" rf="1.e2" npts="1001"/>
      </basisset>
      <sposet basisset="LCAOBSet" name="spo-up" size="85"> 
        <occupation mode="ground"/>
        <coefficient size="85" spindataset="0"/>
      </sposet>
      <sposet basisset="LCAOBSet" name="spo-dn" size="85">
        <occupation mode="ground"/>
        <coefficient size="85" spindataset="0"/>
      </sposet>
    </sposet_collection>
    <determinantset>
      <multideterminant optimize="yes" spo_up="spo-up" spo_dn="spo-dn" algorithm="precomputed_table_method">
        <detlist size="1487" type="DETS" cutoff="1e-20" href="LiH.orbs.h5"/>
      </multideterminant>
    </determinantset>
</wavefunction>
)";
  test_LiH_msd(spo_xml_string1_new, "spo-up", 85, 105, true, true);
}

#ifdef QMC_COMPLEX
void test_Bi_msd(const std::string& spo_xml_string,
                 const std::string& check_sponame,
                 int check_spo_size,
                 int check_basisset_size)
{
  Communicate* c = OHMMS::Controller;

  ParticleSetPool ptcl = ParticleSetPool(c);
  auto ions_uptr       = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr       = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion0");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create(std::vector<int>{1});
  ions_.R[0]           = {0.0, 0.0, 0.0};
  SpeciesSet& ispecies = ions_.getSpeciesSet();
  int LiIdx            = ispecies.addSpecies("Bi");

  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create(std::vector<int>{5});
  elec_.R[0] = {1.592992772, -2.241313928, -0.7315193518};
  elec_.R[1] = {0.07621077199, 0.8497557547, 1.604678718};
  elec_.R[2] = {2.077473445, 0.680621113, -0.5251243321};
  elec_.R[3] = {-1.488849594, 0.7470552741, 0.6659555498};
  elec_.R[4] = {-1.448485879, 0.7337274141, 0.02687190951};

  elec_.spins[0] = 4.882003828;
  elec_.spins[1] = 0.06469299507;
  elec_.spins[2] = 5.392168887;
  elec_.spins[3] = 5.33941214;
  elec_.spins[4] = 3.127416326;
  elec_.setSpinor(true);

  SpeciesSet& tspecies     = elec_.getSpeciesSet();
  int upIdx                = tspecies.addSpecies("u");
  int massIdx              = tspecies.addAttribute("mass");
  tspecies(massIdx, upIdx) = 1.0;
  // Necessary to set mass
  elec_.resetGroups();

  Libxml2Document doc;
  bool okay = doc.parseFromString(spo_xml_string);
  REQUIRE(okay);

  xmlNodePtr ein_xml = doc.getRoot();

  RuntimeOptions runtime_options;
  WaveFunctionFactory wf_factory(elec_, ptcl.getPool(), c);
  auto twf_ptr = wf_factory.buildTWF(ein_xml, runtime_options);

  auto& spo = twf_ptr->getSPOSet(check_sponame);
  CHECK(spo.getOrbitalSetSize() == check_spo_size);

  ions_.update();
  elec_.update();

  auto& twf(*twf_ptr);
  twf.setMassTerm(elec_);
  twf.evaluateLog(elec_);

  //Reference values from QWalk with SOC

  app_log() << "twf.evaluateLog logpsi " << std::setprecision(16) << twf.getLogPsi() << " " << twf.getPhase()
            << std::endl;
  CHECK(std::complex<double>(twf.getLogPsi(), twf.getPhase()) ==
        LogComplexApprox(std::complex<double>(-9.653087, 3.311467)));

  twf.prepareGroup(elec_, 0);
  ParticleSet::ComplexType spingrad_old;
  auto grad_old = twf.evalGradWithSpin(elec_, 1, spingrad_old);
  app_log() << "twf.evalGrad grad_old " << std::setprecision(16) << grad_old << std::endl;
  CHECK(grad_old[0] == ComplexApprox(ValueType(0.060932, -0.285244)).epsilon(1e-4));
  CHECK(grad_old[1] == ComplexApprox(ValueType(-0.401769, 0.180544)).epsilon(1e-4));
  CHECK(grad_old[2] == ComplexApprox(ValueType(0.174010, 0.140642)).epsilon(1e-4));
  CHECK(spingrad_old == ComplexApprox(ValueType(0.6766137, -0.8366186)).epsilon(1e-4));

  PosType delta(0.464586, 0.75017, 1.184383);
  double ds = 0.12;
  elec_.makeMoveWithSpin(0, delta, ds);

  ParticleSet::GradType grad_new;
  ParticleSet::ComplexType spingrad_new;
  auto ratio = twf.calcRatioGradWithSpin(elec_, 0, grad_new, spingrad_new);
  app_log() << "twf.calcRatioGrad ratio " << ratio << " grad_new " << grad_new << std::endl;
  CHECK(ValueType(std::abs(ratio)) == ValueApprox(0.991503).epsilon(1e-4));
  CHECK(grad_new[0] == ComplexApprox(ValueType(-0.631184, -0.136918)).epsilon(1e-4));
  CHECK(grad_new[1] == ComplexApprox(ValueType(0.074214, -0.080204)).epsilon(1e-4));
  CHECK(grad_new[2] == ComplexApprox(ValueType(-0.073180, -0.133539)).epsilon(1e-4));
  CHECK(spingrad_new == ComplexApprox(ValueType(-0.135438, -0.6085006)).epsilon(1e-4));

  ratio = twf.calcRatio(elec_, 0);
  app_log() << "twf.calcRatio ratio " << ratio << std::endl;
  CHECK(ValueType(std::abs(ratio)) == ValueApprox(0.991503).epsilon(1e-4));

  elec_.accept_rejectMove(0, false);

  //now lets test batched interface
  const int num_walkers = 2;
  ResourceCollection pset_res("test_pset_res");
  elec_.createResource(pset_res);
  ParticleSet elec_clone(elec_);
  RefVectorWithLeader<ParticleSet> p_list(elec_, {elec_, elec_clone});
  ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(pset_res, p_list);

  ResourceCollection twf_res("test_twf_res");
  twf.createResource(twf_res);
  auto twf_clone = twf.makeClone(elec_clone);
  RefVectorWithLeader<TrialWaveFunction> twf_list(twf, {twf, *twf_clone});
  ResourceCollectionTeamLock<TrialWaveFunction> mw_twf_lock(twf_res, twf_list);

  ParticleSet::mw_update(p_list);
  TrialWaveFunction::mw_evaluateLog(twf_list, p_list);

  for (int iw = 0; iw < num_walkers; iw++)
  {
    CHECK(std::complex<double>(twf_list[iw].getLogPsi(), twf_list[iw].getPhase()) ==
          LogComplexApprox(std::complex<double>(-9.653087, 3.311467)));
  }

  TrialWaveFunction::mw_prepareGroup(twf_list, p_list, 0);

  int moved_elec_id = 1;
  TWFGrads<CoordsType::POS_SPIN> grads_old(num_walkers);
  TrialWaveFunction::mw_evalGrad(twf_list, p_list, moved_elec_id, grads_old);
  for (int iw = 0; iw < num_walkers; iw++)
  {
    CHECK(grads_old.grads_positions[iw][0] == ComplexApprox(ValueType(0.060932, -0.285244)).epsilon(1e-4));
    CHECK(grads_old.grads_positions[iw][1] == ComplexApprox(ValueType(-0.401769, 0.180544)).epsilon(1e-4));
    CHECK(grads_old.grads_positions[iw][2] == ComplexApprox(ValueType(0.174010, 0.140642)).epsilon(1e-4));
    CHECK(grads_old.grads_spins[iw] == ComplexApprox(ValueType(0.6766137, -0.8366186)).epsilon(1e-4));
  }

  moved_elec_id = 0;
  MCCoords<CoordsType::POS_SPIN> displs(num_walkers);
  displs.positions = {delta, delta};
  displs.spins     = {ds, ds};
  ParticleSet::mw_makeMove(p_list, moved_elec_id, displs);

  std::vector<PsiValueType> ratios(num_walkers);
  TrialWaveFunction::mw_calcRatio(twf_list, p_list, moved_elec_id, ratios);
  for (int iw = 0; iw < num_walkers; iw++)
    CHECK(ValueType(std::abs(ratios[iw])) == ValueApprox(0.991503).epsilon(1e-4));
  std::fill(ratios.begin(), ratios.end(), ValueType(0));

  TWFGrads<CoordsType::POS_SPIN> grads_new(num_walkers);
  TrialWaveFunction::mw_calcRatioGrad(twf_list, p_list, moved_elec_id, ratios, grads_new);
  for (int iw = 0; iw < num_walkers; iw++)
  {
    CHECK(ValueType(std::abs(ratios[iw])) == ValueApprox(0.991503).epsilon(1e-4));
    CHECK(grads_new.grads_positions[iw][0] == ComplexApprox(ValueType(-0.631184, -0.136918)).epsilon(1e-4));
    CHECK(grads_new.grads_positions[iw][1] == ComplexApprox(ValueType(0.074214, -0.080204)).epsilon(1e-4));
    CHECK(grads_new.grads_positions[iw][2] == ComplexApprox(ValueType(-0.073180, -0.133539)).epsilon(1e-4));
    CHECK(grads_new.grads_spins[iw] == ComplexApprox(ValueType(-0.135438, -0.6085006)).epsilon(1e-4));
  }
}

TEST_CASE("Bi-spinor multi Slater dets", "[wavefunction]")
{
  app_log() << "-----------------------------------------------------------------" << std::endl;
  app_log() << "Bi using the table method no precomputation" << std::endl;
  app_log() << "-----------------------------------------------------------------" << std::endl;
  const char* spo_xml_string1 = R"(<wavefunction name="psi0" target="e">
    <sposet_builder name="spinorbuilder" type="molecularorbital" source="ion0" transform="yes" href="Bi.orbs.h5" precision="double">
        <sposet name="myspo" size="16">
            <occupation mode="ground"/>
        </sposet>
    </sposet_builder>
    <determinantset>
        <multideterminant optimize="no" spo_0="myspo" algorithm="table_method">
            <detlist size="4" type="DETS" nc0="0" ne0="5" nstates="16" cutoff="1e-20">
               <ci coeff=" 0.8586" occ0="1110110000000000"/>
               <ci coeff="-0.2040" occ0="1101110000000000"/>
               <ci coeff=" 0.4081" occ0="1110101000000000"/>
               <ci coeff="-0.2340" occ0="1101101000000000"/>
            </detlist>
        </multideterminant>
    </determinantset>
</wavefunction>)";
  test_Bi_msd(spo_xml_string1, "myspo", 16, 123);

  app_log() << "-----------------------------------------------------------------" << std::endl;
  app_log() << "Bi using the table method with new optimization" << std::endl;
  app_log() << "-----------------------------------------------------------------" << std::endl;
  const char* spo_xml_string1_new = R"(<wavefunction name="psi0" target="e">
    <sposet_builder name="spinorbuilder" type="molecularorbital" source="ion0" transform="yes" href="Bi.orbs.h5" precision="double">
        <sposet name="myspo" size="16">
            <occupation mode="ground"/>
        </sposet>
    </sposet_builder>
    <determinantset>
        <multideterminant optimize="no" spo_0="myspo" algorithm="precomputed_table_method">
            <detlist size="4" type="DETS" nc0="0" ne0="5" nstates="16" cutoff="1e-20">
               <ci coeff=" 0.8586" occ0="1110110000000000"/>
               <ci coeff="-0.2040" occ0="1101110000000000"/>
               <ci coeff=" 0.4081" occ0="1110101000000000"/>
               <ci coeff="-0.2340" occ0="1101101000000000"/>
            </detlist>
        </multideterminant>
    </determinantset>
</wavefunction>)";
  test_Bi_msd(spo_xml_string1_new, "myspo", 16, 123);

  app_log() << "-----------------------------------------------------------------" << std::endl;
  app_log() << "Bi using the table method with new optimization, read from hdf5" << std::endl;
  app_log() << "-----------------------------------------------------------------" << std::endl;
  const char* spo_xml_string2_new = R"(<wavefunction name="psi0" target="e">
    <sposet_builder name="spinorbuilder" type="molecularorbital" source="ion0" transform="yes" href="Bi.orbs.h5" precision="double">
        <sposet name="myspo" size="16">
            <occupation mode="ground"/>
        </sposet>
    </sposet_builder>
    <determinantset>
        <multideterminant optimize="no" spo_0="myspo" algorithm="precomputed_table_method">
            <detlist size="4" type="DETS" nc0="0" ne0="5" nstates="16" cutoff="1e-20" href="Bi.orbs.h5"/>
        </multideterminant>
    </determinantset>
</wavefunction>)";
  test_Bi_msd(spo_xml_string2_new, "myspo", 16, 123);
}
#endif
} // namespace qmcplusplus
