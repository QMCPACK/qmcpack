//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Jastrow/PolynomialFunctor3D.h"
#include "QMCWaveFunctions/Jastrow/JeeIOrbitalSoA.h"
#include "QMCWaveFunctions/Jastrow/eeI_JastrowBuilder.h"
#include "ParticleBase/ParticleAttribOps.h"
#include <ResourceCollection.h>
#include "QMCHamiltonians/NLPPJob.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
using RealType     = WaveFunctionComponent::RealType;
using LogValueType = WaveFunctionComponent::LogValueType;
using PsiValueType = WaveFunctionComponent::PsiValueType;

TEST_CASE("PolynomialFunctor3D functor zero", "[wavefunction]")
{
  PolynomialFunctor3D functor("test_functor");

  double r = 1.2;
  double u = functor.evaluate(r, r, r);
  REQUIRE(u == 0.0);
}

void test_J3_polynomial3D(const DynamicCoordinateKind kind_selected)
{
  Communicate* c = OHMMS::Controller;

  const SimulationCell simulation_cell;
  ParticleSet ions_(simulation_cell, kind_selected);
  ParticleSet elec_(simulation_cell, kind_selected);

  ions_.setName("ion");
  ions_.create({2});
  ions_.R[0] = {2.0, 0.0, 0.0};
  ions_.R[1] = {-2.0, 0.0, 0.0};
  SpeciesSet& source_species(ions_.getSpeciesSet());
  source_species.addSpecies("O");
  ions_.update();

  elec_.setName("elec");
  elec_.create({2, 2});
  elec_.R[0] = {1.00, 0.0, 0.0};
  elec_.R[1] = {0.0, 0.0, 0.0};
  elec_.R[2] = {-1.00, 0.0, 0.0};
  elec_.R[3] = {0.0, 0.0, 2.0};
  SpeciesSet& target_species(elec_.getSpeciesSet());
  int upIdx                          = target_species.addSpecies("u");
  int downIdx                        = target_species.addSpecies("d");
  int chargeIdx                      = target_species.addAttribute("charge");
  target_species(chargeIdx, upIdx)   = -1;
  target_species(chargeIdx, downIdx) = -1;
  //elec_.resetGroups();

  const char* particles = R"(<tmp>
    <jastrow name="J3" type="eeI" function="polynomial" source="ion" print="yes">
      <correlation ispecies="O" especies="u" isize="3" esize="3" rcut="10">
        <coefficients id="uuO" type="Array" optimize="yes"> 8.227710241e-06 2.480817653e-06 -5.354068112e-06 -1.112644787e-05 -2.208006078e-06 5.213121933e-06 -1.537865869e-05 8.899030233e-06 6.257255156e-06 3.214580988e-06 -7.716743107e-06 -5.275682077e-06 -1.778457637e-06 7.926231121e-06 1.767406868e-06 5.451359059e-08 2.801423724e-06 4.577282736e-06 7.634608083e-06 -9.510673173e-07 -2.344131575e-06 -1.878777219e-06 3.937363358e-07 5.065353773e-07 5.086724869e-07 -1.358768154e-07</coefficients>
      </correlation>
      <correlation ispecies="O" especies1="u" especies2="d" isize="3" esize="3" rcut="10">
        <coefficients id="udO" type="Array" optimize="yes"> -6.939530224e-06 2.634169299e-05 4.046077477e-05 -8.002682388e-06 -5.396795988e-06 6.697370507e-06 5.433953051e-05 -6.336849668e-06 3.680471431e-05 -2.996059772e-05 1.99365828e-06 -3.222705626e-05 -8.091669063e-06 4.15738535e-06 4.843939112e-06 3.563650208e-07 3.786332474e-05 -1.418336941e-05 2.282691374e-05 1.29239286e-06 -4.93580873e-06 -3.052539228e-06 9.870288001e-08 1.844286407e-06 2.970561871e-07 -4.364303677e-08</coefficients>
      </correlation>
    </jastrow>
</tmp>
)";
  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr jas_eeI = xmlFirstElementChild(root);

  eeI_JastrowBuilder jastrow(c, elec_, ions_);
  std::unique_ptr<WaveFunctionComponent> jas(jastrow.buildComponent(jas_eeI));

  using J3Type              = JeeIOrbitalSoA<PolynomialFunctor3D>;
  auto j3_uptr              = jastrow.buildComponent(jas_eeI);
  WaveFunctionComponent* j3 = dynamic_cast<J3Type*>(j3_uptr.get());
  REQUIRE(j3 != nullptr);

  // update all distance tables
  elec_.update();

  double logpsi_real = std::real(j3->evaluateLog(elec_, elec_.G, elec_.L));
  CHECK(logpsi_real == Approx(-1.193457749)); // note: number not validated

  double KE = -0.5 * (Dot(elec_.G, elec_.G) + Sum(elec_.L));
  CHECK(KE == Approx(-0.058051245)); // note: number not validated

  using ValueType = QMCTraits::ValueType;
  using PosType   = QMCTraits::PosType;

  // set virtutal particle position
  PosType newpos(0.3, 0.2, 0.5);

  elec_.makeVirtualMoves(newpos);
  std::vector<ValueType> ratios(elec_.getTotalNum());
  j3->evaluateRatiosAlltoOne(elec_, ratios);

  CHECK(std::real(ratios[0]) == Approx(0.8744938582));
  CHECK(std::real(ratios[1]) == Approx(1.0357541137));
  CHECK(std::real(ratios[2]) == Approx(0.8302245609));
  CHECK(std::real(ratios[3]) == Approx(0.7987703724));

  elec_.makeMove(0, newpos - elec_.R[0]);
  PsiValueType ratio_0 = j3->ratio(elec_, 0);
  elec_.rejectMove(0);

  elec_.makeMove(1, newpos - elec_.R[1]);
  PsiValueType ratio_1 = j3->ratio(elec_, 1);
  elec_.rejectMove(1);

  elec_.makeMove(2, newpos - elec_.R[2]);
  PsiValueType ratio_2 = j3->ratio(elec_, 2);
  elec_.rejectMove(2);

  elec_.makeMove(3, newpos - elec_.R[3]);
  PsiValueType ratio_3 = j3->ratio(elec_, 3);
  elec_.rejectMove(3);

  CHECK(std::real(ratio_0) == Approx(0.8744938582));
  CHECK(std::real(ratio_1) == Approx(1.0357541137));
  CHECK(std::real(ratio_2) == Approx(0.8302245609));
  CHECK(std::real(ratio_3) == Approx(0.7987703724));

  UniqueOptObjRefs opt_obj_refs;
  j3->extractOptimizableObjectRefs(opt_obj_refs);
  REQUIRE(opt_obj_refs.size() == 2);

  opt_variables_type optvars;
  Vector<WaveFunctionComponent::ValueType> dlogpsi;
  Vector<WaveFunctionComponent::ValueType> dhpsioverpsi;

  for (OptimizableObject& obj : opt_obj_refs)
    obj.checkInVariablesExclusive(optvars);
  optvars.resetIndex();
  const int NumOptimizables(optvars.size());
  j3->checkOutVariables(optvars);
  dlogpsi.resize(NumOptimizables);
  dhpsioverpsi.resize(NumOptimizables);
  j3->evaluateDerivatives(elec_, optvars, dlogpsi, dhpsioverpsi);

  app_log() << std::endl << "reporting dlogpsi and dhpsioverpsi" << std::scientific << std::endl;
  for (int iparam = 0; iparam < NumOptimizables; iparam++)
    app_log() << "param=" << iparam << " : " << dlogpsi[iparam] << "  " << dhpsioverpsi[iparam] << std::endl;
  app_log() << std::endl;

  CHECK(std::real(dlogpsi[43]) == Approx(1.3358726814e+05));
  CHECK(std::real(dhpsioverpsi[43]) == Approx(-2.3246270644e+05));

  Vector<WaveFunctionComponent::ValueType> dlogpsiWF;
  dlogpsiWF.resize(NumOptimizables);
  j3->evaluateDerivativesWF(elec_, optvars, dlogpsiWF);
  for (int i = 0; i < NumOptimizables; i++)
    CHECK(dlogpsi[i] == ValueApprox(dlogpsiWF[i]));

  VirtualParticleSet VP(elec_, 2);
  std::vector<PosType> newpos2(2);
  std::vector<ValueType> ratios2(2);
  newpos2[0] = newpos - elec_.R[1];
  newpos2[1] = PosType(0.2, 0.5, 0.3) - elec_.R[1];
  VP.makeMoves(elec_, 1, newpos2);
  j3->evaluateRatios(VP, ratios2);

  CHECK(std::real(ratios2[0]) == Approx(1.0357541137));
  CHECK(std::real(ratios2[1]) == Approx(1.0257141422));

  std::fill(ratios2.begin(), ratios2.end(), 0);
  Matrix<ValueType> dratio(2, NumOptimizables);
  j3->evaluateDerivRatios(VP, optvars, ratios2, dratio);

  CHECK(std::real(ratios2[0]) == Approx(1.0357541137));
  CHECK(std::real(ratios2[1]) == Approx(1.0257141422));
  CHECK(std::real(dratio[0][43]) == Approx(-1.4282569e+03));

  // testing batched interfaces
  ResourceCollection pset_res("test_pset_res");
  ResourceCollection wfc_res("test_wfc_res");

  elec_.createResource(pset_res);
  j3->createResource(wfc_res);

  // make a clones
  ParticleSet elec_clone(elec_);
  auto j3_clone = j3->makeClone(elec_clone);

  // testing batched interfaces
  RefVectorWithLeader<ParticleSet> p_ref_list(elec_, {elec_, elec_clone});
  RefVectorWithLeader<WaveFunctionComponent> j3_ref_list(*j3, {*j3, *j3_clone});

  ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(pset_res, p_ref_list);
  ResourceCollectionTeamLock<WaveFunctionComponent> mw_wfc_lock(wfc_res, j3_ref_list);

  std::vector<bool> isAccepted(2, true);
  ParticleSet::mw_update(p_ref_list);
  j3->mw_recompute(j3_ref_list, p_ref_list, isAccepted);

  // test NLPP related APIs
  const int nknot = 3;
  VirtualParticleSet vp(elec_, nknot), vp_clone(elec_clone, nknot);
  RefVectorWithLeader<VirtualParticleSet> vp_list(vp, {vp, vp_clone});
  ResourceCollection vp_res("test_vp_res");
  vp.createResource(vp_res);
  ResourceCollectionTeamLock<VirtualParticleSet> mw_vp_lock(vp_res, vp_list);

  const int ei_table_index = elec_.addTable(ions_);
  const auto& ei_table1    = elec_.getDistTableAB(ei_table_index);
  // make virtual move of elec 0, reference ion 1
  NLPPJob<RealType> job1(1, 0, ei_table1.getDistances()[0][1], -ei_table1.getDisplacements()[0][1]);
  const auto& ei_table2 = elec_clone.getDistTableAB(ei_table_index);
  // make virtual move of elec 1, reference ion 3
  NLPPJob<RealType> job2(3, 1, ei_table2.getDistances()[1][3], -ei_table2.getDisplacements()[1][3]);

  std::vector<PosType> deltaV1{{0.1, 0.2, 0.3}, {0.1, 0.3, 0.2}, {0.2, 0.1, 0.3}};
  std::vector<PosType> deltaV2{{0.02, 0.01, 0.03}, {0.02, 0.03, 0.01}, {0.03, 0.01, 0.02}};

  VirtualParticleSet::mw_makeMoves(vp_list, p_ref_list, {deltaV1, deltaV2}, {job1, job2}, false);

  std::vector<std::vector<ValueType>> nlpp_ratios(2);
  nlpp_ratios[0].resize(nknot);
  nlpp_ratios[1].resize(nknot);
  j3->mw_evaluateRatios(j3_ref_list, RefVectorWithLeader<const VirtualParticleSet>(vp, {vp, vp_clone}), nlpp_ratios);

  CHECK(ValueApprox(nlpp_ratios[0][0]) == ValueType(1.0273599625));
  CHECK(ValueApprox(nlpp_ratios[0][1]) == ValueType(1.0227555037));
  CHECK(ValueApprox(nlpp_ratios[0][2]) == ValueType(1.0473958254));
  CHECK(ValueApprox(nlpp_ratios[1][0]) == ValueType(1.0013145208));
  CHECK(ValueApprox(nlpp_ratios[1][1]) == ValueType(1.0011137724));
  CHECK(ValueApprox(nlpp_ratios[1][2]) == ValueType(1.0017225742));
}

TEST_CASE("PolynomialFunctor3D Jastrow", "[wavefunction]")
{
  test_J3_polynomial3D(DynamicCoordinateKind::DC_POS);
  test_J3_polynomial3D(DynamicCoordinateKind::DC_POS_OFFLOAD);
}
} // namespace qmcplusplus
