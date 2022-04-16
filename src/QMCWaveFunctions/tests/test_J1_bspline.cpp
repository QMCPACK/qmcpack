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
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "QMCWaveFunctions/Jastrow/RadialJastrowBuilder.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "QMCWaveFunctions/Jastrow/J1OrbitalSoA.h"
#include <ResourceCollection.h>
#include "QMCHamiltonians/NLPPJob.h"

#include <cstdio>
#include <string>


// Uncomment to print information and values from the underlying functor
//#define PRINT_SPLINE_DATA

using std::string;

namespace qmcplusplus
{
using RealType     = WaveFunctionComponent::RealType;
using PsiValueType = WaveFunctionComponent::PsiValueType;

TEST_CASE("BSpline functor zero", "[wavefunction]")
{
  BsplineFunctor<double> bf;

  double r = 1.2;
  double u = bf.evaluate(r);
  REQUIRE(u == 0.0);
}

TEST_CASE("BSpline functor one", "[wavefunction]")
{
  BsplineFunctor<double> bf;

  bf.resize(1);

  double r = 1.2;
  double u = bf.evaluate(r);
  REQUIRE(u == 0.0);
}

void test_J1_spline(const DynamicCoordinateKind kind_selected)
{
  Communicate* c = OHMMS::Controller;

  const SimulationCell simulation_cell;
  ParticleSet ions_(simulation_cell, kind_selected);
  ParticleSet elec_(simulation_cell);

  ions_.setName("ion");
  ions_.create({1});
  ions_.R[0][0] = 2.0;
  ions_.R[0][1] = 0.0;
  ions_.R[0][2] = 0.0;

  SpeciesSet& ispecies       = ions_.getSpeciesSet();
  int CIdx                   = ispecies.addSpecies("C");
  int ichargeIdx             = ispecies.addAttribute("charge");
  ispecies(ichargeIdx, CIdx) = 4;
  ions_.resetGroups();
  ions_.update();

  elec_.setName("elec");
  elec_.create({1, 1});
  elec_.R[0][0] = 1.00;
  elec_.R[0][1] = 0.0;
  elec_.R[0][2] = 0.0;
  elec_.R[1][0] = 0.0;
  elec_.R[1][1] = 0.0;
  elec_.R[1][2] = 0.0;

  SpeciesSet& tspecies         = elec_.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;
  elec_.resetGroups();

  const char* j1_xml_char = R"XML(<tmp>
   <jastrow type="One-Body" name="J1" function="bspline" source="ion" print="yes" gpu="no">
       <correlation elementType="C" rcut="10" size="8" cusp="0.0">
               <coefficients id="eC" type="Array">
-0.2032153051 -0.1625595974 -0.143124599 -0.1216434956 -0.09919771951 -0.07111729038 -0.04445345869 -0.02135082917
               </coefficients>
            </correlation>
         </jastrow>
</tmp>)XML";

  const char* j1_omptarget_xml_char = R"XML(<tmp>
   <jastrow type="One-Body" name="J1" function="bspline" source="ion" print="yes" gpu="omptarget">
       <correlation elementType="C" rcut="10" size="8" cusp="0.0">
               <coefficients id="eC" type="Array">
-0.2032153051 -0.1625595974 -0.143124599 -0.1216434956 -0.09919771951 -0.07111729038 -0.04445345869 -0.02135082917
               </coefficients>
            </correlation>
         </jastrow>
</tmp>)XML";

  Libxml2Document doc;
  bool okay =
      doc.parseFromString(kind_selected == DynamicCoordinateKind::DC_POS_OFFLOAD ? j1_omptarget_xml_char : j1_xml_char);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr jas1 = xmlFirstElementChild(root);

  RadialJastrowBuilder jastrow(c, elec_, ions_);

  using J1Type = J1OrbitalSoA<BsplineFunctor<RealType>>;
  auto j1_uptr = jastrow.buildComponent(jas1);
  J1Type* j1   = dynamic_cast<J1Type*>(j1_uptr.get());
  REQUIRE(j1);

  // update all distance tables
  elec_.update();

  double logpsi_real = std::real(j1->evaluateLog(elec_, elec_.G, elec_.L));
  REQUIRE(logpsi_real == Approx(0.3160552244)); // note: number not validated

  //Ionic Derivative Test.
  QMCTraits::GradType gsource(0.0);
  TinyVector<ParticleSet::ParticleGradient, OHMMS_DIM> grad_grad_source;
  TinyVector<ParticleSet::ParticleLaplacian, OHMMS_DIM> lapl_grad_source;
  int nelecs = elec_.getTotalNum();

  for (int dim = 0; dim < OHMMS_DIM; dim++)
  {
    grad_grad_source[dim].resize(nelecs);
    lapl_grad_source[dim].resize(nelecs);
  }

  /////////////////////////////////////////
  //Testing the ion gradient w.r.t. ion # 0
  /////////////////////////////////////////
  //NOTE:  All test values in this section are validated against finite differences
  //       for this configuration.

  //First we test evalGradSource(P,ions,ionid);

  gsource = j1->evalGradSource(elec_, ions_, 0);

  //Gradient comparison
  REQUIRE(std::real(gsource[0]) == Approx(-0.04695203659));
  REQUIRE(std::real(gsource[1]) == Approx(0.00000000000));
  REQUIRE(std::real(gsource[2]) == Approx(0.00000000000));

  //Now we test evalGradSource that returns higher order derivatives.
  gsource = j1->evalGradSource(elec_, ions_, 0, grad_grad_source, lapl_grad_source);

  //Gradient comparison
  REQUIRE(std::real(gsource[0]) == Approx(-0.04695203659));
  REQUIRE(std::real(gsource[1]) == Approx(0.00000000000));
  REQUIRE(std::real(gsource[2]) == Approx(0.00000000000));

  //Ion gradient of electron gradient comparison.
  REQUIRE(std::real(grad_grad_source[0][0][0]) == Approx(-0.008883672));
  REQUIRE(std::real(grad_grad_source[0][1][0]) == Approx(-0.002111879));
  REQUIRE(std::real(grad_grad_source[1][0][1]) == Approx(0.028489287));
  REQUIRE(std::real(grad_grad_source[1][1][1]) == Approx(0.009231375));
  REQUIRE(std::real(grad_grad_source[2][0][2]) == Approx(0.028489287));
  REQUIRE(std::real(grad_grad_source[2][1][2]) == Approx(0.009231375));

  //Ion gradient of electron laplacians.
  REQUIRE(std::real(lapl_grad_source[0][0]) == Approx(0.1494918378));
  REQUIRE(std::real(lapl_grad_source[0][1]) == Approx(-0.0056182539));
  REQUIRE(std::real(lapl_grad_source[1][0]) == Approx(0.0000000000));
  REQUIRE(std::real(lapl_grad_source[1][1]) == Approx(0.0000000000));
  REQUIRE(std::real(lapl_grad_source[2][0]) == Approx(0.0000000000));
  REQUIRE(std::real(lapl_grad_source[2][1]) == Approx(0.0000000000));


  // now test evaluateHessian
  WaveFunctionComponent::HessVector grad_grad_psi;
  grad_grad_psi.resize(elec_.getTotalNum());
  grad_grad_psi = 0.0;

  j1->evaluateHessian(elec_, grad_grad_psi);

  std::vector<double> hess_values = {
      0.00888367, 0, 0, 0, -0.0284893, 0, 0, 0, -0.0284893, 0.00211188, 0, 0, 0, -0.00923137, 0, 0, 0, -0.00923137,
  };

  int m = 0;
  for (int n = 0; n < elec_.getTotalNum(); n++)
    for (int i = 0; i < OHMMS_DIM; i++)
      for (int j = 0; j < OHMMS_DIM; j++, m++)
      {
        REQUIRE(std::real(grad_grad_psi[n](i, j)) == Approx(hess_values[m]));
      }

  j1->evaluateLog(elec_, elec_.G, elec_.L); // evaluateHessian has side effects


  struct JValues
  {
    double r;
    double u;
    double du;
    double ddu;
  };

  // Cut and paste from output of gen_bspline_jastrow.py
  const int N     = 20;
  JValues Vals[N] = {{0.00, -0.1896634025, 0, 0.06586224647},
                     {0.60, -0.1804990512, 0.02606308248, 0.02101469513},
                     {1.20, -0.1637586749, 0.0255799351, -0.01568108497},
                     {1.80, -0.1506226948, 0.01922435549, -0.005504180392},
                     {2.40, -0.1394848415, 0.01869442683, 0.001517191423},
                     {3.00, -0.128023472, 0.01946283614, 0.00104417293},
                     {3.60, -0.1161729491, 0.02009651096, 0.001689229059},
                     {4.20, -0.1036884223, 0.02172284322, 0.003731878464},
                     {4.80, -0.08992443283, 0.0240346508, 0.002736384838},
                     {5.40, -0.07519614609, 0.02475121662, -0.000347832122},
                     {6.00, -0.06054074137, 0.02397053075, -0.001842295859},
                     {6.60, -0.04654631918, 0.0225837382, -0.002780345968},
                     {7.20, -0.03347994129, 0.02104406699, -0.00218107833},
                     {7.80, -0.0211986378, 0.01996899618, -0.00173646255},
                     {8.40, -0.01004416026, 0.01635533409, -0.01030907776},
                     {9.00, -0.002594125744, 0.007782377232, -0.01556475446},
                     {9.60, -0.0001660240476, 0.001245180357, -0.006225901786},
                     {10.20, 0, 0, 0},
                     {10.80, 0, 0, 0},
                     {11.40, 0, 0, 0}};

  BsplineFunctor<RealType>* bf = j1->getFunctors()[0];

  for (int i = 0; i < N; i++)
  {
    RealType dv  = 0.0;
    RealType ddv = 0.0;
    RealType val = bf->evaluate(Vals[i].r, dv, ddv);
    REQUIRE(Vals[i].u == Approx(val));
    REQUIRE(Vals[i].du == Approx(dv));
    REQUIRE(Vals[i].ddu == Approx(ddv));
  }

#ifdef PRINT_SPLINE_DATA
  // write out values of the Bspline functor
  //BsplineFunctor<double> *bf = j1->J1Functors[0];
  printf("NumParams = %d\n", bf->NumParams);
  printf("CuspValue = %g\n", bf->CuspValue);
  printf("DeltaR = %g\n", bf->DeltaR);
  printf("SplineCoeffs size = %d\n", bf->SplineCoefs.size());
  for (int j = 0; j < bf->SplineCoefs.size(); j++)
  {
    printf("%d %g\n", j, bf->SplineCoefs[j]);
  }
  printf("\n");

  for (int i = 0; i < 20; i++)
  {
    double r      = 0.6 * i;
    elec_.R[0][0] = r;
    elec_.update();
    double logpsi_real = std::real(j1->evaluateLog(elec_, elec_.G, elec_.L));
    //double alt_val = bf->evaluate(r);
    double dv      = 0.0;
    double ddv     = 0.0;
    double alt_val = bf->evaluate(r, dv, ddv);
    printf("%g %g %g %g %g\n", r, logpsi_real, alt_val, dv, ddv);
  }
#endif

  using ValueType = QMCTraits::ValueType;
  using PosType   = QMCTraits::PosType;

  // set virtutal particle position
  PosType newpos(0.3, 0.2, 0.5);

  elec_.makeVirtualMoves(newpos);
  std::vector<ValueType> ratios(elec_.getTotalNum());
  j1->evaluateRatiosAlltoOne(elec_, ratios);

  REQUIRE(std::real(ratios[0]) == Approx(0.9819208747));
  REQUIRE(std::real(ratios[1]) == Approx(1.0040884258));

  elec_.makeMove(0, newpos - elec_.R[0]);
  PsiValueType ratio_0 = j1->ratio(elec_, 0);
  elec_.rejectMove(0);

  REQUIRE(std::real(ratio_0) == Approx(0.9819208747));

  // test acceptMove results
  elec_.makeMove(1, newpos - elec_.R[1]);
  PsiValueType ratio_1 = j1->ratio(elec_, 1);
  j1->acceptMove(elec_, 1);
  elec_.acceptMove(1);

  REQUIRE(std::real(ratio_1) == Approx(1.0040884258));
  REQUIRE(std::real(j1->get_log_value()) == Approx(0.32013531536));

  // test to make sure that setting cusp for J1 works properly
  const char* j1_xml_char_2 = R"XML(<tmp>
<jastrow type="One-Body" name="J1" function="bspline" source="ion" print="yes" gpu="no">
    <correlation elementType="C" rcut="10" size="8" cusp="2.0">
        <coefficients id="eC" type="Array">
 -0.2032153051 -0.1625595974 -0.143124599 -0.1216434956 -0.09919771951 -0.07111729038 -0.04445345869 -0.02135082917
        </coefficients>
    </correlation>
</jastrow>
</tmp>)XML";

  const char* j1_omptarget_xml_char_2 = R"XML(<tmp>
<jastrow type="One-Body" name="J1" function="bspline" source="ion" print="yes" gpu="omptarget">
    <correlation elementType="C" rcut="10" size="8" cusp="2.0">
        <coefficients id="eC" type="Array">
 -0.2032153051 -0.1625595974 -0.143124599 -0.1216434956 -0.09919771951 -0.07111729038 -0.04445345869 -0.02135082917
        </coefficients>
    </correlation>
</jastrow>
</tmp>)XML";

  Libxml2Document doc2;
  bool okay2 = doc2.parseFromString(kind_selected == DynamicCoordinateKind::DC_POS_OFFLOAD ? j1_omptarget_xml_char_2
                                                                                           : j1_xml_char_2);
  REQUIRE(okay2);

  xmlNodePtr root2 = doc2.getRoot();

  xmlNodePtr jas2 = xmlFirstElementChild(root2);

  RadialJastrowBuilder jastrow2(c, elec_, ions_);

  auto j12_uptr = jastrow2.buildComponent(jas2);
  J1Type* j12   = dynamic_cast<J1Type*>(j12_uptr.get());
  REQUIRE(j12);

  // Cut and paste from output of gen_bspline_jastrow.py
  // note only the first two rows should change from above
  const int N2      = 20;
  JValues Vals2[N2] = {{0.00, -0.9304041433, 2, -3.534137754},
                       {0.60, -0.252599792, 0.4492630825, -1.634985305},
                       {1.20, -0.1637586749, 0.0255799351, -0.01568108497},
                       {1.80, -0.1506226948, 0.01922435549, -0.005504180392},
                       {2.40, -0.1394848415, 0.01869442683, 0.001517191423},
                       {3.00, -0.128023472, 0.01946283614, 0.00104417293},
                       {3.60, -0.1161729491, 0.02009651096, 0.001689229059},
                       {4.20, -0.1036884223, 0.02172284322, 0.003731878464},
                       {4.80, -0.08992443283, 0.0240346508, 0.002736384838},
                       {5.40, -0.07519614609, 0.02475121662, -0.000347832122},
                       {6.00, -0.06054074137, 0.02397053075, -0.001842295859},
                       {6.60, -0.04654631918, 0.0225837382, -0.002780345968},
                       {7.20, -0.03347994129, 0.02104406699, -0.00218107833},
                       {7.80, -0.0211986378, 0.01996899618, -0.00173646255},
                       {8.40, -0.01004416026, 0.01635533409, -0.01030907776},
                       {9.00, -0.002594125744, 0.007782377232, -0.01556475446},
                       {9.60, -0.0001660240476, 0.001245180357, -0.006225901786},
                       {10.20, 0, 0, 0},
                       {10.80, 0, 0, 0},
                       {11.40, 0, 0, 0}};

  BsplineFunctor<RealType>* bf2 = j12->getFunctors()[0];

  for (int i = 0; i < N2; i++)
  {
    RealType dv  = 0.0;
    RealType ddv = 0.0;
    RealType val = bf2->evaluate(Vals2[i].r, dv, ddv);
    REQUIRE(Vals2[i].du == Approx(dv));
    REQUIRE(Vals2[i].u == Approx(val));
    REQUIRE(Vals2[i].ddu == Approx(ddv));
  }

  // testing batched interfaces
  ResourceCollection pset_res("test_pset_res");
  ResourceCollection wfc_res("test_wfc_res");

  elec_.createResource(pset_res);
  j12->createResource(wfc_res);

  // make a clones
  ParticleSet elec_clone(elec_);
  auto j12_clone = j12->makeClone(elec_clone);

  // testing batched interfaces
  RefVectorWithLeader<ParticleSet> p_ref_list(elec_, {elec_, elec_clone});
  RefVectorWithLeader<WaveFunctionComponent> j1_ref_list(*j12, {*j12, *j12_clone});

  ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(pset_res, p_ref_list);
  ResourceCollectionTeamLock<WaveFunctionComponent> mw_wfc_lock(wfc_res, j1_ref_list);

  std::vector<bool> isAccepted(2, true);
  ParticleSet::mw_update(p_ref_list);
  j12->mw_recompute(j1_ref_list, p_ref_list, isAccepted);

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
  NLPPJob<RealType> job1(1, 0, elec_.R[0], ei_table1.getDistances()[0][1], -ei_table1.getDisplacements()[0][1]);
  const auto& ei_table2 = elec_clone.getDistTableAB(ei_table_index);
  // make virtual move of elec 1, reference ion 3
  NLPPJob<RealType> job2(3, 1, elec_clone.R[1], ei_table2.getDistances()[1][3], -ei_table2.getDisplacements()[1][3]);

  std::vector<PosType> deltaV1{{0.1, 0.2, 0.3}, {0.1, 0.3, 0.2}, {0.2, 0.1, 0.3}};
  std::vector<PosType> deltaV2{{0.02, 0.01, 0.03}, {0.02, 0.03, 0.01}, {0.03, 0.01, 0.02}};

  VirtualParticleSet::mw_makeMoves(vp_list, {deltaV1, deltaV2}, {job1, job2}, false);

  std::vector<std::vector<ValueType>> nlpp_ratios(2);
  nlpp_ratios[0].resize(nknot);
  nlpp_ratios[1].resize(nknot);
  j12->mw_evaluateRatios(j1_ref_list, RefVectorWithLeader<const VirtualParticleSet>(vp, {vp, vp_clone}), nlpp_ratios);

  CHECK(ValueApprox(nlpp_ratios[0][0]) == ValueType(1.0016646504));
  CHECK(ValueApprox(nlpp_ratios[0][1]) == ValueType(1.0016646504));
  CHECK(ValueApprox(nlpp_ratios[0][2]) == ValueType(1.011890802));
  CHECK(ValueApprox(nlpp_ratios[1][0]) == ValueType(1.0001773019));
  CHECK(ValueApprox(nlpp_ratios[1][1]) == ValueType(1.000242931));
  CHECK(ValueApprox(nlpp_ratios[1][2]) == ValueType(1.0004189199));
}

TEST_CASE("BSpline builder Jastrow J1", "[wavefunction]")
{
  test_J1_spline(DynamicCoordinateKind::DC_POS);
  test_J1_spline(DynamicCoordinateKind::DC_POS_OFFLOAD);
}
} // namespace qmcplusplus
