//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminantBatched.h"
#include "QMCWaveFunctions/tests/FakeSPO.h"
#include "checkMatrix.hpp"
#include <ResourceCollection.h>

namespace qmcplusplus
{
using RealType     = QMCTraits::RealType;
using ValueType    = QMCTraits::ValueType;
using PosType      = QMCTraits::PosType;
using GradType     = QMCTraits::GradType;
using LogValueType = std::complex<QMCTraits::QTFull::RealType>;
using PsiValueType = QMCTraits::QTFull::ValueType;

template<class DET_ENGINE>
void test_DiracDeterminantBatched_first()
{
  using DetType  = DiracDeterminantBatched<DET_ENGINE>;
  auto spo_init  = std::make_unique<FakeSPO>();
  const int norb = 3;
  spo_init->setOrbitalSetSize(norb);
  DetType ddb(std::move(spo_init), 0, norb);
  auto spo = dynamic_cast<FakeSPO*>(ddb.getPhi());

  // occurs in call to registerData
  ddb.dpsiV.resize(norb);
  ddb.d2psiV.resize(norb);

  const SimulationCell simulation_cell;
  ParticleSet elec(simulation_cell);

  elec.create({3});
  ddb.recompute(elec);

  Matrix<ValueType> b;
  b.resize(3, 3);

  b(0, 0) = 0.6159749342;
  b(0, 1) = -0.2408954682;
  b(0, 2) = -0.1646081192;
  b(1, 0) = 0.07923894288;
  b(1, 1) = 0.1496231042;
  b(1, 2) = -0.1428117337;
  b(2, 0) = -0.2974298429;
  b(2, 1) = -0.04586322768;
  b(2, 2) = 0.3927890292;

  checkMatrix(ddb.get_det_engine().get_ref_psiMinv(), b);

  ParticleSet::GradType grad;
  PsiValueType det_ratio  = ddb.ratioGrad(elec, 0, grad);
  PsiValueType det_ratio1 = 0.178276269185;
  REQUIRE(det_ratio1 == ValueApprox(det_ratio));

  ddb.acceptMove(elec, 0);

  b(0, 0) = 3.455170657;
  b(0, 1) = -1.35124809;
  b(0, 2) = -0.9233316353;
  b(1, 0) = 0.05476311768;
  b(1, 1) = 0.1591951095;
  b(1, 2) = -0.1362710138;
  b(2, 0) = -2.235099338;
  b(2, 1) = 0.7119205298;
  b(2, 2) = 0.9105960265;

  checkMatrix(ddb.get_det_engine().get_ref_psiMinv(), b);

  // set virtutal particle position
  PosType newpos(0.3, 0.2, 0.5);

  elec.makeVirtualMoves(newpos);
  std::vector<ValueType> ratios(elec.getTotalNum());
  ddb.evaluateRatiosAlltoOne(elec, ratios);

  CHECK(std::real(ratios[0]) == Approx(1.2070809985));
  CHECK(std::real(ratios[1]) == Approx(0.2498726439));
  CHECK(std::real(ratios[2]) == Approx(-1.3145695364));

  elec.makeMove(0, newpos - elec.R[0]);
  PsiValueType ratio_0 = ddb.ratio(elec, 0);
  elec.rejectMove(0);

  CHECK(std::real(ratio_0) == Approx(-0.5343861437));

  VirtualParticleSet VP(elec, 2);
  std::vector<PosType> newpos2(2);
  std::vector<ValueType> ratios2(2);
  newpos2[0] = newpos - elec.R[1];
  newpos2[1] = PosType(0.2, 0.5, 0.3) - elec.R[1];
  VP.makeMoves(1, elec.R[1], newpos2);
  ddb.evaluateRatios(VP, ratios2);

  CHECK(std::real(ratios2[0]) == Approx(0.4880285278));
  CHECK(std::real(ratios2[1]) == Approx(0.9308456444));

  //test acceptMove
  elec.makeMove(1, newpos - elec.R[1]);
  PsiValueType ratio_1 = ddb.ratio(elec, 1);
  ddb.acceptMove(elec, 1);
  elec.acceptMove(1);

  CHECK(std::real(ratio_1) == Approx(0.9308456444));
  CHECK(std::real(ddb.get_log_value()) == Approx(1.9891064655));
}

TEST_CASE("DiracDeterminantBatched_first", "[wavefunction][fermion]")
{
#if defined(ENABLE_OFFLOAD) && defined(ENABLE_CUDA)
  test_DiracDeterminantBatched_first<MatrixDelayedUpdateCUDA<ValueType, QMCTraits::QTFull::ValueType>>();
#endif
  test_DiracDeterminantBatched_first<MatrixUpdateOMPTarget<ValueType, QMCTraits::QTFull::ValueType>>();
}

//#define DUMP_INFO

template<class DET_ENGINE>
void test_DiracDeterminantBatched_second()
{
  using DetType  = DiracDeterminantBatched<DET_ENGINE>;
  auto spo_init  = std::make_unique<FakeSPO>();
  const int norb = 4;
  spo_init->setOrbitalSetSize(norb);
  DetType ddb(std::move(spo_init), 0, norb);
  auto spo = dynamic_cast<FakeSPO*>(ddb.getPhi());

  // occurs in call to registerData
  ddb.dpsiV.resize(norb);
  ddb.d2psiV.resize(norb);

  const SimulationCell simulation_cell;
  ParticleSet elec(simulation_cell);

  elec.create({4});
  ddb.recompute(elec);

  Matrix<ValueType> orig_a;
  orig_a.resize(4, 4);
  orig_a = spo->a2;

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < norb; j++)
    {
      orig_a(j, i) = spo->v2(i, j);
    }
  }

  //check_matrix(ddb.getPsiMinv(), b);
  DiracMatrix<ValueType> dm;

  Matrix<ValueType> a_update1, scratchT;
  a_update1.resize(4, 4);
  scratchT.resize(4, 4);
  a_update1 = spo->a2;
  for (int j = 0; j < norb; j++)
  {
    a_update1(j, 0) = spo->v2(0, j);
  }

  Matrix<ValueType> a_update2;
  a_update2.resize(4, 4);
  a_update2 = spo->a2;
  for (int j = 0; j < norb; j++)
  {
    a_update2(j, 0) = spo->v2(0, j);
    a_update2(j, 1) = spo->v2(1, j);
  }

  Matrix<ValueType> a_update3;
  a_update3.resize(4, 4);
  a_update3 = spo->a2;
  for (int j = 0; j < norb; j++)
  {
    a_update3(j, 0) = spo->v2(0, j);
    a_update3(j, 1) = spo->v2(1, j);
    a_update3(j, 2) = spo->v2(2, j);
  }

  ParticleSet::GradType grad;
  PsiValueType det_ratio = ddb.ratioGrad(elec, 0, grad);

  simd::transpose(a_update1.data(), a_update1.rows(), a_update1.cols(), scratchT.data(), scratchT.rows(),
                  scratchT.cols());
  LogValueType det_update1;
  dm.invert_transpose(scratchT, a_update1, det_update1);
  PsiValueType det_ratio1 = LogToValue<ValueType>::convert(det_update1 - ddb.get_log_value());
#ifdef DUMP_INFO
  std::cout << "det 0 = " << std::exp(ddb.get_log_value()) << std::endl;
  std::cout << "det 1 = " << std::exp(det_update1) << std::endl;
  std::cout << "det ratio 1 = " << det_ratio1 << std::endl;
#endif
  //double det_ratio1 = 0.178276269185;

  REQUIRE(det_ratio1 == ValueApprox(det_ratio));

  ddb.acceptMove(elec, 0);

  PsiValueType det_ratio2 = ddb.ratioGrad(elec, 1, grad);
  LogValueType det_update2;
  simd::transpose(a_update2.data(), a_update2.rows(), a_update2.cols(), scratchT.data(), scratchT.rows(),
                  scratchT.cols());
  dm.invert_transpose(scratchT, a_update2, det_update2);
  PsiValueType det_ratio2_val = LogToValue<ValueType>::convert(det_update2 - det_update1);
#ifdef DUMP_INFO
  std::cout << "det 1 = " << std::exp(ddb.get_log_value()) << std::endl;
  std::cout << "det 2 = " << std::exp(det_update2) << std::endl;
  std::cout << "det ratio 2 = " << det_ratio2 << std::endl;
#endif
  //double det_ratio2_val = 0.178276269185;
  REQUIRE(det_ratio2 == ValueApprox(det_ratio2_val));

  ddb.acceptMove(elec, 1);

  PsiValueType det_ratio3 = ddb.ratioGrad(elec, 2, grad);
  LogValueType det_update3;
  simd::transpose(a_update3.data(), a_update3.rows(), a_update3.cols(), scratchT.data(), scratchT.rows(),
                  scratchT.cols());
  dm.invert_transpose(scratchT, a_update3, det_update3);
  PsiValueType det_ratio3_val = LogToValue<ValueType>::convert(det_update3 - det_update2);
#ifdef DUMP_INFO
  std::cout << "det 2 = " << std::exp(ddb.get_log_value()) << std::endl;
  std::cout << "det 3 = " << std::exp(det_update3) << std::endl;
  std::cout << "det ratio 3 = " << det_ratio3 << std::endl;
#endif
  REQUIRE(det_ratio3 == ValueApprox(det_ratio3_val));
  //check_value(det_ratio3, det_ratio3_val);

  ddb.acceptMove(elec, 2);

  simd::transpose(orig_a.data(), orig_a.rows(), orig_a.cols(), scratchT.data(), scratchT.rows(), scratchT.cols());
  dm.invert_transpose(scratchT, orig_a, det_update3);

#ifdef DUMP_INFO
  std::cout << "original " << std::endl;
  std::cout << orig_a << std::endl;
  std::cout << "block update " << std::endl;
  std::cout << ddb.getPsiMinv() << std::endl;
#endif

  checkMatrix(ddb.get_det_engine().get_ref_psiMinv(), orig_a);
}

TEST_CASE("DiracDeterminantBatched_second", "[wavefunction][fermion]")
{
#if defined(ENABLE_OFFLOAD) && defined(ENABLE_CUDA)
  test_DiracDeterminantBatched_second<MatrixDelayedUpdateCUDA<ValueType, QMCTraits::QTFull::ValueType>>();
#endif
  test_DiracDeterminantBatched_second<MatrixUpdateOMPTarget<ValueType, QMCTraits::QTFull::ValueType>>();
}

template<class DET_ENGINE>
void test_DiracDeterminantBatched_delayed_update(int delay_rank, DetMatInvertor matrix_inverter_kind)
{
  using DetType  = DiracDeterminantBatched<DET_ENGINE>;
  auto spo_init  = std::make_unique<FakeSPO>();
  const int norb = 4;
  spo_init->setOrbitalSetSize(norb);
  DetType ddc(std::move(spo_init), 0, norb, delay_rank, matrix_inverter_kind);
  auto spo = dynamic_cast<FakeSPO*>(ddc.getPhi());

  // occurs in call to registerData
  ddc.dpsiV.resize(norb);
  ddc.d2psiV.resize(norb);

  const SimulationCell simulation_cell;
  ParticleSet elec(simulation_cell);

  elec.create({4});
  ddc.recompute(elec);

  Matrix<ValueType> orig_a;
  orig_a.resize(4, 4);
  orig_a = spo->a2;

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < norb; j++)
    {
      orig_a(j, i) = spo->v2(i, j);
    }
  }

  //check_matrix(ddc.getPsiMinv(), b);
  DiracMatrix<ValueType> dm;

  Matrix<ValueType> a_update1, scratchT;
  scratchT.resize(4, 4);
  a_update1.resize(4, 4);
  a_update1 = spo->a2;
  for (int j = 0; j < norb; j++)
  {
    a_update1(j, 0) = spo->v2(0, j);
  }

  Matrix<ValueType> a_update2;
  a_update2.resize(4, 4);
  a_update2 = spo->a2;
  for (int j = 0; j < norb; j++)
  {
    a_update2(j, 0) = spo->v2(0, j);
    a_update2(j, 1) = spo->v2(1, j);
  }

  Matrix<ValueType> a_update3;
  a_update3.resize(4, 4);
  a_update3 = spo->a2;
  for (int j = 0; j < norb; j++)
  {
    a_update3(j, 0) = spo->v2(0, j);
    a_update3(j, 1) = spo->v2(1, j);
    a_update3(j, 2) = spo->v2(2, j);
  }


  ParticleSet::GradType grad;
  PsiValueType det_ratio = ddc.ratioGrad(elec, 0, grad);

  simd::transpose(a_update1.data(), a_update1.rows(), a_update1.cols(), scratchT.data(), scratchT.rows(),
                  scratchT.cols());
  LogValueType det_update1;
  dm.invert_transpose(scratchT, a_update1, det_update1);
  PsiValueType det_ratio1 = LogToValue<ValueType>::convert(det_update1 - ddc.get_log_value());
#ifdef DUMP_INFO
  std::cout << "det 0 = " << std::exp(ddc.get_log_value()) << std::endl;
  std::cout << "det 1 = " << std::exp(det_update1) << std::endl;
  std::cout << "det ratio 1 = " << det_ratio1 << std::endl;
#endif
  //double det_ratio1 = 0.178276269185;

  REQUIRE(det_ratio1 == ValueApprox(det_ratio));

  // update of Ainv in ddc is delayed
  ddc.acceptMove(elec, 0, true);
  // force update Ainv in ddc using SM-1 code path
  ddc.completeUpdates();

  checkMatrix(ddc.get_det_engine().get_ref_psiMinv(), a_update1);

  grad = ddc.evalGrad(elec, 1);

  PsiValueType det_ratio2 = ddc.ratioGrad(elec, 1, grad);
  simd::transpose(a_update2.data(), a_update2.rows(), a_update2.cols(), scratchT.data(), scratchT.rows(),
                  scratchT.cols());
  LogValueType det_update2;
  dm.invert_transpose(scratchT, a_update2, det_update2);
  PsiValueType det_ratio2_val = LogToValue<ValueType>::convert(det_update2 - det_update1);
#ifdef DUMP_INFO
  std::cout << "det 1 = " << std::exp(ddc.get_log_value()) << std::endl;
  std::cout << "det 2 = " << std::exp(det_update2) << std::endl;
  std::cout << "det ratio 2 = " << det_ratio2 << std::endl;
#endif
  // check ratio computed directly and the one computed by ddc with no delay
  //double det_ratio2_val = 0.178276269185;
  REQUIRE(det_ratio2 == ValueApprox(det_ratio2_val));

  // update of Ainv in ddc is delayed
  ddc.acceptMove(elec, 1, true);

  grad = ddc.evalGrad(elec, 2);

  PsiValueType det_ratio3 = ddc.ratioGrad(elec, 2, grad);
  simd::transpose(a_update3.data(), a_update3.rows(), a_update3.cols(), scratchT.data(), scratchT.rows(),
                  scratchT.cols());
  LogValueType det_update3;
  dm.invert_transpose(scratchT, a_update3, det_update3);
  PsiValueType det_ratio3_val = LogToValue<ValueType>::convert(det_update3 - det_update2);
#ifdef DUMP_INFO
  std::cout << "det 2 = " << std::exp(ddc.get_log_value()) << std::endl;
  std::cout << "det 3 = " << std::exp(det_update3) << std::endl;
  std::cout << "det ratio 3 = " << det_ratio3 << std::endl;
#endif
  // check ratio computed directly and the one computed by ddc with 1 delay
  REQUIRE(det_ratio3 == ValueApprox(det_ratio3_val));
  //check_value(det_ratio3, det_ratio3_val);

  // maximal delay reached and Ainv is updated fully
  ddc.acceptMove(elec, 2, true);
  ddc.completeUpdates();

  // fresh invert orig_a
  simd::transpose(orig_a.data(), orig_a.rows(), orig_a.cols(), scratchT.data(), scratchT.rows(), scratchT.cols());
  dm.invert_transpose(scratchT, orig_a, det_update3);

#ifdef DUMP_INFO
  std::cout << "original " << std::endl;
  std::cout << orig_a << std::endl;
  std::cout << "delayed update " << std::endl;
  std::cout << ddc.getPsiMinv() << std::endl;
#endif

  // compare all the elements of get_ref_psiMinv() in ddc and orig_a
  checkMatrix(ddc.get_det_engine().get_ref_psiMinv(), orig_a);

  // testing batched interfaces
  ResourceCollection pset_res("test_pset_res");
  ResourceCollection wfc_res("test_wfc_res");

  elec.createResource(pset_res);
  ddc.createResource(wfc_res);

  // make a clones
  ParticleSet elec_clone(elec);
  std::unique_ptr<WaveFunctionComponent> ddc_clone(ddc.makeCopy(ddc.getPhi()->makeClone()));
  auto& ddc_clone_ref = dynamic_cast<DetType&>(*ddc_clone);

  // testing batched interfaces
  RefVectorWithLeader<ParticleSet> p_ref_list(elec, {elec, elec_clone});
  RefVectorWithLeader<WaveFunctionComponent> ddc_ref_list(ddc, {ddc, *ddc_clone});

  ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(pset_res, p_ref_list);
  ResourceCollectionTeamLock<WaveFunctionComponent> mw_wfc_lock(wfc_res, ddc_ref_list);

  std::vector<bool> isAccepted(2, true);
  ParticleSet::mw_update(p_ref_list);
  ddc.mw_recompute(ddc_ref_list, p_ref_list, isAccepted);

  std::vector<PsiValueType> ratios(2);
  std::vector<GradType> grad_new(2);
  ddc.mw_ratioGrad(ddc_ref_list, p_ref_list, 0, ratios, grad_new);

  CHECK(det_ratio1 == ValueApprox(ratios[0]));
  CHECK(det_ratio1 == ValueApprox(ratios[1]));

  ddc.mw_accept_rejectMove(ddc_ref_list, p_ref_list, 0, isAccepted, true);
  ddc.mw_completeUpdates(ddc_ref_list);

  checkMatrix(ddc.get_det_engine().get_ref_psiMinv(), a_update1);
  checkMatrix(ddc_clone_ref.get_det_engine().get_ref_psiMinv(), a_update1);

  ddc.mw_evalGrad(ddc_ref_list, p_ref_list, 1, grad_new);
  ddc.mw_ratioGrad(ddc_ref_list, p_ref_list, 1, ratios, grad_new);

  CHECK(det_ratio2 == ValueApprox(ratios[0]));
  CHECK(det_ratio2 == ValueApprox(ratios[1]));

  ddc.mw_accept_rejectMove(ddc_ref_list, p_ref_list, 1, isAccepted, true);
  ddc.mw_evalGrad(ddc_ref_list, p_ref_list, 2, grad_new);
  ddc.mw_ratioGrad(ddc_ref_list, p_ref_list, 2, ratios, grad_new);

  CHECK(det_ratio3 == ValueApprox(ratios[0]));
  CHECK(det_ratio3 == ValueApprox(ratios[1]));

  ddc.mw_accept_rejectMove(ddc_ref_list, p_ref_list, 2, isAccepted, true);
  ddc.mw_completeUpdates(ddc_ref_list);

  checkMatrix(ddc.get_det_engine().get_ref_psiMinv(), orig_a);
  checkMatrix(ddc_clone_ref.get_det_engine().get_ref_psiMinv(), orig_a);
}

TEST_CASE("DiracDeterminantBatched_delayed_update", "[wavefunction][fermion]")
{
  // maximum delay 2
#if defined(ENABLE_OFFLOAD) && defined(ENABLE_CUDA)
  test_DiracDeterminantBatched_delayed_update<
      MatrixDelayedUpdateCUDA<ValueType, QMCTraits::QTFull::ValueType>>(2, DetMatInvertor::ACCEL);
  test_DiracDeterminantBatched_delayed_update<
      MatrixDelayedUpdateCUDA<ValueType, QMCTraits::QTFull::ValueType>>(2, DetMatInvertor::HOST);
#endif
  test_DiracDeterminantBatched_delayed_update<
      MatrixUpdateOMPTarget<ValueType, QMCTraits::QTFull::ValueType>>(2, DetMatInvertor::ACCEL);
  test_DiracDeterminantBatched_delayed_update<
      MatrixUpdateOMPTarget<ValueType, QMCTraits::QTFull::ValueType>>(2, DetMatInvertor::HOST);
}
} // namespace qmcplusplus
