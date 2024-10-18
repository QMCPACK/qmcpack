//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminant.h"
#include "QMCWaveFunctions/tests/FakeSPO.h"
#include "QMCWaveFunctions/ElectronGas/FreeOrbital.h"
#include "QMCWaveFunctions/SpinorSet.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
using RealType    = QMCTraits::RealType;
using ValueType   = QMCTraits::ValueType;
using ComplexType = QMCTraits::ComplexType;
using PosType     = QMCTraits::PosType;
using LogValue    = std::complex<QMCTraits::QTFull::RealType>;
using PsiValue    = QMCTraits::QTFull::ValueType;

template<typename T1, typename T2>
void check_matrix(Matrix<T1>& a, Matrix<T2>& b)
{
  REQUIRE(a.size() == b.size());
  for (int i = 0; i < a.rows(); i++)
  {
    for (int j = 0; j < a.cols(); j++)
    {
      CHECK(a(i, j) == ValueApprox(b(i, j)));
    }
  }
}

template<typename DET>
void test_DiracDeterminant_first(const DetMatInvertor inverter_kind)
{
  auto spo_init  = std::make_unique<FakeSPO>();
  const int norb = 3;
  spo_init->setOrbitalSetSize(norb);
  DET ddb(std::move(spo_init), 0, norb, 1, inverter_kind);
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

  check_matrix(ddb.psiM, b);

  ParticleSet::GradType grad;
  PsiValue det_ratio  = ddb.ratioGrad(elec, 0, grad);
  PsiValue det_ratio1 = 0.178276269185;
  CHECK(det_ratio1 == ValueApprox(det_ratio));

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

  check_matrix(ddb.psiM, b);

  // set virtutal particle position
  PosType newpos(0.3, 0.2, 0.5);

  elec.makeVirtualMoves(newpos);
  std::vector<ValueType> ratios(elec.getTotalNum());
  ddb.evaluateRatiosAlltoOne(elec, ratios);

  CHECK(std::real(ratios[0]) == Approx(1.2070809985));
  CHECK(std::real(ratios[1]) == Approx(0.2498726439));
  CHECK(std::real(ratios[2]) == Approx(-1.3145695364));

  elec.makeMove(0, newpos - elec.R[0]);
  PsiValue ratio_0 = ddb.ratio(elec, 0);
  elec.rejectMove(0);

  CHECK(std::real(ratio_0) == Approx(-0.5343861437));

  VirtualParticleSet VP(elec, 2);
  std::vector<PosType> newpos2(2);
  std::vector<ValueType> ratios2(2);
  newpos2[0] = newpos - elec.R[1];
  newpos2[1] = PosType(0.2, 0.5, 0.3) - elec.R[1];
  VP.makeMoves(elec, 1, newpos2);
  ddb.evaluateRatios(VP, ratios2);

  CHECK(std::real(ratios2[0]) == Approx(0.4880285278));
  CHECK(std::real(ratios2[1]) == Approx(0.9308456444));

  //test acceptMove
  elec.makeMove(1, newpos - elec.R[1]);
  PsiValue ratio_1 = ddb.ratio(elec, 1);
  ddb.acceptMove(elec, 1);
  elec.acceptMove(1);

  CHECK(std::real(ratio_1) == Approx(0.9308456444));
  CHECK(std::real(ddb.get_log_value()) == Approx(1.9891064655));
}

TEST_CASE("DiracDeterminant_first", "[wavefunction][fermion]")
{
  test_DiracDeterminant_first<DiracDeterminant<>>(DetMatInvertor::HOST);
  test_DiracDeterminant_first<DiracDeterminant<>>(DetMatInvertor::ACCEL);
#if defined(ENABLE_CUDA)
  test_DiracDeterminant_first<DiracDeterminant<DelayedUpdateCUDA<ValueType, QMCTraits::QTFull::ValueType>>>(
      DetMatInvertor::HOST);
  test_DiracDeterminant_first<DiracDeterminant<DelayedUpdateCUDA<ValueType, QMCTraits::QTFull::ValueType>>>(
      DetMatInvertor::ACCEL);
#elif defined(ENABLE_SYCL)
  test_DiracDeterminant_first<DiracDeterminant<DelayedUpdateSYCL<ValueType, QMCTraits::QTFull::ValueType>>>(
      DetMatInvertor::HOST);
#endif
}

//#define DUMP_INFO

template<typename DET>
void test_DiracDeterminant_second(const DetMatInvertor inverter_kind)
{
  auto spo_init  = std::make_unique<FakeSPO>();
  const int norb = 4;
  spo_init->setOrbitalSetSize(norb);
  DET ddb(std::move(spo_init), 0, norb, 1, inverter_kind);
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

  //check_matrix(ddb.psiM, b);
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
  PsiValue det_ratio = ddb.ratioGrad(elec, 0, grad);

  simd::transpose(a_update1.data(), a_update1.rows(), a_update1.cols(), scratchT.data(), scratchT.rows(),
                  scratchT.cols());
  LogValue det_update1;
  dm.invert_transpose(scratchT, a_update1, det_update1);
  PsiValue det_ratio1 = LogToValue<ValueType>::convert(det_update1 - ddb.get_log_value());
#ifdef DUMP_INFO
  app_log() << "det 0 = " << std::exp(ddb.get_log_value()) << std::endl;
  app_log() << "det 1 = " << std::exp(det_update1) << std::endl;
  app_log() << "det ratio 1 = " << det_ratio1 << std::endl;
#endif
  //double det_ratio1 = 0.178276269185;

  CHECK(det_ratio1 == ValueApprox(det_ratio));

  ddb.acceptMove(elec, 0);

  PsiValue det_ratio2 = ddb.ratioGrad(elec, 1, grad);
  LogValue det_update2;
  simd::transpose(a_update2.data(), a_update2.rows(), a_update2.cols(), scratchT.data(), scratchT.rows(),
                  scratchT.cols());
  dm.invert_transpose(scratchT, a_update2, det_update2);
  PsiValue det_ratio2_val = LogToValue<ValueType>::convert(det_update2 - det_update1);
#ifdef DUMP_INFO
  app_log() << "det 1 = " << std::exp(ddb.get_log_value()) << std::endl;
  app_log() << "det 2 = " << std::exp(det_update2) << std::endl;
  app_log() << "det ratio 2 = " << det_ratio2 << std::endl;
#endif
  //double det_ratio2_val = 0.178276269185;
  CHECK(det_ratio2 == ValueApprox(det_ratio2_val));

  ddb.acceptMove(elec, 1);

  PsiValue det_ratio3 = ddb.ratioGrad(elec, 2, grad);
  LogValue det_update3;
  simd::transpose(a_update3.data(), a_update3.rows(), a_update3.cols(), scratchT.data(), scratchT.rows(),
                  scratchT.cols());
  dm.invert_transpose(scratchT, a_update3, det_update3);
  PsiValue det_ratio3_val = LogToValue<ValueType>::convert(det_update3 - det_update2);
#ifdef DUMP_INFO
  app_log() << "det 2 = " << std::exp(ddb.get_log_value()) << std::endl;
  app_log() << "det 3 = " << std::exp(det_update3) << std::endl;
  app_log() << "det ratio 3 = " << det_ratio3 << std::endl;
#endif
  CHECK(det_ratio3 == ValueApprox(det_ratio3_val));
  //check_value(det_ratio3, det_ratio3_val);

  ddb.acceptMove(elec, 2);

  simd::transpose(orig_a.data(), orig_a.rows(), orig_a.cols(), scratchT.data(), scratchT.rows(), scratchT.cols());
  dm.invert_transpose(scratchT, orig_a, det_update3);

#ifdef DUMP_INFO
  app_log() << "original " << std::endl;
  app_log() << orig_a << std::endl;
  app_log() << "block update " << std::endl;
  app_log() << ddb.psiM << std::endl;
#endif

  check_matrix(orig_a, ddb.psiM);
}

TEST_CASE("DiracDeterminant_second", "[wavefunction][fermion]")
{
  test_DiracDeterminant_second<DiracDeterminant<>>(DetMatInvertor::HOST);
  test_DiracDeterminant_second<DiracDeterminant<>>(DetMatInvertor::ACCEL);
#if defined(ENABLE_CUDA)
  test_DiracDeterminant_second<DiracDeterminant<DelayedUpdateCUDA<ValueType, QMCTraits::QTFull::ValueType>>>(
      DetMatInvertor::HOST);
  test_DiracDeterminant_second<DiracDeterminant<DelayedUpdateCUDA<ValueType, QMCTraits::QTFull::ValueType>>>(
      DetMatInvertor::ACCEL);
#elif defined(ENABLE_SYCL)
  test_DiracDeterminant_second<DiracDeterminant<DelayedUpdateSYCL<ValueType, QMCTraits::QTFull::ValueType>>>(
      DetMatInvertor::HOST);
#endif
}

template<typename DET>
void test_DiracDeterminant_delayed_update(const DetMatInvertor inverter_kind)
{
  auto spo_init  = std::make_unique<FakeSPO>();
  const int norb = 4;
  spo_init->setOrbitalSetSize(norb);
  // maximum delay 2
  DET ddc(std::move(spo_init), 0, norb, 2, inverter_kind);
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

  //check_matrix(ddc.psiM, b);
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
  PsiValue det_ratio = ddc.ratioGrad(elec, 0, grad);

  simd::transpose(a_update1.data(), a_update1.rows(), a_update1.cols(), scratchT.data(), scratchT.rows(),
                  scratchT.cols());
  LogValue det_update1;
  dm.invert_transpose(scratchT, a_update1, det_update1);
  PsiValue det_ratio1 = LogToValue<ValueType>::convert(det_update1 - ddc.get_log_value());
#ifdef DUMP_INFO
  app_log() << "det 0 = " << std::exp(ddc.get_log_value()) << std::endl;
  app_log() << "det 1 = " << std::exp(det_update1) << std::endl;
  app_log() << "det ratio 1 = " << det_ratio1 << std::endl;
#endif
  //double det_ratio1 = 0.178276269185;

  CHECK(det_ratio1 == ValueApprox(det_ratio));

  // update of Ainv in ddc is delayed
  ddc.acceptMove(elec, 0, true);
  // force update Ainv in ddc using SM-1 code path
  // also force population of dpsiM or you will have invalid grad below
  ddc.completeUpdates();
  check_matrix(a_update1, ddc.psiM);

  grad = ddc.evalGrad(elec, 1);

  PsiValue det_ratio2 = ddc.ratioGrad(elec, 1, grad);
  simd::transpose(a_update2.data(), a_update2.rows(), a_update2.cols(), scratchT.data(), scratchT.rows(),
                  scratchT.cols());
  LogValue det_update2;
  dm.invert_transpose(scratchT, a_update2, det_update2);
  PsiValue det_ratio2_val = LogToValue<ValueType>::convert(det_update2 - det_update1);
#ifdef DUMP_INFO
  app_log() << "det 1 = " << std::exp(ddc.get_log_value()) << std::endl;
  app_log() << "det 2 = " << std::exp(det_update2) << std::endl;
  app_log() << "det ratio 2 = " << det_ratio2 << std::endl;
#endif
  // check ratio computed directly and the one computed by ddc with no delay
  //double det_ratio2_val = 0.178276269185;
  CHECK(det_ratio2 == ValueApprox(det_ratio2_val));

  // update of Ainv in ddc is delayed
  ddc.acceptMove(elec, 1, true);

  grad = ddc.evalGrad(elec, 2);

  PsiValue det_ratio3 = ddc.ratioGrad(elec, 2, grad);
  simd::transpose(a_update3.data(), a_update3.rows(), a_update3.cols(), scratchT.data(), scratchT.rows(),
                  scratchT.cols());
  LogValue det_update3;
  dm.invert_transpose(scratchT, a_update3, det_update3);
  PsiValue det_ratio3_val = LogToValue<ValueType>::convert(det_update3 - det_update2);
#ifdef DUMP_INFO
  app_log() << "det 2 = " << std::exp(ddc.get_log_value()) << std::endl;
  app_log() << "det 3 = " << std::exp(det_update3) << std::endl;
  app_log() << "det ratio 3 = " << det_ratio3 << std::endl;
#endif
  // check ratio computed directly and the one computed by ddc with 1 delay
  CHECK(det_ratio3 == ValueApprox(det_ratio3_val));
  //check_value(det_ratio3, det_ratio3_val);

  // maximal delay reached and Ainv is updated fully
  ddc.acceptMove(elec, 2, true);
  ddc.completeUpdates();

  // fresh invert orig_a
  simd::transpose(orig_a.data(), orig_a.rows(), orig_a.cols(), scratchT.data(), scratchT.rows(), scratchT.cols());
  dm.invert_transpose(scratchT, orig_a, det_update3);

#ifdef DUMP_INFO
  app_log() << "original " << std::endl;
  app_log() << orig_a << std::endl;
  app_log() << "delayed update " << std::endl;
  app_log() << ddc.psiM << std::endl;
#endif

  // compare all the elements of psiM in ddc and orig_a
  check_matrix(orig_a, ddc.psiM);
}

TEST_CASE("DiracDeterminant_delayed_update", "[wavefunction][fermion]")
{
  test_DiracDeterminant_delayed_update<DiracDeterminant<>>(DetMatInvertor::HOST);
  test_DiracDeterminant_delayed_update<DiracDeterminant<>>(DetMatInvertor::ACCEL);
#if defined(ENABLE_CUDA)
  test_DiracDeterminant_delayed_update<DiracDeterminant<DelayedUpdateCUDA<ValueType, QMCTraits::QTFull::ValueType>>>(
      DetMatInvertor::HOST);
  test_DiracDeterminant_delayed_update<DiracDeterminant<DelayedUpdateCUDA<ValueType, QMCTraits::QTFull::ValueType>>>(
      DetMatInvertor::ACCEL);
#elif defined(ENABLE_SYCL)
  test_DiracDeterminant_delayed_update<DiracDeterminant<DelayedUpdateSYCL<ValueType, QMCTraits::QTFull::ValueType>>>(
      DetMatInvertor::HOST);
#endif
}

#ifdef QMC_COMPLEX
template<typename DET>
void test_DiracDeterminant_spinor_update(const DetMatInvertor inverter_kind)
{
  using ValueType         = QMCTraits::ValueType;
  using PosType           = QMCTraits::PosType;
  using GradType          = QMCTraits::GradType;
  using LogValue          = WaveFunctionComponent::LogValue;
  using ParticlePos       = ParticleSet::ParticlePos;
  using ParticleGradient  = ParticleSet::ParticleGradient;
  using ParticleLaplacian = ParticleSet::ParticleLaplacian;

  // O2 test example from pwscf non-collinear calculation.
  ParticleSet::ParticleLayout lattice;
  lattice.R = {5.10509515, -3.23993545, 0.00000000, 5.10509515, 3.23993545,
               0.00000000, -6.49690625, 0.00000000, 7.08268015};

  //Shamelessly stealing this from test_einset.cpp.  3 particles though.
  const SimulationCell simulation_cell(lattice);
  ParticleSet ions_(simulation_cell);
  ParticleSet elec_(simulation_cell);
  ions_.setName("ion");
  ions_.create({2});

  ions_.R[0] = {0.00000000, 0.00000000, 1.08659253};
  ions_.R[1] = {0.00000000, 0.00000000, -1.08659253};
  elec_.setName("elec");
  elec_.create({3});
  elec_.R[0]     = {0.1, -0.3, 1.0};
  elec_.R[1]     = {-0.1, 0.3, 1.0};
  elec_.R[2]     = {0.1, 0.2, 0.3};
  elec_.spins[0] = 0.0;
  elec_.spins[1] = 0.2;
  elec_.spins[2] = 0.4;
  elec_.setSpinor(true);

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;

  elec_.addTable(ions_);
  elec_.resetGroups();
  elec_.update();
  // </steal>


  const auto nelec = elec_.R.size();
  //Our test case is going to be three electron gas orbitals distinguished by 3 different kpoints.
  //Independent SPO's for the up and down channels.
  //
  std::vector<PosType> kup, kdn;
  std::vector<RealType> k2up, k2dn;


  kup.resize(nelec);
  kup[0] = PosType(0, 0, 0);
  kup[1] = PosType(0.1, 0.2, 0.3);
  kup[2] = PosType(0.4, 0.5, 0.6);

  kdn.resize(nelec);
  kdn[0] = PosType(0, 0, 0);
  kdn[1] = PosType(-0.1, 0.2, -0.3);
  kdn[2] = PosType(0.4, -0.5, 0.6);

  auto spo_up = std::make_unique<FreeOrbital>("free_orb_up", kup);
  auto spo_dn = std::make_unique<FreeOrbital>("free_orb_up", kdn);

  auto spinor_set = std::make_unique<SpinorSet>("free_orb_spinor");
  spinor_set->set_spos(std::move(spo_up), std::move(spo_dn));

  DET dd(std::move(spinor_set), 0, nelec, 1, inverter_kind);
  app_log() << " nelec=" << nelec << std::endl;

  ParticleGradient G;
  ParticleLaplacian L;
  ParticleAttrib<ComplexType> SG;

  G.resize(nelec);
  L.resize(nelec);
  SG.resize(nelec);

  G  = 0.0;
  L  = 0.0;
  SG = 0.0;

  PosType dr(0.1, -0.05, 0.2);
  RealType ds = 0.3;

  app_log() << " BEFORE\n";
  app_log() << " R = " << elec_.R << std::endl;
  app_log() << " s = " << elec_.spins << std::endl;

  //In this section, we're going to test that values and various derivatives come out
  //correctly at the reference configuration.

  LogValue logref = dd.evaluateLog(elec_, G, L);

  CHECK(logref == ComplexApprox(ValueType(-1.1619939279564413, 0.8794794652468605)));
  CHECK(G[0][0] == ComplexApprox(ValueType(0.13416635, 0.2468612)));
  CHECK(G[0][1] == ComplexApprox(ValueType(-1.1165475, 0.71497753)));
  CHECK(G[0][2] == ComplexApprox(ValueType(0.0178403, 0.08212244)));
  CHECK(G[1][0] == ComplexApprox(ValueType(1.00240841, 0.12371593)));
  CHECK(G[1][1] == ComplexApprox(ValueType(1.62679698, -0.41080777)));
  CHECK(G[1][2] == ComplexApprox(ValueType(1.81324632, 0.78589013)));
  CHECK(G[2][0] == ComplexApprox(ValueType(-1.10994555, 0.15525902)));
  CHECK(G[2][1] == ComplexApprox(ValueType(-0.46335602, -0.50809713)));
  CHECK(G[2][2] == ComplexApprox(ValueType(-1.751199, 0.10949589)));
  CHECK(L[0] == ComplexApprox(ValueType(-2.06554158, 1.18145239)));
  CHECK(L[1] == ComplexApprox(ValueType(-5.06340536, 0.82126749)));
  CHECK(L[2] == ComplexApprox(ValueType(-4.82375261, -1.97943258)));

  //This is a workaround for the fact that I haven't implemented
  // evaluateLogWithSpin().  Shouldn't be needed unless we do drifted all-electron moves...
  for (int iat = 0; iat < nelec; iat++)
    dd.evalGradWithSpin(elec_, iat, SG[iat]);

  CHECK(SG[0] == ComplexApprox(ValueType(-1.05686704, -2.01802154)));
  CHECK(SG[1] == ComplexApprox(ValueType(1.18922259, 2.80414598)));
  CHECK(SG[2] == ComplexApprox(ValueType(-0.62617675, -0.51093984)));

  GradType g_singleeval(0.0);
  g_singleeval = dd.evalGrad(elec_, 1);

  CHECK(g_singleeval[0] == ComplexApprox(G[1][0]));
  CHECK(g_singleeval[1] == ComplexApprox(G[1][1]));
  CHECK(g_singleeval[2] == ComplexApprox(G[1][2]));


  //And now we're going to propose a trial spin+particle move and check the ratio and gradients at the
  //new location.
  //
  elec_.makeMoveAndCheckWithSpin(1, dr, ds);

  ValueType ratio_new;
  ValueType spingrad_new;
  GradType grad_new;

  //This tests ratio only evaluation.  Indirectly a call to evaluate(P,iat)
  ratio_new = dd.ratio(elec_, 1);
  CHECK(ratio_new == ComplexApprox(ValueType(1.7472917722050971, 1.1900872950904169)));

  ratio_new = dd.ratioGrad(elec_, 1, grad_new);
  CHECK(ratio_new == ComplexApprox(ValueType(1.7472917722050971, 1.1900872950904169)));
  CHECK(grad_new[0] == ComplexApprox(ValueType(0.5496675534224996, -0.07968022499097227)));
  CHECK(grad_new[1] == ComplexApprox(ValueType(0.4927399293808675, -0.29971549854643653)));
  CHECK(grad_new[2] == ComplexApprox(ValueType(1.2792642963632226, 0.12110307514989149)));

  grad_new     = 0;
  spingrad_new = 0;
  ratio_new    = dd.ratioGradWithSpin(elec_, 1, grad_new, spingrad_new);
  CHECK(ratio_new == ComplexApprox(ValueType(1.7472917722050971, 1.1900872950904169)));
  CHECK(grad_new[0] == ComplexApprox(ValueType(0.5496675534224996, -0.07968022499097227)));
  CHECK(grad_new[1] == ComplexApprox(ValueType(0.4927399293808675, -0.29971549854643653)));
  CHECK(grad_new[2] == ComplexApprox(ValueType(1.2792642963632226, 0.12110307514989149)));
  CHECK(spingrad_new == ComplexApprox(ValueType(1.164708841479661, 0.9576425115390172)));


  //Cool.  Now we test the transition between rejecting a move and accepting a move.
  //Reject the move first.  We want to see if everything stays the same.  evalGrad and evalSpinGrad for ease of use.

  elec_.rejectMove(1);
  //Going to check evalGrad and evalGradWithSpin for simplicity.
  g_singleeval = dd.evalGrad(elec_, 1);
  CHECK(g_singleeval[0] == ComplexApprox(G[1][0]));
  CHECK(g_singleeval[1] == ComplexApprox(G[1][1]));
  CHECK(g_singleeval[2] == ComplexApprox(G[1][2]));

  ValueType spingrad_old_test;
  g_singleeval = dd.evalGradWithSpin(elec_, 1, spingrad_old_test);

  CHECK(spingrad_old_test == ComplexApprox(SG[1]));
  CHECK(g_singleeval[0] == ComplexApprox(G[1][0]));
  CHECK(g_singleeval[1] == ComplexApprox(G[1][1]));
  CHECK(g_singleeval[2] == ComplexApprox(G[1][2]));

  //Now we test what happens if we accept a move...
  elec_.makeMoveAndCheckWithSpin(1, dr, ds);
  elec_.acceptMove(1);

  LogValue lognew(0.0);
  G      = 0.0; //evalauteLog += onto the G and L arguments.  So we zero them out.
  L      = 0.0;
  SG     = 0.0;
  lognew = dd.evaluateLog(elec_, G, L);

  for (int iat = 0; iat < nelec; iat++)
    dd.evalGradWithSpin(elec_, iat, SG[iat]);
  //logval for the new configuration has been computed with python.
  //The others reference values are computed earlier in this section.  New values equal the previous
  // "new values" associated with the previous trial moves.
  CHECK(lognew == ComplexApprox(ValueType(-0.41337396772929913, 1.4774106123071726)));
  CHECK(G[1][0] == ComplexApprox(grad_new[0]));
  CHECK(G[1][1] == ComplexApprox(grad_new[1]));
  CHECK(G[1][2] == ComplexApprox(grad_new[2]));
  CHECK(SG[1] == ComplexApprox(spingrad_new));
}

TEST_CASE("DiracDeterminant_spinor_update", "[wavefunction][fermion]")
{
  test_DiracDeterminant_spinor_update<DiracDeterminant<>>(DetMatInvertor::HOST);
  test_DiracDeterminant_spinor_update<DiracDeterminant<>>(DetMatInvertor::ACCEL);
#ifdef ENABLE_CUDA
  test_DiracDeterminant_spinor_update<DiracDeterminant<DelayedUpdateCUDA<ValueType, QMCTraits::QTFull::ValueType>>>(
      DetMatInvertor::HOST);
  test_DiracDeterminant_spinor_update<DiracDeterminant<DelayedUpdateCUDA<ValueType, QMCTraits::QTFull::ValueType>>>(
      DetMatInvertor::ACCEL);
#endif
}
#endif

} // namespace qmcplusplus
