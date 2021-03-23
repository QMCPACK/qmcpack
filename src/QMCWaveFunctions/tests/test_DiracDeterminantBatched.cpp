//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminantBatched.h"
#include "QMCWaveFunctions/tests/FakeSPO.h"
#include "QMCWaveFunctions/tests/CheckMatrix.hpp"

#ifdef QMC_COMPLEX //This is for the spinor test.
#include "QMCWaveFunctions/ElectronGas/ElectronGasComplexOrbitalBuilder.h"
#endif

#include "QMCWaveFunctions/SpinorSet.h"

#include <stdio.h>
#include <string>

using std::string;


namespace qmcplusplus
{
using RealType     = QMCTraits::RealType;
using ValueType    = QMCTraits::ValueType;
using LogValueType = std::complex<QMCTraits::QTFull::RealType>;
using PsiValueType = QMCTraits::QTFull::ValueType;

#if defined(ENABLE_OFFLOAD)
#if defined(ENABLE_CUDA)
typedef DiracDeterminantBatched<MatrixDelayedUpdateCUDA<ValueType, QMCTraits::QTFull::ValueType>> DetType;
#else
typedef DiracDeterminantBatched<MatrixUpdateOMPTarget<ValueType, QMCTraits::QTFull::ValueType>> DetType;
#endif
#else
typedef DiracDeterminantBatched<> DetType;
#endif

namespace testing
{
// checks if clas has a member handles
template<class T, typename = void>
struct has_handles : std::false_type
{};
// Other handle using det engines could be added.
template<class T>
struct has_handles<T, decltype(std::declval<T>().cuda_handles_, void())> : std::true_type
{};


class DiracDeterminantBatchedTest
{
public:
  template<class DET_ENGINE>
  DiracDeterminantBatchedMultiWalkerResource<DET_ENGINE>& get_mw_res(DiracDeterminantBatched<DET_ENGINE>& ddb)
  {
    return *(ddb.mw_res_);
  }
  template<class DET_ENGINE>
  void guardMultiWalkerRes(DiracDeterminantBatched<DET_ENGINE>& ddb)
  {
    return ddb.guardMultiWalkerRes();
  }

  template<class DET_ENGINE, typename = std::enable_if<testing::has_handles<DET_ENGINE>::value>>
  typename DET_ENGINE::Handles& getHandles(DiracDeterminantBatched<DET_ENGINE>& ddb)
  {
    return ddb.get_det_engine().getHandles();
  }
};
}

// I've observed downstream unit tests fail due to small issues in resource create/acquire/release.
// unit test these.
TEST_CASE("DiracDeterminantBatched_resources", "[wavefunction][fermion]")
{
  int count = 0;

  auto spo_init = std::make_unique<FakeSPO>();
  spo_init->setOrbitalSetSize(3);
  DetType ddb(std::move(spo_init));
  auto spo = dynamic_cast<FakeSPO*>(ddb.getPhi());

  ResourceCollection res_col("test resources");
  ddb.createResource(res_col);
  ddb.acquireResource(res_col);
  res_col.rewind();
  ddb.releaseResource(res_col);
}
  
TEST_CASE("DiracDeterminantBatched_first", "[wavefunction][fermion]")
{
  auto spo_init = std::make_unique<FakeSPO>();
  spo_init->setOrbitalSetSize(3);
  DetType ddb(std::move(spo_init));
  auto spo = dynamic_cast<FakeSPO*>(ddb.getPhi());

  int norb = 3;
  ddb.set(0, norb);

  // occurs in call to registerData
  ddb.dpsiV.resize(norb);
  ddb.d2psiV.resize(norb);


  ParticleSet elec;
  elec.create(3);

  testing::DiracDeterminantBatchedTest ddbt;
  ddbt.guardMultiWalkerRes(ddb);

  ResourceCollection res_col("test_determinant");
  ddb.createResource(res_col);
  ddb.acquireResource(res_col);
  
  ddb.recompute(ddbt.get_mw_res(ddb), elec);
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

  checkMatrix(ddb.get_det_engine().get_psiMinv(), b, std::string("ddb.psiMinv in"  __FILE__), __LINE__);


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

  checkMatrix(ddb.get_det_engine().get_psiMinv(), b, std::string("bad ddb.psiMinv after accept move in"  __FILE__), __LINE__);
}

//#define DUMP_INFO

TEST_CASE("DiracDeterminantBatched_second", "[wavefunction][fermion]")
{
  auto spo_init = std::make_unique<FakeSPO>();
  spo_init->setOrbitalSetSize(4);
  DetType ddb(std::move(spo_init));
  auto spo = dynamic_cast<FakeSPO*>(ddb.getPhi());

  int norb = 4;
  ddb.set(0, norb);

  // occurs in call to registerData
  ddb.dpsiV.resize(norb);
  ddb.d2psiV.resize(norb);


  ParticleSet elec;

  elec.create(4);

  ResourceCollection res_col("test_determinant");
  ddb.createResource(res_col);
  ddb.acquireResource(res_col);

  testing::DiracDeterminantBatchedTest ddbt;
  ddb.recompute(ddbt.get_mw_res(ddb), elec);

  using ValueMatrix = DetType::DetEngine_t::OffloadPinnedValueMatrix_t;
  ValueMatrix orig_a;
  orig_a.resize(4, 4);
  orig_a = spo->a2;

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < norb; j++)
    {
      orig_a(j, i) = spo->v2(i, j);
    }
  }

  //checkMatrix(ddb.psiMinv, b);
  DetType::DetEngine_t::DiracMatrixCompute& dm = ddb.get_det_engine().get_det_inverter();

  ValueMatrix a_update1, scratchT;
  a_update1.resize(4, 4);
  scratchT.resize(4, 4);
  a_update1 = spo->a2;
  for (int j = 0; j < norb; j++)
  {
    a_update1(j, 0) = spo->v2(0, j);
  }

  ValueMatrix a_update2;
  a_update2.resize(4, 4);
  a_update2 = spo->a2;
  for (int j = 0; j < norb; j++)
  {
    a_update2(j, 0) = spo->v2(0, j);
    a_update2(j, 1) = spo->v2(1, j);
  }

  ValueMatrix a_update3;
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
  // waiting for C++17
  // if constexpr (std::is_same<decltype(dm), DiracMatrixComputeCUDA<QMCTraits::QTFull::ValueType>>)
  //   dm.invert_transpose(ddbt.get_handles(ddb), scratchT, a_update1, det_update1);
  // else
  //   dm.invert_transpose(scratchT, a_update1, det_update1);

  auto& mw_res = ddbt.get_mw_res(ddb);

#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
  dm.invert_transpose(ddbt.getHandles(ddb), scratchT, a_update1, mw_res.log_values);
  det_update1 = mw_res.log_values[0];
#else
  dm.invert_transpose(scratchT, a_update1, det_update1);
#endif


  PsiValueType det_ratio1 = LogToValue<ValueType>::convert(det_update1 - ddb.get_log_value());
#ifdef DUMP_INFO
  std::cout << "det 0 = " << std::exp(ddb.LogValue) << std::endl;
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

#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
  dm.invert_transpose(ddbt.getHandles(ddb), scratchT, a_update2, mw_res.log_values);
  det_update2 = mw_res.log_values[0];

#else
  dm.invert_transpose(scratchT, a_update2, det_update2);
#endif
  PsiValueType det_ratio2_val = LogToValue<ValueType>::convert(det_update2 - det_update1);
#ifdef DUMP_INFO
  std::cout << "det 1 = " << std::exp(ddb.LogValue) << std::endl;
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

  #if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
  dm.invert_transpose(ddbt.getHandles(ddb), scratchT, a_update3, mw_res.log_values);
  det_update3 = mw_res.log_values[0];
#else
  dm.invert_transpose(scratchT, a_update3, det_update3);
#endif

  PsiValueType det_ratio3_val = LogToValue<ValueType>::convert(det_update3 - det_update2);
#ifdef DUMP_INFO
  std::cout << "det 2 = " << std::exp(ddb.LogValue) << std::endl;
  std::cout << "det 3 = " << std::exp(det_update3) << std::endl;
  std::cout << "det ratio 3 = " << det_ratio3 << std::endl;
#endif
  REQUIRE(det_ratio3 == ValueApprox(det_ratio3_val));
  //check_value(det_ratio3, det_ratio3_val);

  ddb.acceptMove(elec, 2);

  simd::transpose(orig_a.data(), orig_a.rows(), orig_a.cols(), scratchT.data(), scratchT.rows(), scratchT.cols());
#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
  dm.invert_transpose(ddbt.getHandles(ddb), scratchT, orig_a, mw_res.log_values);
  det_update3 = mw_res.log_values[0];
#else
  dm.invert_transpose(scratchT, orig_a, det_update3);
#endif

#ifdef DUMP_INFO
  std::cout << "original " << std::endl;
  std::cout << orig_a << std::endl;
  std::cout << "block update " << std::endl;
  std::cout << ddb.psiMinv << std::endl;
#endif

  checkMatrix(ddb.get_det_engine().get_psiMinv(), orig_a);
}

TEST_CASE("DiracDeterminantBatched_delayed_update", "[wavefunction][fermion]")
{
  auto spo_init = std::make_unique<FakeSPO>();
  spo_init->setOrbitalSetSize(4);
  DetType ddc(std::move(spo_init));
  auto spo = dynamic_cast<FakeSPO*>(ddc.getPhi());

  int norb = 4;
  // maximum delay 2
  ddc.set(0, norb, 2);

  // occurs in call to registerData
  ddc.dpsiV.resize(norb);
  ddc.d2psiV.resize(norb);


  ParticleSet elec;

  elec.create(4);

  ResourceCollection res_col("test_determinant");
  ddc.createResource(res_col);
  ddc.acquireResource(res_col);

  testing::DiracDeterminantBatchedTest ddbt;
  ddc.recompute(ddbt.get_mw_res(ddc), elec);

  using ValueMatrix = DetType::DetEngine_t::OffloadPinnedValueMatrix_t;
  ValueMatrix orig_a;
  orig_a.resize(4, 4);
  orig_a = spo->a2;

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < norb; j++)
    {
      orig_a(j, i) = spo->v2(i, j);
    }
  }

  DetType::DetEngine_t::DiracMatrixCompute& dm = ddc.get_det_engine().get_det_inverter();

  ValueMatrix a_update1, scratchT;
  scratchT.resize(4, 4);
  a_update1.resize(4, 4);
  a_update1 = spo->a2;
  for (int j = 0; j < norb; j++)
  {
    a_update1(j, 0) = spo->v2(0, j);
  }

  ValueMatrix a_update2;
  a_update2.resize(4, 4);
  a_update2 = spo->a2;
  for (int j = 0; j < norb; j++)
  {
    a_update2(j, 0) = spo->v2(0, j);
    a_update2(j, 1) = spo->v2(1, j);
  }

  ValueMatrix a_update3;
  a_update3.resize(4, 4);
  a_update3 = spo->a2;
  for (int j = 0; j < norb; j++)
  {
    a_update3(j, 0) = spo->v2(0, j);
    a_update3(j, 1) = spo->v2(1, j);
    a_update3(j, 2) = spo->v2(2, j);
  }


  auto& mw_res = ddbt.get_mw_res(ddc);
  
  ParticleSet::GradType grad;
  PsiValueType det_ratio = ddc.ratioGrad(elec, 0, grad);

  simd::transpose(a_update1.data(), a_update1.rows(), a_update1.cols(), scratchT.data(), scratchT.rows(),
                  scratchT.cols());
  LogValueType det_update1;
#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
  dm.invert_transpose(ddbt.getHandles(ddc), scratchT, a_update1, mw_res.log_values);
  det_update1 = mw_res.log_values[0];
#else
  dm.invert_transpose(scratchT, a_update1, det_update1);
#endif

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

  checkMatrix(ddc.get_det_engine().get_psiMinv(), a_update1);

  grad = ddc.evalGrad(elec, 1);

  PsiValueType det_ratio2 = ddc.ratioGrad(elec, 1, grad);
  simd::transpose(a_update2.data(), a_update2.rows(), a_update2.cols(), scratchT.data(), scratchT.rows(),
                  scratchT.cols());
  LogValueType det_update2;

  #if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
  dm.invert_transpose(ddbt.getHandles(ddc), scratchT, a_update2, mw_res.log_values);
  det_update2 = mw_res.log_values[0];
#else
  dm.invert_transpose(scratchT, a_update2, det_update2);
#endif
  PsiValueType det_ratio2_val = LogToValue<ValueType>::convert(det_update2 - det_update1);
#ifdef DUMP_INFO
  std::cout << "det 1 = " << std::exp(ddc.LogValue) << std::endl;
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
  #if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
  dm.invert_transpose(ddbt.getHandles(ddc), scratchT, a_update3, mw_res.log_values);
  det_update3 = mw_res.log_values[0];
#else
  dm.invert_transpose(scratchT, a_update3, det_update3);
#endif
  PsiValueType det_ratio3_val = LogToValue<ValueType>::convert(det_update3 - det_update2);
#ifdef DUMP_INFO
  std::cout << "det 2 = " << std::exp(ddc.LogValue) << std::endl;
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
#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
  dm.invert_transpose(ddbt.getHandles(ddc), scratchT, orig_a, mw_res.log_values);
  det_update3 = mw_res.log_values[0];
#else
  dm.invert_transpose(scratchT, orig_a, det_update3);
#endif

#ifdef DUMP_INFO
  std::cout << "original " << std::endl;
  std::cout << orig_a << std::endl;
  std::cout << "delayed update " << std::endl;
  std::cout << ddc.psiMinv << std::endl;
#endif

  // compare all the elements of psiMinv in ddc and orig_a
  checkMatrix(ddc.get_det_engine().get_psiMinv(), orig_a);
}

TEST_CASE("DiracDeterminantBatched_mw_delayed_update", "[wavefunction][fermion]")
{
  UPtrVector<DetType> dets;
  UPtrVector<ParticleSet> psets;
  int nw = 4;
  int norb = 4;

  for(int i = 0; i < nw; ++i){
    auto spo = std::make_unique<FakeSPO>();
    spo->setOrbitalSetSize(4);
    dets.emplace_back(std::make_unique<DetType>(std::move(spo)));
    // maximum delay 2
    dets.back()->set(0, norb, 2);
    dets.back()->dpsiV.resize(norb);
    dets.back()->d2psiV.resize(norb);
    psets.emplace_back(std::make_unique<ParticleSet>());
    psets.back()->create(4);
  }

  testing::DiracDeterminantBatchedTest ddbt;
  ResourceCollection collection("test resources");;
  dets[0]->createResource(collection);
  dets[0]->acquireResource(collection);
  auto& mw_res = ddbt.get_mw_res(*(dets[0]));
  mw_res.log_values.resize(4);

  RefVectorWithLeader<WaveFunctionComponent> det_refs(*dets[0]);
  RefVectorWithLeader<ParticleSet> pset_refs(*psets[0]);
  for(int iw = 0; iw < nw; ++iw){
    det_refs.push_back(*dets[iw]);
    pset_refs.push_back(*psets[iw]);
  }
  auto& det_leader = det_refs.getCastedLeader<DetType>();
  std::vector<bool> recompute_mask(det_refs.size(), true);
  det_leader.mw_recompute(det_refs, pset_refs, recompute_mask);

  std::vector<PsiValueType> ratios(nw, 0.0);
  std::vector<ParticleSet::GradType> grads_new(nw, 0.0);
  det_leader.mw_ratioGrad(det_refs, pset_refs, 0, ratios, grads_new);
  
  
//   ParticleSet elec;

//   elec.create(4);

//   testing::DiracDeterminantBatchedTest ddbt;
//   ddbt.guardMultiWalkerRes(ddc);
//   ddc.recompute(ddbt.get_mw_res(ddc), elec);

//   using ValueMatrix = DetType::DetEngine_t::OffloadPinnedValueMatrix_t;
//   ValueMatrix orig_a;
//   orig_a.resize(4, 4);
//   orig_a = spo->a2;

//   for (int i = 0; i < 3; i++)
//   {
//     for (int j = 0; j < norb; j++)
//     {
//       orig_a(j, i) = spo->v2(i, j);
//     }
//   }

//   DetType::DetEngine_t::DiracMatrixCompute& dm = ddc.get_det_engine().get_det_inverter();

//   ValueMatrix a_update1, scratchT;
//   scratchT.resize(4, 4);
//   a_update1.resize(4, 4);
//   a_update1 = spo->a2;
//   for (int j = 0; j < norb; j++)
//   {
//     a_update1(j, 0) = spo->v2(0, j);
//   }

//   ValueMatrix a_update2;
//   a_update2.resize(4, 4);
//   a_update2 = spo->a2;
//   for (int j = 0; j < norb; j++)
//   {
//     a_update2(j, 0) = spo->v2(0, j);
//     a_update2(j, 1) = spo->v2(1, j);
//   }

//   ValueMatrix a_update3;
//   a_update3.resize(4, 4);
//   a_update3 = spo->a2;
//   for (int j = 0; j < norb; j++)
//   {
//     a_update3(j, 0) = spo->v2(0, j);
//     a_update3(j, 1) = spo->v2(1, j);
//     a_update3(j, 2) = spo->v2(2, j);
//   }


//   auto& mw_res = ddbt.get_mw_res(ddc);
  
//   ParticleSet::GradType grad;
//   PsiValueType det_ratio = ddc.ratioGrad(elec, 0, grad);

//   simd::transpose(a_update1.data(), a_update1.rows(), a_update1.cols(), scratchT.data(), scratchT.rows(),
//                   scratchT.cols());
//   LogValueType det_update1;
// #if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
//   dm.invert_transpose(ddbt.getHandles(ddc), scratchT, a_update1, mw_res.log_values);
//   det_update1 = mw_res.log_values[0];
// #else
//   dm.invert_transpose(scratchT, a_update1, det_update1);
// #endif

//   PsiValueType det_ratio1 = LogToValue<ValueType>::convert(det_update1 - ddc.get_log_value());
// #ifdef DUMP_INFO
//   std::cout << "det 0 = " << std::exp(ddc.get_log_value()) << std::endl;
//   std::cout << "det 1 = " << std::exp(det_update1) << std::endl;
//   std::cout << "det ratio 1 = " << det_ratio1 << std::endl;
// #endif
//   //double det_ratio1 = 0.178276269185;

//   REQUIRE(det_ratio1 == ValueApprox(det_ratio));

//   // update of Ainv in ddc is delayed
//   ddc.acceptMove(elec, 0, true);
//   // force update Ainv in ddc using SM-1 code path
//   ddc.completeUpdates();

//   checkMatrix(ddc.get_det_engine().get_psiMinv(), a_update1);

//   grad = ddc.evalGrad(elec, 1);

//   PsiValueType det_ratio2 = ddc.ratioGrad(elec, 1, grad);
//   simd::transpose(a_update2.data(), a_update2.rows(), a_update2.cols(), scratchT.data(), scratchT.rows(),
//                   scratchT.cols());
//   LogValueType det_update2;

//   #if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
//   dm.invert_transpose(ddbt.getHandles(ddc), scratchT, a_update2, mw_res.log_values);
//   det_update2 = mw_res.log_values[0];
// #else
//   dm.invert_transpose(scratchT, a_update2, det_update2);
// #endif
//   PsiValueType det_ratio2_val = LogToValue<ValueType>::convert(det_update2 - det_update1);
// #ifdef DUMP_INFO
//   std::cout << "det 1 = " << std::exp(ddc.LogValue) << std::endl;
//   std::cout << "det 2 = " << std::exp(det_update2) << std::endl;
//   std::cout << "det ratio 2 = " << det_ratio2 << std::endl;
// #endif
//   // check ratio computed directly and the one computed by ddc with no delay
//   //double det_ratio2_val = 0.178276269185;
//   REQUIRE(det_ratio2 == ValueApprox(det_ratio2_val));

//   // update of Ainv in ddc is delayed
//   ddc.acceptMove(elec, 1, true);

//   grad = ddc.evalGrad(elec, 2);

//   PsiValueType det_ratio3 = ddc.ratioGrad(elec, 2, grad);
//   simd::transpose(a_update3.data(), a_update3.rows(), a_update3.cols(), scratchT.data(), scratchT.rows(),
//                   scratchT.cols());
//   LogValueType det_update3;
//   #if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
//   dm.invert_transpose(ddbt.getHandles(ddc), scratchT, a_update3, mw_res.log_values);
//   det_update3 = mw_res.log_values[0];
// #else
//   dm.invert_transpose(scratchT, a_update3, det_update3);
// #endif
//   PsiValueType det_ratio3_val = LogToValue<ValueType>::convert(det_update3 - det_update2);
// #ifdef DUMP_INFO
//   std::cout << "det 2 = " << std::exp(ddc.LogValue) << std::endl;
//   std::cout << "det 3 = " << std::exp(det_update3) << std::endl;
//   std::cout << "det ratio 3 = " << det_ratio3 << std::endl;
// #endif
//   // check ratio computed directly and the one computed by ddc with 1 delay
//   REQUIRE(det_ratio3 == ValueApprox(det_ratio3_val));
//   //check_value(det_ratio3, det_ratio3_val);

//   // maximal delay reached and Ainv is updated fully
//   ddc.acceptMove(elec, 2, true);
//   ddc.completeUpdates();

//   // fresh invert orig_a
//   simd::transpose(orig_a.data(), orig_a.rows(), orig_a.cols(), scratchT.data(), scratchT.rows(), scratchT.cols());
// #if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
//   dm.invert_transpose(ddbt.getHandles(ddc), scratchT, orig_a, mw_res.log_values);
//   det_update3 = mw_res.log_values[0];
// #else
//   dm.invert_transpose(scratchT, orig_a, det_update3);
// #endif

// #ifdef DUMP_INFO
//   std::cout << "original " << std::endl;
//   std::cout << orig_a << std::endl;
//   std::cout << "delayed update " << std::endl;
//   std::cout << ddc.psiMinv << std::endl;
// #endif

//   // compare all the elements of psiMinv in ddc and orig_a
//   checkMatrix(ddc.get_det_engine().get_psiMinv(), orig_a);
}

} // namespace qmcplusplus
