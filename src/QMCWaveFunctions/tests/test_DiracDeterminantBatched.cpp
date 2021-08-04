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
#include "Utilities/for_testing/checkMatrix.hpp"

#ifdef QMC_COMPLEX //This is for the spinor test.
#include "QMCWaveFunctions/ElectronGas/ElectronGasComplexOrbitalBuilder.h"
#endif

#include "QMCWaveFunctions/SpinorSet.h"

#include <type_traits>
#include <stdio.h>
#include <string>

using std::string;

#define DUMP_INFO

namespace qmcplusplus
{
using RealType     = QMCTraits::RealType;
using ValueType    = QMCTraits::ValueType;
using LogValueType = std::complex<QMCTraits::QTFull::RealType>;
using PsiValueType = QMCTraits::QTFull::ValueType;

#if defined(ENABLE_CUDA)
using DetType = DiracDeterminantBatched<MatrixDelayedUpdateCUDA<ValueType, QMCTraits::QTFull::ValueType>>;
#elif defined(ENABLE_OFFLOAD)
using DetType = DiracDeterminantBatched<MatrixUpdateOMPTarget<ValueType, QMCTraits::QTFull::ValueType>>;
#else
using DetType = DiracDeterminantBatched<MatrixUpdateOMPTarget<ValueType, QMCTraits::QTFull::ValueType>>;
#endif

namespace testing
{
template<class T, typename = void>
struct has_handles : std::false_type
{};

// For det engines using different handle member handle accessors for Ts and
// checking for that instead of cuda_handles_ would work.
//struct impl_has_handles<T, void_t<decltype(T::cuda_handles_)>> : public std::true_type
template<class T>
struct has_handles<T, decltype(std::declval<T>().cuda_handles_, void())> : std::true_type
{};

template<bool B>
struct HandleTraits;

class DiracDeterminantBatchedTest
{
public:
  template<typename DT>
  using OffloadPinnedAllocator = OMPallocator<DT, PinnedAlignedAllocator<DT>>;
  template<typename T>
  using OffloadPinnedVector = Vector<T, OffloadPinnedAllocator<T>>;
  template<typename T>
  using OffloadPinnedMatrix = Matrix<T, OffloadPinnedAllocator<T>>;

  template<class DET_ENGINE>
  DiracDeterminantBatchedMultiWalkerResource& get_mw_res(DiracDeterminantBatched<DET_ENGINE>& ddb)
  {
    return *(ddb.mw_res_);
  }
  template<class DET_ENGINE>
  void invertPsiM(DiracDeterminantBatched<DET_ENGINE>& ddb,
                  OffloadPinnedMatrix<typename DiracDeterminantBatched<DET_ENGINE>::ValueType>& logdetT,
                  OffloadPinnedMatrix<typename DiracDeterminantBatched<DET_ENGINE>::ValueType>& a_inv)
  {
    ddb.invertPsiM(*ddb.mw_res_, logdetT, a_inv);
  }
  template<class DET_ENGINE>
  Resource& getHandles(DET_ENGINE& det_engine)
  {
    return HandleTraits<has_handles<DET_ENGINE>::value>::getHandle(*this, det_engine);
  }

  DummyResource dummy_;
};

template<>
struct HandleTraits<true>
{
  template<class DET_ENGINE>
  static Resource& getHandle(DiracDeterminantBatchedTest& ddbt, DET_ENGINE& det_engine)
  {
    return det_engine.getHandles();
  }
};

template<>
struct HandleTraits<false>
{
  template<class DET_ENGINE>
  static Resource& getHandle(DiracDeterminantBatchedTest& ddbt, DET_ENGINE& ddb)
  {
    return ddbt.dummy_;
  }
};

} // namespace testing

// I've observed downstream unit tests fail due to small issues in resource create/acquire/release.
// unit test these.
TEST_CASE("DiracDeterminantBatched_resources", "[wavefunction][fermion]")
{
  int count = 0;

  auto spo_init = std::make_unique<FakeSPO>();
  spo_init->setOrbitalSetSize(3);
  DetType ddb(std::move(spo_init));
  auto spo = dynamic_cast<FakeSPO*>(ddb.getPhi());

  ResourceCollection det_res("test resources");
  ddb.createResource(det_res);
  RefVectorWithLeader<WaveFunctionComponent> det_ref_list{ddb, {ddb}};
  ResourceCollectionTeamLock<WaveFunctionComponent> lock{det_res, det_ref_list};
}

TEST_CASE("DiracDeterminantBatched_mw_delayed_update", "[wavefunction][fermion]")
{
  UPtrVector<DetType> dets;
  UPtrVector<ParticleSet> psets;

  int nw   = 4;
  int norb = 4;

  for (int i = 0; i < nw; ++i)
  {
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
  ResourceCollection det_res("test resources");
  dets[0]->createResource(det_res);
  RefVectorWithLeader<WaveFunctionComponent> det_ref_list{*dets[0], convertUPtrToRefVector<WaveFunctionComponent>(dets)};
  ResourceCollectionTeamLock<WaveFunctionComponent> det_lock{det_res, det_ref_list};
  auto& mw_res = ddbt.get_mw_res(*(dets[0]));
  mw_res.log_values.resize(4);

  RefVectorWithLeader<WaveFunctionComponent> det_refs(*dets[0]);
  RefVectorWithLeader<ParticleSet> pset_refs(*psets[0]);
  for (int iw = 0; iw < nw; ++iw)
  {
    det_refs.push_back(*dets[iw]);
    pset_refs.push_back(*psets[iw]);
  }
  auto& det_leader = det_refs.getCastedLeader<DetType>();
  std::vector<bool> recompute_mask(det_refs.size(), true);
  det_leader.mw_recompute(det_refs, pset_refs, recompute_mask);

  std::vector<PsiValueType> ratios(nw, 0.0);
  std::vector<ParticleSet::GradType> grads_new(nw);
  det_leader.mw_ratioGrad(det_refs, pset_refs, 0, ratios, grads_new);

  // using ValueMatrix = DetType::DetEngine_t::OffloadPinnedValueMatrix_t;
  // ValueMatrix orig_a;
  // orig_a.resize(4, 4);
  // orig_a = spo->a2;

  // for (int i = 0; i < 3; i++)
  // {
  //   for (int j = 0; j < norb; j++)
  //   {
  //     orig_a(j, i) = spo->v2(i, j);
  //   }
  // }

  // DetType::DetEngine_t::DiracMatrixCompute& dm = ddc.get_det_engine().get_det_inverter();

  // ValueMatrix a_update1, scratchT;
  // scratchT.resize(4, 4);
  // a_update1.resize(4, 4);
  // a_update1 = spo->a2;
  // for (int j = 0; j < norb; j++)
  // {
  //   a_update1(j, 0) = spo->v2(0, j);
  // }

  // ValueMatrix a_update2;
  // a_update2.resize(4, 4);
  // a_update2 = spo->a2;
  // for (int j = 0; j < norb; j++)
  // {
  //   a_update2(j, 0) = spo->v2(0, j);
  //   a_update2(j, 1) = spo->v2(1, j);
  // }

  // ValueMatrix a_update3;
  // a_update3.resize(4, 4);
  // a_update3 = spo->a2;
  // for (int j = 0; j < norb; j++)
  // {
  //   a_update3(j, 0) = spo->v2(0, j);
  //   a_update3(j, 1) = spo->v2(1, j);
  //   a_update3(j, 2) = spo->v2(2, j);
  // }


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
