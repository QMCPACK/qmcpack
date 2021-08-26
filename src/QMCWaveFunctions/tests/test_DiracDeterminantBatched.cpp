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
  typename DiracDeterminantBatched<DET_ENGINE>::DiracDeterminantBatchedMultiWalkerResource& get_mw_res(DiracDeterminantBatched<DET_ENGINE>& ddb)
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
  auto spo_init = std::make_unique<FakeSPO>();
  spo_init->setOrbitalSetSize(3);
  DetType ddb(std::move(spo_init));

  ResourceCollection det_res("test resources");
  ddb.createResource(det_res);
  RefVectorWithLeader<WaveFunctionComponent> det_ref_list{ddb, {ddb}};
  ResourceCollectionTeamLock<WaveFunctionComponent> lock{det_res, det_ref_list};
}

} // namespace qmcplusplus
