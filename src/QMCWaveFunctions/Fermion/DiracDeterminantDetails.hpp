//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DIRAC_DETERMINANT_DETAIL_HPP
#define QMCPLUSPLUS_DIRAC_DETERMINANT_DETAIL_HPP

#include "Configuration.h"
#include "Platforms/PinnedAllocator.h"
#include "OMPTarget/OMPallocator.hpp"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "WaveFunctionComponent.h"
#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
#include "MatrixDelayedUpdateCUDA.h"
#endif

namespace qmcplusplus
{

  /** cutout class for explicit specialization of virtual methods in DiracDeterminantBatched
   *
   *  perhaps something better can be done in C++17 or C++20
   */
template<class DET, class DET_ENGINE>
class DiracDeterminantDetails
{
public:
  using ValueType = QMCTraits::QTFull::ValueType;
  template<typename DT>
  using OffloadPinnedAllocator     = OMPallocator<DT, PinnedAlignedAllocator<DT>>;
  using OffloadPinnedValueVector_t = Vector<ValueType, OffloadPinnedAllocator<ValueType>>;
  using OffloadPinnedValueMatrix_t = Matrix<ValueType, OffloadPinnedAllocator<ValueType>>;

  static void mw_recomputeDispatch(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                   const RefVectorWithLeader<ParticleSet>& p_list,
                                   const std::vector<bool>& recompute);
  // template<class DET_ENGINE>
  // static void mw_invertPsiM(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
  //                           const RefVector<const OffloadPinnedValueMatrix_t>& logdetT_list);
};
  
} // namespace qmcplusplus

#endif
