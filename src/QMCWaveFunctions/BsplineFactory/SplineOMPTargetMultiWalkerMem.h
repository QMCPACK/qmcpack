//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_OFFLOADSHAREDMEM_H
#define QMCPLUSPLUS_OFFLOADSHAREDMEM_H

#include "OMPTarget/OffloadAlignedAllocators.hpp"
#include "ResourceCollection.h"

namespace qmcplusplus
{

template<typename ST, typename TT>
struct SplineOMPTargetMultiWalkerMem : public Resource
{
  ///team private ratios for reduction, numVP x numTeams
  Matrix<TT, OffloadPinnedAllocator<TT>> mw_ratios_private;
  ///team private ratios and grads for reduction, numVP x numTeams
  Matrix<TT, OffloadPinnedAllocator<TT>> rg_private;
  ///offload scratch space, dynamically resized to the maximal need
  Vector<ST, OffloadPinnedAllocator<ST>> mw_offload_scratch;
  ///result scratch space, dynamically resized to the maximal need
  Vector<TT, OffloadPinnedAllocator<TT>> mw_results_scratch;
  ///position scratch space, used to avoid allocation on the fly and faster transfer
  Vector<ST, OffloadPinnedAllocator<ST>> mw_pos_copy;
  ///multi purpose H2D buffer for mw_evaluateVGLandDetRatioGrads
  Matrix<char, OffloadPinnedAllocator<char>> buffer_H2D;
  ///multi purpose H2D buffer for mw_evaluateDetRatios
  Vector<char, OffloadPinnedAllocator<char>> det_ratios_buffer_H2D;

  SplineOMPTargetMultiWalkerMem() : Resource("SplineOMPTargetMultiWalkerMem") {}

  SplineOMPTargetMultiWalkerMem(const SplineOMPTargetMultiWalkerMem&) : SplineOMPTargetMultiWalkerMem() {}

  std::unique_ptr<Resource> makeClone() const override
  {
    return std::make_unique<SplineOMPTargetMultiWalkerMem>(*this);
  }
};
} // namespace qmcplusplus
#endif
