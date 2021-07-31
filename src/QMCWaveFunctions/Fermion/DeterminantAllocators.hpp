//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// Refactored from: OMPallocator.hpp
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file
 *  These allocators are problematic to maintain consistency of
 */

#ifndef QMCPLUSPLUS_DETERMINANT_ALLOCATORS_HPP
#define QMCPLUSPLUS_DETERMINANT_ALLOCATORS_HPP

#include "Platforms/PinnedAllocator.h"

#if defined(ENABLE_CUDA) && ! defined(ENABLE_OFFLOAD)
#include "DualAllocator.hpp"
namespace qmcplusplus
{
  template<typename T>
  using UnpinnedDualAllocator = DualAllocator<T, CUDAAllocator<T>, aligned_allocator<T>>;
  template<typename T>
  using PinnedDualAllocator = DualAllocator<T, CUDAAllocator<T>, PinnedAlignedAllocator<T>>;
}
#else
#include "OMPTarget/OffloadAlignedAllocators.hpp"
namespace qmcplusplus
{
  template<typename T>
  using UnpinnedDualAllocator = OffloadAllocator<T>;
  template<typename T>
  using PinnedDualAllocator = OffloadPinnedAllocator<T>;
}
#endif

#endif
