//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// Filef developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file
 *  These allocators are to make code that should be 
 *  generic with the respect to accelerator code flavor actually so,
 *  but only through configuration time switches.
 *  A DualAllocator as in DualAllocator.hpp constructed of a Host and Device allocator
 *  or a OMPallocator which leverages the OMP runtime magic to map host and implicit device data
 *  if offload is enabled or is just a host allocator otherwise.
 */

#ifndef QMCPLUSPLUS_DUAL_ALLOCATOR_ALIASES_HPP
#define QMCPLUSPLUS_DUAL_ALLOCATOR_ALIASES_HPP

#include "PinnedAllocator.h"
#if (defined(ENABLE_CUDA) || defined(ENABLE_SYCL)) && !defined(ENABLE_OFFLOAD)
#include "DualAllocator.hpp"
#if defined(ENABLE_CUDA)
namespace qmcplusplus
{
  template<typename T>
  using UnpinnedDualAllocator = DualAllocator<T, CUDAAllocator<T>, aligned_allocator<T>>;
  template<typename T>
  using PinnedDualAllocator = DualAllocator<T, CUDAAllocator<T>, PinnedAlignedAllocator<T>>;
}
#elif defined(ENABLE_SYCL)
namespace qmcplusplus
{
  template<typename T>
  using UnpinnedDualAllocator = DualAllocator<T, SYCLAllocator<T>, aligned_allocator<T>>;
  template<typename T>
  using PinnedDualAllocator = DualAllocator<T, SYCLAllocator<T>, PinnedAlignedAllocator<T>>;
}
#else
#error unhandled platform
#endif

#else // ENABLE_OFFLOAD or no CUDA or SYCL
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
