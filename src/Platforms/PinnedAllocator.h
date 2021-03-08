//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_PINNED_ALLOCATOR_H
#define QMCPLUSPLUS_PINNED_ALLOCATOR_H

#include <memory>
#include "CPU/SIMD/aligned_allocator.hpp"
#ifdef ENABLE_CUDA
#include "CUDA/CUDAallocator.hpp"
#endif

namespace qmcplusplus
{

template<typename T>
#ifdef ENABLE_CUDA
using PinnedAllocator = CUDALockedPageAllocator<T>;
#else
using PinnedAllocator = std::allocator<T>;
#endif

template<typename T, size_t ALIGN = QMC_SIMD_ALIGNMENT>
#ifdef ENABLE_CUDA
using PinnedAlignedAllocator = CUDALockedPageAllocator<T, aligned_allocator<T, ALIGN>>;
#else
using PinnedAlignedAllocator = aligned_allocator<T, ALIGN>;
#endif

} // namespace qmcplusplus

#endif
