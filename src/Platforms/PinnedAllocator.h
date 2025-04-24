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
#include "VendorKind.hpp"
#include "MemManageAlias.hpp"

namespace qmcplusplus
{
/** The fact that the pinned allocators are not always pinned hurts readability elsewhere. */
#if defined(ENABLE_CUDA) || defined(ENABLE_SYCL)
template<typename T>
using PinnedAllocator = compute::MemManage<VendorKind>::PageLockedAllocator<T>;
template<typename T, size_t ALIGN = QMC_SIMD_ALIGNMENT>
using PinnedAlignedAllocator = compute::MemManage<VendorKind>::PageLockedAllocator<T, aligned_allocator<T, ALIGN>>;
#else
template<typename T>
using PinnedAllocator = std::allocator<T>;
template<typename T, size_t ALIGN = QMC_SIMD_ALIGNMENT>
using PinnedAlignedAllocator = aligned_allocator<T, ALIGN>;
#endif
} // namespace qmcplusplus

#endif
