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


#ifndef QMCPLUSPLUS_OMPTARGET_ALIGNED_ALLOCATOR_H
#define QMCPLUSPLUS_OMPTARGET_ALIGNED_ALLOCATOR_H

#include <CPU/SIMD/aligned_allocator.hpp>
#include "OMPallocator.hpp"
#include "PinnedAllocator.h"

namespace qmcplusplus
{
template<typename T>
using OffloadAllocator = OMPallocator<T, aligned_allocator<T>>;
template<typename T>
using OffloadPinnedAllocator = OMPallocator<T, PinnedAlignedAllocator<T>>;
} // namespace qmcplusplus

#endif
