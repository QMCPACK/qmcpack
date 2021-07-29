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
// -*- C++ -*-
/** @file
 */
#ifndef QMCPLUSPLUS_OMPTARGET_ALIGNED_ALLOCATOR_H
#define QMCPLUSPLUS_OMPTARGET_ALIGNED_ALLOCATOR_H


#include "OMPallocator.hpp"
#include "PinnedAllocator.h"

namespace qmcplusplus
{
  template<typename DT>
  using OffloadAllocator = OMPallocator<DT, aligned_allocator<DT>>;
  template<typename DT>
  using OffloadPinnedAllocator = OMPallocator<DT, PinnedAlignedAllocator<DT>>;
} // namespace qmcplusplus
#endif
