//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// Refactored from: OMPallocator.cpp
//////////////////////////////////////////////////////////////////////////////////////


#include <cstddef>
#include <atomic>

namespace qmcplusplus
{
  std::atomic<size_t> dual_device_mem_allocated(0);
}
