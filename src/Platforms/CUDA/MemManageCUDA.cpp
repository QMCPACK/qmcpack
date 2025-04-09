//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "MemManageCUDA.hpp"

namespace qmcplusplus
{
std::atomic<size_t> CUDAallocator_device_mem_allocated(0);
namespace compute
{
template class MemManage<PlatformKind::CUDA>;
} // namespace compute
} // namespace qmcplusplus
