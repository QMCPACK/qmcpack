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


#include "MemManageSYCL.hpp"

namespace qmcplusplus
{
namespace compute
{
MemoryUsageAccount MemManage<PlatformKind::SYCL>::device_mem_usage_;
template class MemManage<PlatformKind::SYCL>;
} // namespace compute
} // namespace qmcplusplus
