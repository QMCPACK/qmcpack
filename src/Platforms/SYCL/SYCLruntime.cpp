//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include <CL/sycl.hpp>
#include "SYCLDeviceManager.h"
#include "SYCLruntime.hpp"

namespace qmcplusplus
{
sycl::queue getSYCLDefaultDeviceDefaultQueue() { return SYCLDeviceManager::getDefaultDeviceQueue(); }
} // namespace qmcplusplus
