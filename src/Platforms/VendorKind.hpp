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


#ifndef QMCPLUSPLUS_VENDOR_KIND_H
#define QMCPLUSPLUS_VENDOR_KIND_H

#include "config.h"
#include "Common/PlatformKinds.hpp"

namespace qmcplusplus
{
#if defined(ENABLE_SYCL)
  constexpr auto VendorKind = PlatformKind::SYCL;
#elif defined(ENABLE_CUDA)
  constexpr auto VendorKind = PlatformKind::CUDA;
#endif
} // namespace qmcplusplus
#endif
