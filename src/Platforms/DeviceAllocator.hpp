//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Youngjun Lee, leey@anl.gov, Argonne National Laboratory
//
// File created by: Youngjun Lee, leey@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file
 *  This file provides a generic interface for device allocators
 *  that can be used with different accelerator backends (CUDA, SYCL, OMPTarget).
 */

#ifndef QMCPLUSPLUS_DEVICEALLOCATOR_HPP
#define QMCPLUSPLUS_DEVICEALLOCATOR_HPP

#include <config.h>

#if (defined(ENABLE_CUDA) || defined(ENABLE_SYCL)) && !defined(ENABLE_OFFLOAD)

#if defined(ENABLE_CUDA)
#include "MemManageAlias.hpp"
namespace qmcplusplus
{
template<typename T>
using DeviceAllocator = CUDAAllocator<T>;
}

#elif defined(ENABLE_SYCL)
#include "SYCL/SYCLallocator.hpp"
namespace qmcplusplus
{
template<typename T>
using DeviceAllocator = SYCLAllocator<T>;
}

#else
#error unhandled platform
#endif

#else // ENABLE_OFFLOAD or no CUDA or SYCL

#include "OMPTarget/OffloadAlignedAllocators.hpp"
namespace qmcplusplus
{
template<typename T>
using DeviceAllocator = OffloadDeviceAllocator<T>;
}
#endif

#endif
