//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file ccl_wrapper.hpp
 * @brief Vendor-neutral C++ type and function aliases for NCCL, RCCL, and oneCCL C-APIs.
 */
#ifndef QMCPLUSPLUS_CCL_WRAPPER_HPP
#define QMCPLUSPLUS_CCL_WRAPPER_HPP

#include "config.h"
#if defined(ENABLE_CUDA)
#include <nccl.h>
#if NCCL_VERSION_CODE >= NCCL_VERSION(2, 31, 0)
#include <nccl_device/host.h>
#else
#include <nccl_device/core.h>
#endif
#elif defined(ENABLE_SYCL)
#include <oneapi/ccl.h>
#endif

namespace qmcplusplus
{
#if defined(ENABLE_CUDA)
using cclComm_t                        = ncclComm_t;
using cclWindow_t                      = ncclWindow_t;
using cclUniqueId                      = ncclUniqueId;
using cclResult_t                      = ncclResult_t;
constexpr auto cclSuccess              = ncclSuccess;
constexpr auto CCL_WIN_DEFAULT         = NCCL_WIN_DEFAULT;
constexpr auto cclGetErrorString       = ncclGetErrorString;
inline cclResult_t cclSetDevice(int) { return ncclSuccess; }
constexpr auto cclGetUniqueId          = ncclGetUniqueId;
constexpr auto cclCommInitRank         = ncclCommInitRank;
constexpr auto cclCommDestroy          = ncclCommDestroy;
constexpr auto cclMemAlloc             = ncclMemAlloc;
constexpr auto cclMemFree              = ncclMemFree;
constexpr auto cclCommWindowRegister   = ncclCommWindowRegister;
constexpr auto cclCommWindowDeregister = ncclCommWindowDeregister;
constexpr auto cclGetPeerDevicePointer = ncclGetPeerDevicePointer;
#elif defined(ENABLE_SYCL)
using cclComm_t                        = onecclComm_t;
using cclWindow_t                      = onecclWindow_t;
using cclUniqueId                      = onecclUniqueId;
using cclResult_t                      = onecclResult_t;
constexpr auto cclSuccess              = onecclSuccess;
constexpr auto CCL_WIN_DEFAULT         = ONECCL_WINDOW_DEFAULT;
constexpr auto cclGetErrorString       = onecclGetErrorString;
constexpr auto cclSetDevice            = onecclSetDevice;
constexpr auto cclGetUniqueId          = onecclGetUniqueId;
constexpr auto cclCommInitRank         = onecclCommInitRank;
constexpr auto cclCommDestroy          = onecclCommDestroy;
constexpr auto cclMemAlloc             = onecclMemAlloc;
constexpr auto cclMemFree              = onecclMemFree;
constexpr auto cclCommWindowRegister   = onecclCommWindowRegister;
constexpr auto cclCommWindowDeregister = onecclCommWindowDeregister;
constexpr auto cclGetPeerDevicePointer = onecclGetPeerDevicePointer;
#endif
} // namespace qmcplusplus

#endif
