//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ACCELBLAS_H
#define QMCPLUSPLUS_ACCELBLAS_H

#include "config.h"
#if defined(ENABLE_CUDA)
#include "CUDA/AccelBLAS_CUDA.hpp"
#endif
#if defined(ENABLE_SYCL)
#include "SYCL/AccelBLAS_SYCL.hpp"
#endif
#include "OMPTarget/AccelBLAS_OMPTarget.hpp"

#endif
