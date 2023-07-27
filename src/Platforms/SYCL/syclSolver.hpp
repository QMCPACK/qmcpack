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

#ifndef QMCPLUSPLUS_SYCL_MKL_SOLVER_H
#define QMCPLUSPLUS_SYCL_MKL_SOLVER_H

#include <oneapi/mkl/lapack.hpp>
#include <mkl_service.h>

namespace qmcplusplus
{
namespace syclSolver
{
using namespace oneapi::mkl::lapack;

inline void freeBuffer() { mkl_free_buffers(); }
} // namespace syclSolver
} // namespace qmcplusplus

#endif
