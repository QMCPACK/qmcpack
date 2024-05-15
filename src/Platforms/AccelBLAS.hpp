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

#include "PlatformKinds.hpp"
#include "Queue.hpp"

namespace qmcplusplus
{

namespace compute
{

template<PlatformKind PL>
class BLASHandle;

namespace BLAS
{
template<PlatformKind PL, typename T>
void gemm(BLASHandle<PL>& queue,
const char transa,
                           const char transb,
                           int m,
                           int n,
                           int k,
                           const T* alpha,
                           const T* A,
                           int lda,
                           const T* B,
                           int ldb,
                           const T* beta,
                           T* C,
                           int ldc);
}
}

} // namespace qmcplusplus

#endif
