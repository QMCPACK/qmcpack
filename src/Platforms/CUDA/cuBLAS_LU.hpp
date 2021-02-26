//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_CUBLAS_LU_HPP
#define QMCPLUSPLUS_CUBLAS_LU_HPP

#include <complex>
#include <cuda_runtime_api.h>

namespace qmcplusplus
{
namespace QMC_CUDA
{
void computeInverseAndDetLog_batched(cublasHandle_t& h_cublas,
                                     const int n,
                                     const int lda,
                                     double** Ms,
                                     double* LU_diags,
                                     int* pivots,
                                     int* infos,
                                     double* log_dets,
                                     const int batch_size);
} // namespace QMC_CUDA
} // namespace qmcplusplus
#endif
