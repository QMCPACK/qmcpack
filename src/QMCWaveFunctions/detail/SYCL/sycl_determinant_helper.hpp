//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SYCL_DETERMINANT_HELPER_H
#define QMCPLUSPLUS_SYCL_DETERMINANT_HELPER_H

#include <complex>
#include <vector>
#include <sycl/sycl.hpp>

namespace qmcplusplus
{
template<typename T>
sycl::event applyW_stageV_sycl(sycl::queue& aq,
                               const std::vector<sycl::event>& dependencies,
                               const int* delay_list_gpu,
                               const int delay_count,
                               T* temp_gpu,
                               const int numorbs,
                               const int ndelay,
                               T* V_gpu,
                               const T* Ainv);

template<typename T, typename TMAT, typename INDEX>
std::complex<T> computeLogDet_sycl(sycl::queue& aq,
                                   int n,
                                   int lda,
                                   const TMAT* a,
                                   const INDEX* pivot,
                                   const std::vector<sycl::event>& dependencies = {});

} // namespace qmcplusplus
#endif
