//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SYCL_DETERMINANT_HELPER_H
#define QMCPLUSPLUS_SYCL_DETERMINANT_HELPER_H

#include <complex>
#include <vector>
#include <CL/sycl.hpp>

namespace qmcplusplus
{

template<typename T>
sycl::event applyW_stageV_sycl(sycl::queue& aq,
                               const std::vector<cl::sycl::event>& dependencies,
                               const int* restrict delay_list_gpu,
                               const int delay_count,
                               T* restrict temp_gpu,
                               const int numorbs,
                               const int ndelay,
                               T* restrict V_gpu,
                               const T* restrict Ainv);
}
#endif
