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
#include "QMCWaveFunctions/detail/SYCL/sycl_determinant_helper.hpp"

namespace qmcplusplus
{
  template<>
    sycl::event 
    applyW_stageV_sycl(sycl::queue& aq,
                       const int* restrict delay_list_gpu, const int delay_count,
                       float* restrict temp_gpu, const int numorbs, const int ndelay,
                       float* restrict V_gpu, const float* restrict Ainv);

  template<>
    sycl::event 
    applyW_stageV_sycl(sycl::queue& aq,
                       const int* restrict delay_list_gpu, const int delay_count,
                       double* restrict temp_gpu, const int numorbs, const int ndelay,
                       double* restrict V_gpu, const double* restrict Ainv);

  template<>
    sycl::event 
    applyW_stageV_sycl(sycl::queue& aq,
                       const int* restrict delay_list_gpu, const int delay_count,
                       std::complex<float>* restrict temp_gpu, const int numorbs, const int ndelay,
                       std::complex<float>* restrict V_gpu, const std::complex<float>* restrict Ainv);

  template<>
    sycl::event 
    applyW_stageV_sycl(sycl::queue& aq,
                       const int* restrict delay_list_gpu, const int delay_count,
                       std::complex<double>* restrict temp_gpu, const int numorbs, const int ndelay,
                       std::complex<double>* restrict V_gpu, const std::complex<double>* restrict Ainv);
}
