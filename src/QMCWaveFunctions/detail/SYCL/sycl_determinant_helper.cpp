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


#include "sycl_determinant_helper.hpp"

namespace qmcplusplus
{
template<typename T>
sycl::event applyW_stageV_sycl(sycl::queue& aq,
                               const std::vector<sycl::event>& dependencies,
                               const int* restrict delay_list_gpu,
                               const int delay_count,
                               T* restrict temp_gpu,
                               const int numorbs,
                               const int ndelay,
                               T* restrict V_gpu,
                               const T* restrict Ainv)
{
  const size_t BS = 128;
  const size_t NB = (numorbs + BS - 1) / BS;

  return aq.parallel_for(sycl::nd_range<1>{{BS * NB}, {BS}}, dependencies, [=](sycl::nd_item<1> item) {
    int col = item.get_global_id(0);

    // move rows of Ainv to V
    for (int row = 0; row < delay_count; row++)
    {
      const T* Ainv_row = Ainv + numorbs * delay_list_gpu[row];
      T* V_row          = V_gpu + numorbs * row;
      if (col < numorbs)
        V_row[col] = Ainv_row[col];
    }

    // apply W to temp
    if (col < delay_count)
      temp_gpu[ndelay * delay_list_gpu[col] + col] -= T(1);
  });
}

template sycl::event applyW_stageV_sycl(sycl::queue& aq,
                                        const std::vector<sycl::event>& dependencies,
                                        const int* restrict delay_list_gpu,
                                        const int delay_count,
                                        float* restrict temp_gpu,
                                        const int numorbs,
                                        const int ndelay,
                                        float* restrict V_gpu,
                                        const float* restrict Ainv);

template sycl::event applyW_stageV_sycl(sycl::queue& aq,
                                        const std::vector<sycl::event>& dependencies,
                                        const int* restrict delay_list_gpu,
                                        const int delay_count,
                                        double* restrict temp_gpu,
                                        const int numorbs,
                                        const int ndelay,
                                        double* restrict V_gpu,
                                        const double* restrict Ainv);

template sycl::event applyW_stageV_sycl(sycl::queue& aq,
                                        const std::vector<sycl::event>& dependencies,
                                        const int* restrict delay_list_gpu,
                                        const int delay_count,
                                        std::complex<float>* restrict temp_gpu,
                                        const int numorbs,
                                        const int ndelay,
                                        std::complex<float>* restrict V_gpu,
                                        const std::complex<float>* restrict Ainv);

template sycl::event applyW_stageV_sycl(sycl::queue& aq,
                                        const std::vector<sycl::event>& dependencies,
                                        const int* restrict delay_list_gpu,
                                        const int delay_count,
                                        std::complex<double>* restrict temp_gpu,
                                        const int numorbs,
                                        const int ndelay,
                                        std::complex<double>* restrict V_gpu,
                                        const std::complex<double>* restrict Ainv);


template<typename T, typename TMAT, typename INDEX>
std::complex<T> computeLogDet_sycl(sycl::queue& aq,
                                   int n,
                                   int lda,
                                   const TMAT* restrict a,
                                   const INDEX* restrict pivot,
                                   const std::vector<sycl::event>& dependencies)
{
  constexpr size_t COLBS = 128;

  std::complex<T> result{};
  {
    sycl::buffer<std::complex<T>, 1> abuff(&result, {1});
    aq.submit([&](sycl::handler& cgh) {
      cgh.depends_on(dependencies);

      size_t n_max = ((n + COLBS - 1) / COLBS) * COLBS;
      sycl::global_ptr<const TMAT> A{a};
      sycl::global_ptr<const INDEX> Pivot{pivot};
      cgh.parallel_for(sycl::range<1>{n_max}, sycl::reduction(abuff, cgh, {T{}, T{}}, std::plus<std::complex<T>>()),
                       [=](sycl::id<1> i, auto& sum) {
                         std::complex<T> val{};
                         if (i < n)
                           val = (Pivot[i] == i + 1) ? std::log(std::complex<T>(A[i * lda + i]))
                                                     : std::log(std::complex<T>(-A[i * lda + i]));
                         sum.combine(val);
                       });
    });
  } //synchronous
  return result;
}

template std::complex<double> computeLogDet_sycl(sycl::queue& aq,
                                                 int n,
                                                 int lda,
                                                 const double* restrict a,
                                                 const std::int64_t* restrict pivot,
                                                 const std::vector<sycl::event>& dependencies);

template std::complex<double> computeLogDet_sycl(sycl::queue& aq,
                                                 int n,
                                                 int lda,
                                                 const std::complex<double>* restrict a,
                                                 const std::int64_t* restrict pivot,
                                                 const std::vector<sycl::event>& dependencies);
} // namespace qmcplusplus
