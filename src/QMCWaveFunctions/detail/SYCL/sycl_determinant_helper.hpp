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

namespace qmcplusplus
{
template<typename T, typename TMAT, typename index_t>
  inline std::complex<T> computeLogDet(int n, 
                                       int lda, 
                                       const TMAT* restrict inv_mat, 
                                       const index_t* restrict pivot)
  {
    std::complex<T> logdet{};
    for (size_t i = 0; i < n; i++)
      logdet += std::log(std::complex<T>((pivot[i] == i + 1) ? inv_mat[i*lda+i] : -inv_mat[i*lda+i]));
    return logdet;
  }

/** compute logdet of a single array using SYCL::reduction 
 */
template<typename T, typename TMAT, typename Index_t>
inline std::complex<T> computeLogDet(sycl::queue& aq, 
                                     int n, 
                                     int lda, 
                                     const TMAT* restrict a, 
                                     const Index_t* restrict pivot)
  {
    constexpr size_t COLBS=32;

    std::complex<T> result{};
    {
      sycl::buffer<std::complex<T>,1> abuff(&result,{1});
      aq.submit([&](sycl::handler& cgh) {
          size_t n_max=((n+COLBS-1)/COLBS)*COLBS;
          sycl::global_ptr<const TMAT>  A{a};
          sycl::global_ptr<const Index_t>  Pivot{pivot};
          cgh.parallel_for(sycl::range<1>{n_max},
              sycl::reduction(abuff,cgh,{T{},T{}},std::plus<std::complex<T>>()), 
              [=](sycl::id<1> i, auto& sum)
              {
                  std::complex<T> val{};
                  //if(i<n) val= std::log(std::complex<T>((Pivot[i] == i + 1) ? A[i*lda+i] : -A[i*lda+i]));
                  if(i<n) 
                  val = (Pivot[i] == i + 1) ? 
                    std::log(std::complex<T>(A[i*lda+i])) : std::log(std::complex<T>(-A[i*lda+i]));
                  sum.combine(val);
              });
          });
    } //synchronous
    return  result;
  }

template<typename T, typename TMAT, typename Index_t>
  inline std::complex<T> computeLogDetNDR(sycl::queue& aq, 
      int n, int lda, 
      const TMAT* restrict a, 
      const Index_t* restrict pivot,
      const size_t COLBS=256)
  {
    std::complex<T> result{};
    {
      sycl::buffer<std::complex<T>,1> abuff(&result,{1});
      aq.submit([&](sycl::handler& cgh) {
          size_t n_max=((n+COLBS-1)/COLBS)*COLBS;
          sycl::global_ptr<const TMAT>  A{a};
          sycl::global_ptr<const Index_t>  Pivot{pivot};
          cgh.parallel_for(sycl::nd_range<1>{{size_t(n_max)},{COLBS}}, 
              sycl::reduction(abuff,cgh,{T{},T{}},std::plus<std::complex<T>>()), 
              [=](sycl::nd_item<1> item, auto& sum)
              {
                  const unsigned i=item.get_global_id(0);
                  std::complex<T> val{};
                  //if(i<n) val= std::log(std::complex<T>((Pivot[i] == i + 1) ? A[i*lda+i] : -A[i*lda+i]));
                  if(i<n) 
                  val = (Pivot[i] == i + 1) ? 
                    std::log(std::complex<T>(A[i*lda+i])) : std::log(std::complex<T>(-A[i*lda+i]));
                  sum.combine(val);
              });
          });
    } 
    return  result;
  }

template<typename TMAT, typename T, typename Index_t>
  inline sycl::event computeLogDet_batched(sycl::queue& aq, 
                                   int n, 
                                   int lda, 
                                   const TMAT* restrict mat_lus, 
                                   const Index_t* restrict pivots,
                                   std::complex<T>* restrict logdets,
                                   const int batch_size,
                                   const size_t COLBS=128)
  {
    return aq.submit([&](sycl::handler& cgh) {

        sycl::accessor<std::complex<T>, 1, sycl::access::mode::write, sycl::access::target::local> 
        logdet_vals(sycl::range<1>{COLBS}, cgh);

        cgh.parallel_for(sycl::nd_range<1>{{batch_size*COLBS},{COLBS}}, 
            [=](sycl::nd_item<1> item) {
            const unsigned iw  = item.get_group(0);
            const unsigned tid = item.get_local_id(0);
            const Index_t* restrict pv_iw = pivots+iw*lda;
            const TMAT* restrict    lu_iw = mat_lus+iw*n*lda;

            std::complex<T> val{};
            const unsigned  num_col_blocks = (n + COLBS - 1) / COLBS;
            for (unsigned block_num = 0; block_num < num_col_blocks; block_num++)
            {
              const unsigned i = tid + block_num * COLBS;
              if(i<n)
                val += (pv_iw[i] == i + 1) ? 
                  std::log(std::complex<T>(lu_iw[i*lda+i])) : std::log(std::complex<T>(-lu_iw[i*lda+i]));
            }

            logdet_vals[tid]=val;

            for (unsigned iend = COLBS / 2; iend > 0; iend /= 2)
            {
              item.barrier(sycl::access::fence_space::local_space);
              if (tid < iend)
              {
                 logdet_vals[tid] += logdet_vals[tid + iend];
              }
            }
            if (tid == 0)
              logdets[iw] = logdet_vals[0];
            });
    });
  }

template<typename TMAT, typename T, typename Index_t>
  inline sycl::event computeLogDetGroup(sycl::queue& aq, 
                                   int n, 
                                   int lda, 
                                   const TMAT* restrict mat_lus, 
                                   const Index_t* restrict pivots,
                                   std::complex<T>* restrict logdets,
                                   const size_t batch_size,
                                   const size_t COLBS=128)
  {
    return aq.parallel_for(sycl::nd_range<1>{{batch_size*COLBS},{COLBS}}, 
            [=](sycl::nd_item<1> item) {
            const unsigned iw  = item.get_group(0);
            const unsigned tid = item.get_local_id(0);
            const Index_t* restrict pv_iw = pivots+iw*lda;
            const TMAT* restrict    lu_iw = mat_lus+iw*n*lda;

            std::complex<T> val{};
            const unsigned  num_col_blocks = (n + COLBS - 1) / COLBS;
            for (unsigned block_num = 0; block_num < num_col_blocks; block_num++)
            {
              const unsigned i = tid + block_num * COLBS;
              if(i<n)
                val += (pv_iw[i] == i + 1) ? 
                  std::log(std::complex<T>(lu_iw[i*lda+i])) : std::log(std::complex<T>(-lu_iw[i*lda+i]));
            }

            val = sycl::reduce_over_group(item.get_group(),val,{T{},T{}},std::plus<std::complex<T>>());

            if(iw==0)
              logdets[iw] = val;
            });
  }

  template<typename T, typename Index_t>
  inline sycl::event 
  applyW_stageV(sycl::queue& aq, T* restrict UV ,const size_t norb, 
      const Index_t* restrict delay_list, const size_t delay_count)
  {
    constexpr T mone(-1);

    if(delay_count < 16)
      return aq.submit([&](sycl::handler& cgh) {
          cgh.single_task([=](){
              for(int i=0; i<delay_count; ++i)
              UV[delay_list[i]*delay_count+i] += mone;
              }); 
          });
    else
      return aq.parallel_for(sycl::nd_range<1>{{delay_count},{delay_count}}, 
          [=](sycl::nd_item<1> item) {
          const unsigned i=item.get_global_id(0);
          UV[delay_list[i]*delay_count+i] += mone;
          });
  }

  template<typename T, typename Index_t>
  inline sycl::event 
  applyW_batched(sycl::queue& aq, 
                 const Index_t** restrict delay_list,
                 const int norb, 
                 T** restrict UV,
                 const int delay_count, 
                 const size_t batch_count)
  {
    constexpr T mone(-1);

    return aq.parallel_for(sycl::nd_range<1>{{batch_count},{batch_count}}, 
        [=](sycl::nd_item<1> item) {
        const int  iw                     = item.get_local_id(0);
        const int* restrict delay_list_iw = delay_list[iw];
        T* restrict uv_iw=UV[iw];
        for(int j=0; j<delay_count; j++)
        {
          uv_iw[delay_list_iw[j]*delay_count+j] += mone;
        }
        });
  }

template<typename T>
inline sycl::event copyAinvRow_saveGL(sycl::queue& aq,
                               const int rowchanged,
                               const int n,
                               const T* const Ainv[],
                               const int lda,
                               T* const temp[],
                               T* const rcopy[],
                               const T* const dphi_in[],
                               const T* const d2phi_in[],
                               T* const dphi_out[],
                               T* const d2phi_out[],
                               const size_t batch_count,
                               const size_t COLBS=128)
{
  return aq.parallel_for(sycl::nd_range<1>{{batch_count*COLBS},{COLBS}}, 
      [=](sycl::nd_item<1> item) {
      const int iw                      = item.get_group(0);
      const T* __restrict__ Ainv_iw     = Ainv[iw];
      T* __restrict__ temp_iw           = temp[iw];
      T* __restrict__ rcopy_iw          = rcopy[iw];
      const T* __restrict__ dphi_in_iw  = dphi_in[iw];
      const T* __restrict__ d2phi_in_iw = d2phi_in[iw];
      T* __restrict__ dphi_out_iw       = dphi_out[iw];
      T* __restrict__ d2phi_out_iw      = d2phi_out[iw];

      const unsigned tid = item.get_local_id(0);
      if (tid == 0)
        temp_iw[rowchanged] = -temp_iw[rowchanged];

      const unsigned num_col_blocks = (n + COLBS - 1) / COLBS;
      for (unsigned ib = 0; ib < num_col_blocks; ib++)
      {
        const unsigned col_id = ib * COLBS + tid; 
        if (col_id < n)
        {
          rcopy_iw[col_id] = Ainv_iw[rowchanged * lda + col_id];

          // the following copying data on the device is not part of SM-1
          // it is intended to copy dphiV and d2phiV from temporary to final without a separate kernel.
          dphi_out_iw[col_id * 3]     = dphi_in_iw[col_id * 3];
          dphi_out_iw[col_id * 3 + 1] = dphi_in_iw[col_id * 3 + 1];
          dphi_out_iw[col_id * 3 + 2] = dphi_in_iw[col_id * 3 + 2];
          d2phi_out_iw[col_id]        = d2phi_in_iw[col_id];
        }
      }
      });
}

template<typename T, unsigned DIM=3>
inline sycl::event calcGradients(sycl::queue& aq,
                                 const int n,
                                 const T* const Ainvrow[],
                                 const T* const dpsiMrow[],
                                 T* const grads_now,
                                 size_t batch_count,
                                 const size_t COLBS=128)
{
  //translation of calcGradients_cuda+calcGradients_kernel
  return aq.submit([&](sycl::handler& cgh) {

      sycl::accessor<T, 1, sycl::access::mode::write, sycl::access::target::local> 
         sum(sycl::range<1>{DIM*COLBS}, cgh);
      cgh.parallel_for(sycl::nd_range<1>{{batch_count*COLBS},{COLBS}}, 
      [=](sycl::nd_item<1> item) {
      const unsigned iw                    = item.get_group(0);
      const T* __restrict__ invRow    = Ainvrow[iw];
      const T* __restrict__ dpsiM_row = dpsiMrow[iw];

      const unsigned tid = item.get_local_id(0);
      for (unsigned idim = 0; idim < DIM; idim++)
        sum[idim * COLBS + tid] = T{};

      const unsigned num_col_blocks = (n + COLBS - 1) / COLBS;
      for (unsigned ib = 0; ib < num_col_blocks; ib++)
      {
         const unsigned col_id = ib * COLBS + tid;
         if (col_id < n)
           for (unsigned idim = 0; idim < DIM; idim++)
             sum[idim * COLBS + tid] += invRow[col_id] * dpsiM_row[col_id * DIM + idim];
      }

       for (unsigned iend = COLBS / 2; iend > 0; iend /= 2)
       {
         item.barrier(sycl::access::fence_space::local_space);
         if (tid < iend)
           for (unsigned idim = 0; idim < DIM; idim++)
             sum[idim * COLBS + tid] += sum[idim * COLBS + tid + iend];
       }
     
       if (tid == 0)
         for (int idim = 0; idim < DIM; idim++)
           grads_now[iw * DIM + idim] = sum[idim * COLBS];
      });
  });
}

template<typename T>
sycl::event add_delay_list_save_sigma_VGL(sycl::queue& aq,
                                          int* const delay_list[],
                                          const int rowchanged,
                                          const int delay_count,
                                          T* const binv[],
                                          const int binv_lda,
                                          const T* const ratio_inv,
                                          const T* const phi_in[],
                                          const T* const dphi_in[],
                                          const T* const d2phi_in[],
                                          T* const phi_out[],
                                          T* const dphi_out[],
                                          T* const d2phi_out[],
                                          const int norb,
                                          const int n_accepted,
                                          const size_t batch_count,
                                          const size_t COLBS=64) 
{
  return aq.parallel_for(sycl::nd_range<1>{{batch_count*COLBS},{COLBS}}, 
      [=](sycl::nd_item<1> item) {
  const unsigned tid = item.get_local_id(0);
  const unsigned iw  = item.get_group(0);

  if (iw < n_accepted)
  {
    // real accept, settle y and Z
    int* __restrict__ delay_list_iw   = delay_list[iw];
    T* __restrict__ binvrow_iw        = binv[iw] + delay_count * binv_lda;
    const T* __restrict__ phi_in_iw   = phi_in[iw];
    const T* __restrict__ dphi_in_iw  = dphi_in[iw];
    const T* __restrict__ d2phi_in_iw = d2phi_in[iw];
    T* __restrict__ phi_out_iw        = phi_out[iw];
    T* __restrict__ dphi_out_iw       = dphi_out[iw];
    T* __restrict__ d2phi_out_iw      = d2phi_out[iw];

    if (tid == 0)
    {
      delay_list_iw[delay_count] = rowchanged;
      binvrow_iw[delay_count]    = ratio_inv[iw];
    }

    const unsigned num_delay_count_col_blocks = (delay_count + COLBS - 1) / COLBS;
    for (unsigned ib = 0; ib < num_delay_count_col_blocks; ib++)
    {
      const unsigned col_id = ib * COLBS + tid;
      if (col_id < delay_count)
        binvrow_iw[col_id] *= ratio_inv[iw];
    }

    const unsigned num_col_blocks = (norb + COLBS - 1) / COLBS;
    for (unsigned ib = 0; ib < num_col_blocks; ib++)
    {
      const unsigned col_id = ib * COLBS + tid;
      if (col_id < norb)
      {
        // copy phiV, dphiV and d2phiV from temporary to final without a separate kernel.
        phi_out_iw[col_id]          = phi_in_iw[col_id];
        dphi_out_iw[col_id * 3]     = dphi_in_iw[col_id * 3];
        dphi_out_iw[col_id * 3 + 1] = dphi_in_iw[col_id * 3 + 1];
        dphi_out_iw[col_id * 3 + 2] = dphi_in_iw[col_id * 3 + 2];
        d2phi_out_iw[col_id]        = d2phi_in_iw[col_id];
      }
    }
  }
  else
  {
    // fake accept. Set Y, Z with zero and x with 1
    T* __restrict__ Urow_iw   = phi_out[iw];
    const unsigned num_blocks_norb = (norb + COLBS - 1) / COLBS;
    for (unsigned ib = 0; ib < num_blocks_norb; ib++)
    {
      const unsigned col_id = ib * COLBS + tid;
      if (col_id < norb)
        Urow_iw[col_id] = T(0);
    }

    T* __restrict__ binv_iw          = binv[iw];
    const unsigned num_blocks_delay_count = (delay_count + COLBS - 1) / COLBS;
    for (unsigned ib = 0; ib < num_blocks_delay_count; ib++)
    {
      const unsigned col_id = ib * COLBS + tid;
      if (col_id < delay_count)
        binv_iw[delay_count * binv_lda + col_id] = binv_iw[delay_count + binv_lda * col_id] = T(0);
    }

    int* __restrict__ delay_list_iw = delay_list[iw];
    if (tid == 0)
    {
      binv_iw[delay_count * binv_lda + delay_count] = T(1);
      delay_list_iw[delay_count]                    = -1;
    }
  }
      });
}


  /** utilities for debugging */
  inline double inverse_gflops(size_t N, double t)
  {
    double dn=N;
    return 2.0e-9 * (4./3.)*dn*dn*dn/t;
  }

  template<typename T>
  sycl::event applyW_stageV_sycl(sycl::queue& aq,
                                 const int* restrict delay_list_gpu, const int delay_count,
                                 T* restrict temp_gpu, const int numorbs, const int ndelay,
                                 T* restrict V_gpu, const T* restrict Ainv)
  {
    const size_t BS = 128;
    const size_t NB = (numorbs+BS-1)/BS;

    return aq.parallel_for(sycl::nd_range<1>{{BS*NB},{BS}}, 
        [=](sycl::nd_item<1> item) {
        int col = item.get_global_id(0);     

        // move rows of Ainv to V
        for(int row=0; row<delay_count; row++)
        {
          const T* Ainv_row = Ainv + numorbs * delay_list_gpu[row];
          T* V_row = V_gpu + numorbs * row;
          if( col<numorbs ) V_row[col] = Ainv_row[col];
        }

        // apply W to temp
        if( col<delay_count )
           temp_gpu[ndelay*delay_list_gpu[col] + col] -= T(1);
        });
  }


#if 0
  template<>
    extern sycl::event 
    applyW_stageV_sycl(sycl::queue& aq,
                       const int* restrict delay_list_gpu, const int delay_count,
                       float* restrict temp_gpu, const int numorbs, const int ndelay,
                       float* restrict V_gpu, const float* restrict Ainv);

  template<>
    extern sycl::event 
    applyW_stageV_sycl(sycl::queue& aq,
                       const int* restrict delay_list_gpu, const int delay_count,
                       double* restrict temp_gpu, const int numorbs, const int ndelay,
                       double* restrict V_gpu, const double* restrict Ainv);

  template<>
    extern sycl::event 
    applyW_stageV_sycl(sycl::queue& aq,
                       const int* restrict delay_list_gpu, const int delay_count,
                       std::complex<float>* restrict temp_gpu, const int numorbs, const int ndelay,
                       std::complex<float>* restrict V_gpu, const std::complex<float>* restrict Ainv);

  template<>
    extern sycl::event 
    applyW_stageV_sycl(sycl::queue& aq,
                       const int* restrict delay_list_gpu, const int delay_count,
                       std::complex<double>* restrict temp_gpu, const int numorbs, const int ndelay,
                       std::complex<double>* restrict V_gpu, const std::complex<double>* restrict Ainv);
#endif
}
#endif
