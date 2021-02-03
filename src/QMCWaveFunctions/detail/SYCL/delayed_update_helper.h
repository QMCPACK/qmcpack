//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Thomas Applencourt, apl@anl.gov, Argonne National Laboratory
//
// File created by: Thomas Applencourt, apl@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef SYCL_DELAYED_UPDATE_HELPER_H
#define SYCL_DELAYED_UPDATE_HELPER_H

#include <CL/sycl.hpp>
#include <complex>
namespace sycl = cl::sycl;

#include <iostream>

/** helper kernel for delayed update algorithm
 * W matrix is applied and copy selected rows of Ainv into V
 */
template<typename T, int BS>
void applyW_stageV_kernel(const int* delay_list_gpu,
                          const int delay_count,
                          T* temp_gpu,
                          const int numorbs,
                          const int ndelay,
                          T* V_gpu,
                          const T* Ainv,
                          sycl::nd_item<1> item_ct1)
{
  int col = item_ct1.get_local_id(0) + item_ct1.get_group(0) * BS;
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
    temp_gpu[ndelay * delay_list_gpu[col] + col] = (temp_gpu[ndelay * delay_list_gpu[col] + col]) - T(1);
}

template<typename T>
void applyW_stageV_sycl(const int* delay_list_gpu,
                        const int delay_count,
                        T* temp_gpu,
                        const int numorbs,
                        const int ndelay,
                        T* V_gpu,
                        const T* Ainv,
                        sycl::queue q)
{
  const int BS = 128;
  const int NB = (numorbs + BS - 1) / BS;
  sycl::range<1> dimBlock(BS);
  sycl::range<1> dimGrid(NB);


  q.submit([&](sycl::handler& cgh) {
    cgh.parallel_for(sycl::nd_range<1>(dimBlock * dimGrid, dimBlock), [=](sycl::nd_item<1> item_ct1) {
      applyW_stageV_kernel<T, BS>(delay_list_gpu, delay_count, temp_gpu, numorbs, ndelay, V_gpu, Ainv, item_ct1);
    });
  });
}

/** Identidy matrix
 */
template<typename T, int BS>
void make_identity_matrix_kernel(const int nrows, T* mat, const int lda, sycl::nd_item<2> item_ct1)
{
  int col = item_ct1.get_local_id(0) + item_ct1.get_group(0) * BS;
  if (col < nrows)
  {
    for (int row = item_ct1.get_group(1) * BS;
         row < sycl::min((unsigned int)((item_ct1.get_group(1) + 1) * BS), (unsigned int)nrows); row++)
      mat[row * lda + col] = T(0);
    if (item_ct1.get_group(0) == item_ct1.get_group(1))
      mat[col * lda + col] = T(1);
  }
}

template<typename T>
void make_identity_matrix_sycl(const int nrows, T* mat, const int lda, sycl::queue q)
{
  const int BS = 128;
  const int NB = (nrows + BS - 1) / BS;
  sycl::range<2> dimBlock(BS, 1);
  sycl::range<2> dimGrid(NB, NB);

  q.submit([&](sycl::handler& cgh) {
    cgh.parallel_for(sycl::nd_range<2>(dimGrid * dimBlock, dimBlock),
                     [=](sycl::nd_item<2> item_ct1) { make_identity_matrix_kernel<T, BS>(nrows, mat, lda, item_ct1); });
  });
}

/** extract matrix diagonal
 */
template<typename T, int BS>
void extract_matrix_diagonal_kernel(const int nrows, const T* mat, const int lda, T* diag, sycl::nd_item<1> item_ct1)
{
  int col = item_ct1.get_local_id(0) + item_ct1.get_group(0) * BS;
  if (col < nrows)
    diag[col] = mat[col * lda + col];
}

template<typename T>
void extract_matrix_diagonal_sycl(const int nrows, const T* mat, const int lda, T* diag, sycl::queue q)
{
  const int BS = 128;
  const int NB = (nrows + BS - 1) / BS;
  sycl::range<1> dimBlock(BS);
  sycl::range<1> dimGrid(NB);

  q.submit([&](sycl::handler& cgh) {
    cgh.parallel_for(sycl::nd_range<1>(dimGrid * dimBlock, dimBlock), [=](sycl::nd_item<1> item_ct1) {
      extract_matrix_diagonal_kernel<T, BS>(nrows, mat, lda, diag, item_ct1);
    });
  });
}

/** copy matrix with precision difference
 */

template<typename T_IN, typename T_OUT, int BS>
void copy_matrix_kernel(const int nrows,
                        const int ncols,
                        const T_IN* mat_in,
                        const int lda,
                        T_OUT* mat_out,
                        const int ldb,
                        sycl::nd_item<2> item_ct1)
{
  int col = item_ct1.get_local_id(0) + item_ct1.get_group(0) * BS;
  if (col < ncols)
  {
    for (int row = item_ct1.get_group(1) * BS;
         row < sycl::min((unsigned int)((item_ct1.get_group(1) + 1) * BS), (unsigned int)nrows); row++)
      mat_out[row * ldb + col] = (T_OUT)mat_in[row * lda + col];
  }
}

template<typename T_IN, typename T_OUT>
void copy_matrix_sycl(int nrows, int ncols, T_IN* mat_in, int lda, T_OUT* mat_out, int ldb, sycl::queue q)
{
  const int BS  = 128;
  const int NB1 = (ncols + BS - 1) / BS;
  const int NB2 = (nrows + BS - 1) / BS;
  sycl::range<2> dimBlock(BS, 1);
  sycl::range<2> dimGrid(NB1, NB2);

  q.submit([&](sycl::handler& cgh) {
    cgh.parallel_for(sycl::nd_range<2>(dimGrid * dimBlock, dimBlock), [=](sycl::nd_item<2> item_ct1) {
      copy_matrix_kernel<T_IN, T_OUT, BS>(nrows, ncols, mat_in, lda, mat_out, ldb, item_ct1);
    });
  });
}

#endif
