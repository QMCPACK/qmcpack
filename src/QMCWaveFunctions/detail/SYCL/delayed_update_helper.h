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

#include <complex>

/** helper function for delayed update algorithm
 * W matrix is applied and copy selected rows of Ainv into V
 */

void applyW_stageV_sycl(const int *delay_list_gpu, const int delay_count,
                        float* temp_gpu, const int numorbs, const int ndelay,
                        float* V_gpu, const float* Ainv,
                        sycl::queue q);

void applyW_stageV_sycl(const int *delay_list_gpu, const int delay_count,
                        std::complex<float>* temp_gpu, const int numorbs, const int ndelay,
                        std::complex<float>* V_gpu, const std::complex<float>* Ainv,
                        sycl::queue q);

void applyW_stageV_sycl(const int *delay_list_gpu, const int delay_count,
                        double* temp_gpu, const int numorbs, const int ndelay,
                        double* V_gpu, const double* Ainv,
                        sycl::queue q);

void applyW_stageV_sycl(const int *delay_list_gpu, const int delay_count,
                        std::complex<double>* temp_gpu, const int numorbs, const int ndelay,
                        std::complex<double>* V_gpu, const std::complex<double>* Ainv,
                        sycl::queue q);

/** create identity matrix on the device
 */
void make_identity_matrix_sycl(const int nrows, double* mat, const int lda, sycl::queue q);

void make_identity_matrix_sycl(const int nrows, std::complex<double>* mat, const int lda, sycl::queue q);

/** extract matrix diagonal
 */
void extract_matrix_diagonal_sycl(const int nrows, const double* mat, const int lda, double* diag, sycl::queue q);

void extract_matrix_diagonal_sycl(const int nrows, const std::complex<double>* mat, const int lda, std::complex<double>* diag, sycl::queue q);

/** copy matrix with precision difference
 */
void copy_matrix_sycl(const int nrows, const int ncols, const double* mat_in, const int lda, float* mat_out, const int ldb, sycl::queue q);

void copy_matrix_sycl(const int nrows, const int ncols, const float* mat_in, const int lda, double* mat_out, const int ldb, sycl::queue q);

void copy_matrix_sycl(const int nrows, const int ncols, const std::complex<double>* mat_in, const int lda, std::complex<float>* mat_out, const int ldb, sycl::queue q);

void copy_matrix_sycl(const int nrows, const int ncols, const std::complex<float>* mat_in, const int lda, std::complex<double>* mat_out, const int ldb, sycl::queue q);

#endif
