///////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Fionn Malone, malone14@llnl.gov, LLNL
//
// File created by: Fionn Malone, malone14@llnl.gov, LLNL
////////////////////////////////////////////////////////////////////////////////

#ifndef KERNELS_TERM_BY_TERM_OPERATIONS_H
#define KERNELS_TERM_BY_TERM_OPERATIONS_H

#include <cassert>
#include <complex>

namespace kernels
{
void term_by_term_mat_vec_plus(int dim,
                               int nrow,
                               int ncol,
                               std::complex<double>* A,
                               int lda,
                               std::complex<double>* x,
                               int incx);
void term_by_term_mat_vec_minus(int dim,
                                int nrow,
                                int ncol,
                                std::complex<double>* A,
                                int lda,
                                std::complex<double>* x,
                                int incx);
void term_by_term_mat_vec_mult(int dim,
                               int nrow,
                               int ncol,
                               std::complex<double>* A,
                               int lda,
                               std::complex<double>* x,
                               int incx);
void term_by_term_mat_vec_div(int dim,
                              int nrow,
                              int ncol,
                              std::complex<double>* A,
                              int lda,
                              std::complex<double>* x,
                              int incx);
void term_by_term_mat_vec_plus(int dim, int nrow, int ncol, std::complex<double>* A, int lda, double* x, int incx);
void term_by_term_mat_vec_minus(int dim, int nrow, int ncol, std::complex<double>* A, int lda, double* x, int incx);
void term_by_term_mat_vec_mult(int dim, int nrow, int ncol, std::complex<double>* A, int lda, double* x, int incx);
void term_by_term_mat_vec_div(int dim, int nrow, int ncol, std::complex<double>* A, int lda, double* x, int incx);

} // namespace kernels

#endif
