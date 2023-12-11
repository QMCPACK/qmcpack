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

#ifndef AFQMC_DETERMINANT_KERNELS_HPP
#define AFQMC_DETERMINANT_KERNELS_HPP

#include <cassert>
#include <complex>

namespace kernels
{
double determinant_from_getrf_gpu(int N, double* m, int lda, int* piv, double LogOverlapFactor);
std::complex<double> determinant_from_getrf_gpu(int N,
                                                std::complex<double>* m,
                                                int lda,
                                                int* piv,
                                                std::complex<double> LogOverlapFactor);

void determinant_from_getrf_gpu(int N, double* m, int lda, int* piv, double LogOverlapFactor, double* res);
void determinant_from_getrf_gpu(int N,
                                std::complex<double>* m,
                                int lda,
                                int* piv,
                                std::complex<double> LogOverlapFactor,
                                std::complex<double>* res);

void strided_determinant_from_getrf_gpu(int N,
                                        double* m,
                                        int lda,
                                        int mstride,
                                        int* piv,
                                        int pstride,
                                        double LogOverlapFactor,
                                        double* res,
                                        int nbatch);
void strided_determinant_from_getrf_gpu(int N,
                                        std::complex<double>* m,
                                        int lda,
                                        int mstride,
                                        int* piv,
                                        int pstride,
                                        std::complex<double> LogOverlapFactor,
                                        std::complex<double>* res,
                                        int nbatch);

void batched_determinant_from_getrf_gpu(int N,
                                        double** m,
                                        int lda,
                                        int* piv,
                                        int pstride,
                                        double LogOverlapFactor,
                                        double* res,
                                        int nbatch);
void batched_determinant_from_getrf_gpu(int N,
                                        std::complex<double>** m,
                                        int lda,
                                        int* piv,
                                        int pstride,
                                        std::complex<double> LogOverlapFactor,
                                        std::complex<double>* res,
                                        int nbatch);

std::complex<double> determinant_from_geqrf_gpu(int N, double* m, int lda, double* piv, double LogOverlapFactor);
std::complex<double> determinant_from_geqrf_gpu(int N,
                                                std::complex<double>* m,
                                                int lda,
                                                std::complex<double>* piv,
                                                std::complex<double> LogOverlapFactor);

void determinant_from_geqrf_gpu(int N, double* m, int lda, double* piv);
void determinant_from_geqrf_gpu(int N, std::complex<double>* m, int lda, std::complex<double>* piv);

void scale_columns(int n, int m, double* A, int lda, double* scl);
void scale_columns(int n, int m, std::complex<double>* A, int lda, std::complex<double>* scl);


} // namespace kernels

#endif
