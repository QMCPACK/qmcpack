//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_AUWN_BUN_CUW_KERNELS_HPP
#define AFQMC_AUWN_BUN_CUW_KERNELS_HPP

#include <complex>

namespace kernels
{
// C[u][w] = alpha * sum_a A[u][w][a] * B[u][a]
void Auwn_Bun_Cuw(int nu,
                  int nw,
                  int na,
                  std::complex<double> alpha,
                  std::complex<double> const* A,
                  std::complex<double> const* B,
                  std::complex<double>* C);
void Auwn_Bun_Cuw(int nu,
                  int nw,
                  int na,
                  std::complex<float> alpha,
                  std::complex<float> const* A,
                  std::complex<float> const* B,
                  std::complex<float>* C);

// C[u][w] = alpha * sum_i A[w][i][u] * B[i][u]
void Awiu_Biu_Cuw(int nu,
                  int nw,
                  int ni,
                  std::complex<double> alpha,
                  std::complex<double> const* A,
                  double const* B,
                  int ldb,
                  std::complex<double>* C,
                  int ldc);
void Awiu_Biu_Cuw(int nu,
                  int nw,
                  int ni,
                  std::complex<float> alpha,
                  std::complex<float> const* A,
                  float const* B,
                  int ldb,
                  std::complex<float>* C,
                  int ldc);
void Awiu_Biu_Cuw(int nu,
                  int nw,
                  int ni,
                  std::complex<double> alpha,
                  std::complex<double> const* A,
                  std::complex<double> const* B,
                  int ldb,
                  std::complex<double>* C,
                  int ldc);
void Awiu_Biu_Cuw(int nu,
                  int nw,
                  int ni,
                  std::complex<float> alpha,
                  std::complex<float> const* A,
                  std::complex<float> const* B,
                  int ldb,
                  std::complex<float>* C,
                  int ldc);

// C[i][k] = sum_i A[i][j][k] * B[k][j]
void Aijk_Bkj_Cik(int ni,
                  int nj,
                  int nk,
                  std::complex<double> const* A,
                  int lda,
                  int stride,
                  std::complex<double> const* B,
                  int ldb,
                  std::complex<double>* C,
                  int ldc);
void Aijk_Bkj_Cik(int ni,
                  int nj,
                  int nk,
                  std::complex<double> const* A,
                  int lda,
                  int stride,
                  double const* B,
                  int ldb,
                  std::complex<double>* C,
                  int ldc);
void Aijk_Bkj_Cik(int ni,
                  int nj,
                  int nk,
                  std::complex<float> const* A,
                  int lda,
                  int stride,
                  std::complex<float> const* B,
                  int ldb,
                  std::complex<float>* C,
                  int ldc);
void Aijk_Bkj_Cik(int ni,
                  int nj,
                  int nk,
                  std::complex<float> const* A,
                  int lda,
                  int stride,
                  float const* B,
                  int ldb,
                  std::complex<float>* C,
                  int ldc);

// A[w][i][j] = B[i][w][j]
void viwj_vwij(int nw, int ni, int i0, int iN, std::complex<double> const* B, std::complex<double>* A);
void viwj_vwij(int nw, int ni, int i0, int iN, std::complex<double> const* B, std::complex<float>* A);
void viwj_vwij(int nw, int ni, int i0, int iN, std::complex<float> const* B, std::complex<double>* A);
void viwj_vwij(int nw, int ni, int i0, int iN, std::complex<float> const* B, std::complex<float>* A);

// element-wise C[k][i][j] = A[i][j] * B[j][k]
void element_wise_Aij_Bjk_Ckij(char transA,
                               int ni,
                               int nj,
                               int nk,
                               double const* A,
                               int lda,
                               std::complex<double> const* B,
                               int ldb,
                               std::complex<double>* C,
                               int ldc1,
                               int ldc2);
void element_wise_Aij_Bjk_Ckij(char transA,
                               int ni,
                               int nj,
                               int nk,
                               float const* A,
                               int lda,
                               std::complex<float> const* B,
                               int ldb,
                               std::complex<float>* C,
                               int ldc1,
                               int ldc2);
void element_wise_Aij_Bjk_Ckij(char transA,
                               int ni,
                               int nj,
                               int nk,
                               std::complex<double> const* A,
                               int lda,
                               std::complex<double> const* B,
                               int ldb,
                               std::complex<double>* C,
                               int ldc1,
                               int ldc2);
void element_wise_Aij_Bjk_Ckij(char transA,
                               int ni,
                               int nj,
                               int nk,
                               std::complex<float> const* A,
                               int lda,
                               std::complex<float> const* B,
                               int ldb,
                               std::complex<float>* C,
                               int ldc1,
                               int ldc2);

// element-wise C[k][j][i] = A[i][j] * B[j][k]
void element_wise_Aij_Bjk_Ckji(int ni,
                               int nj,
                               int nk,
                               double const* A,
                               int lda,
                               std::complex<double> const* B,
                               int ldb,
                               std::complex<double>* C,
                               int ldc,
                               int stride);
void element_wise_Aij_Bjk_Ckji(int ni,
                               int nj,
                               int nk,
                               float const* A,
                               int lda,
                               std::complex<float> const* B,
                               int ldb,
                               std::complex<float>* C,
                               int ldc,
                               int stride);
void element_wise_Aij_Bjk_Ckji(int ni,
                               int nj,
                               int nk,
                               std::complex<double> const* A,
                               int lda,
                               std::complex<double> const* B,
                               int ldb,
                               std::complex<double>* C,
                               int ldc,
                               int stride);
void element_wise_Aij_Bjk_Ckji(int ni,
                               int nj,
                               int nk,
                               std::complex<float> const* A,
                               int lda,
                               std::complex<float> const* B,
                               int ldb,
                               std::complex<float>* C,
                               int ldc,
                               int stride);

} // namespace kernels

#endif
