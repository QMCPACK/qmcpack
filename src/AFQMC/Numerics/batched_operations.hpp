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

#ifndef AFQMC_NUMERICS_HELPERS_BATCHED_OPERATIONS_HPP
#define AFQMC_NUMERICS_HELPERS_BATCHED_OPERATIONS_HPP

#include <cassert>
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
#include "AFQMC/Memory/custom_pointers.hpp"
#include "AFQMC/Numerics/device_kernels.hpp"
#endif

namespace ma
{
/** Contract 3D tensor over for use in exchange calculation.
 * E_x[w] ~ sum_{labn} T1[l,w,a,b,n] T2[l,w,b,a,n]
 * l is the batching parameter (like Q vectors in k-point code)
 *
 * \param[in]     nbatch Number of tensors to batch over.
 * \param[in]     nwalk  Number of walkers.
 * \param[in]     nocc   Number of electrons.
 * \param[in]     nchol  Number of Cholesky vectors.
 * \param[in]     alpha  Pointer to array of nbatch scaling factors.
 * \param[in]     Tab    Pointer to packed tensors {T1[l0,w0,a,b,n],T2[[l0,w0,a,b,n]...}.
 *                       Should be allocated on device.
 * \param[in,out] y      Pointer to accumulator array. Typically E[w].
 * \param[in]     incy   Stride for y.
*/
template<typename T, typename Q>
void batched_dot_wabn_wban(int nbatch,
                           int nwalk,
                           int nocc,
                           int nchol,
                           std::complex<Q> const* alpha,
                           std::complex<Q> const* Tab,
                           std::complex<T>* y,
                           int incy)
{
  int nocc2nc = nocc * nocc * nchol;
  for (int batch = 0; batch < nbatch; ++batch)
  {
    for (int w = 0; w < nwalk; ++w)
    {
      std::complex<Q> E_(0.0);
      auto A_(Tab + (2 * batch * nwalk + w) * nocc2nc);
      auto B_(Tab + ((2 * batch + 1) * nwalk + w) * nocc2nc);
      using ma::dot;
      for (int a = 0; a < nocc; ++a)
        for (int b = 0; b < nocc; ++b)
          E_ += ma::dot(nchol, A_ + (a * nocc + b) * nchol, 1, B_ + (b * nocc + a) * nchol, 1);
      y[w * incy] += static_cast<std::complex<T>>(alpha[batch] * E_);
    }
  }
}

template<typename T, typename Q>
void batched_dot_wanb_wbna(int nbatch,
                           int nwalk,
                           int nocc,
                           int nchol,
                           std::complex<Q> const* alpha,
                           std::complex<Q> const* Tab,
                           std::complex<T>* y,
                           int incy)
{
  int nocc2nc = nocc * nocc * nchol;
  for (int batch = 0; batch < nbatch; ++batch)
  {
    for (int w = 0; w < nwalk; ++w)
    {
      std::complex<Q> E_(0.0);
      auto A_(Tab + (2 * batch * nwalk + w) * nocc2nc);
      auto B_(Tab + ((2 * batch + 1) * nwalk + w) * nocc2nc);
      using ma::dot;
      for (int a = 0; a < nocc; ++a)
        for (int b = 0; b < nocc; ++b)
          E_ += ma::dot(nchol, A_ + a * nocc * nchol + b, nocc, B_ + b * nocc * nchol + a, nocc);
      y[w * incy] += static_cast<std::complex<T>>(alpha[batch] * E_);
    }
  }
}

template<typename T, typename Q>
void dot_wabn(int nwalk,
              int nocc,
              int nchol,
              std::complex<Q> alpha,
              std::complex<Q> const* Tab,
              std::complex<T>* y,
              int incy)
{
  int nocc2nc = nocc * nocc * nchol;
  for (int w = 0; w < nwalk; ++w)
  {
    std::complex<Q> E_(0.0);
    auto A_(Tab + w * nocc2nc);
    using ma::dot;
    for (int a = 0; a < nocc; ++a)
      for (int b = 0; b < nocc; ++b)
        E_ += ma::dot(nchol, A_ + (a * nocc + b) * nchol, 1, A_ + (b * nocc + a) * nchol, 1);
    y[w * incy] += static_cast<std::complex<T>>(alpha * E_);
  }
}

/** Construct generalized Fock matrix.
 *
 * \param[in]     nwalk  Number of walkers.
 * \param[in]     nmo    Number of basis functions.
 * \param[in]     nchol  Number of Cholesky vectors.
 * \param[in]     alpha  Scale factor.
 * \param[in]     Tab    Pointer to packed tensors {T1[w0,a,b,n],T2[[w0,a,b,n]...}.
 *                       Should be allocated on device.
 * \param[out]    F      Pointer to buffer for generalised Fock matrix.
*/
template<typename T, typename Q>
void dot_wpan_waqn_Fwpq(int nwalk,
                        int nmo,
                        int nchol,
                        std::complex<Q> alpha,
                        std::complex<Q> const* Tab,
                        std::complex<T>* F)
{
  std::complex<T> alp(static_cast<std::complex<T>>(alpha));
  for (int w = 0; w < nwalk; ++w)
  {
    auto A_(Tab + w * nmo * nmo * nchol);
    using ma::dot;
    for (int p = 0; p < nmo; ++p)
      for (int q = 0; q < nmo; ++q, ++F)
        for (int a = 0; a < nmo; ++a)
          *F += alp *
              static_cast<std::complex<T>>(
                    ma::dot(nchol, A_ + (p * nmo + a) * nchol, 1, A_ + (a * nmo + q) * nchol, 1));
  }
}

template<typename T, typename Q>
void dot_wanb(int nwalk,
              int nocc,
              int nchol,
              std::complex<Q> alpha,
              std::complex<Q> const* Tab,
              std::complex<T>* y,
              int incy)
{
  int nocc2nc = nocc * nchol * nocc;
  for (int w = 0; w < nwalk; ++w)
  {
    std::complex<Q> E_(0.0);
    auto A_(Tab + w * nocc2nc);
    using ma::dot;
    for (int a = 0; a < nocc; ++a)
      for (int b = 0; b < nocc; ++b)
        E_ += ma::dot(nchol, A_ + a * nocc * nchol + b, nocc, A_ + b * nocc * nchol + a, nocc);
    y[w * incy] += static_cast<std::complex<T>>(alpha * E_);
  }
}


/** .
 *
 * \param[in]     nters      Number terms batching over.
 * \param[in]     nwalk      Number of walkers.
 * \param[in]     nocc       Number of electrons.
 * \param[in]     nchol_max  Max number of Cholesky vectors per kpoint.
 * \param[in]     nchol_tot  Total number of Cholesky.
 * \param[in]     ncholQ     Number of Cholesky for Q vector.
 * \param[in]     ncholQ0    Number of Cholesky for Q vector.
 * \param[in]     kdiag      Pointer to array of number of k point pairs for each batch.
 * \param[in]     Tab        Pointer to buffer containing Tab, packed.
 * \param[out]    Kl         Pointer to buffer containing Kl.
 * \param[out]    Kr         Pointer to buffer containing Kr.
*/
template<typename T, typename Q>
void batched_Tab_to_Klr(int nterms,
                        int nwalk,
                        int nocc,
                        int nchol_max,
                        int nchol_tot,
                        int ncholQ,
                        int ncholQ0,
                        int* kdiag,
                        Q const* Tab,
                        T* Kl,
                        T* Kr)
{
  for (int w = 0; w < nwalk; ++w)
  {
    for (int k = 0; k < nterms; k++)
    {
      int batch = kdiag[k];
      for (int a = 0; a < nocc; a++)
      {
        auto Tba_(Tab + batch * nwalk * nocc * nocc * nchol_max + ((w * nocc + a) * nocc + a) * nchol_max);
        auto Kr_(Kr + w * nchol_tot + ncholQ0);
        for (int c = 0; c < ncholQ; ++c)
          Kr_[c] += static_cast<T>(Tba_[c]);
      }
    }
    for (int k = 0; k < nterms; k++)
    {
      int batch = kdiag[k];
      for (int a = 0; a < nocc; a++)
      {
        auto Tab_(Tab + (batch + 1) * nwalk * nocc * nocc * nchol_max + ((w * nocc + a) * nocc + a) * nchol_max);
        auto Kl_(Kl + w * nchol_tot + ncholQ0);
        for (int c = 0; c < ncholQ; ++c)
          Kl_[c] += static_cast<T>(Tab_[c]);
      }
    }
  }
}

// Not used.
template<typename T, typename Q>
void batched_Tanb_to_Klr(int nterms,
                         int nwalk,
                         int nocc,
                         int nchol_max,
                         int nchol_tot,
                         int ncholQ,
                         int ncholQ0,
                         int* kdiag,
                         Q const* Tab,
                         T* Kl,
                         T* Kr)
{
  for (int w = 0; w < nwalk; ++w)
  {
    for (int k = 0; k < nterms; k++)
    {
      int batch = kdiag[k];
      for (int a = 0; a < nocc; a++)
      {
        auto Tba_(Tab + batch * nwalk * nocc * nocc * nchol_max + ((w * nocc + a) * nocc) * nchol_max + a);
        auto Kr_(Kr + w * nchol_tot + ncholQ0);
        for (int c = 0; c < ncholQ; ++c)
          Kr_[c] += static_cast<T>(Tba_[c * nocc]);
      }
    }
    for (int k = 0; k < nterms; k++)
    {
      int batch = kdiag[k];
      for (int a = 0; a < nocc; a++)
      {
        auto Tab_(Tab + (batch + 1) * nwalk * nocc * nocc * nchol_max + ((w * nocc + a) * nocc) * nchol_max + a);
        auto Kl_(Kl + w * nchol_tot + ncholQ0);
        for (int c = 0; c < ncholQ; ++c)
          Kl_[c] += static_cast<T>(Tab_[c * nocc]);
      }
    }
  }
}

template<typename T, typename Q>
void Tab_to_Kl(int nwalk, int nocc, int nchol, Q const* Tab, T* Kl)
{
  for (int w = 0; w < nwalk; ++w)
  {
    for (int a = 0; a < nocc; a++)
    {
      auto Tab_(Tab + ((w * nocc + a) * nocc + a) * nchol);
      auto Kl_(Kl + w * nchol);
      for (int c = 0; c < nchol; ++c)
        Kl_[c] += static_cast<T>(Tab_[c]);
    }
  }
}

template<typename T, typename Q>
void Tanb_to_Kl(int nwalk, int nocc, int nchol, int nchol_tot, Q const* Tab, T* Kl)
{
  for (int w = 0; w < nwalk; ++w)
  {
    for (int a = 0; a < nocc; a++)
    {
      auto Tab_(Tab + ((w * nocc + a) * nocc) * nchol + a);
      auto Kl_(Kl + w * nchol_tot);
      for (int c = 0; c < nchol; ++c)
        Kl_[c] += static_cast<T>(Tab_[c * nocc]);
    }
  }
}

template<typename T, typename T1>
void vbias_from_v1(int nwalk,
                   int nkpts,
                   int nchol_max,
                   int* Qsym,
                   int* kminus,
                   int* ncholpQ,
                   int* ncholpQ0,
                   std::complex<T> const alpha,
                   std::complex<T1> const* v1,
                   std::complex<T>* vb)
{
  for (int Q = 0; Q < nkpts; Q++)
  {
    if (Qsym[Q] < 0)
      return;
    int Qm   = kminus[Q];
    int nc0  = ncholpQ0[Q];
    int nc   = ncholpQ[Q];
    int ncm  = ncholpQ[Qm];
    int Qm_  = Qm;
    int ntot = nc * nwalk;
    if (Qsym[Q] > 0)
      Qm_ = nkpts + Qsym[Q] - 1;

    // v+
    auto vb_(vb + nc0 * nwalk);
    auto v1_(v1 + Q * nchol_max * nwalk);
    auto v2_(v1 + Qm_ * nchol_max * nwalk);
    // v+ = a*(v[Q]+v[-Q])
    for (int n = 0; n < ntot; ++n)
      vb_[n] += alpha * static_cast<std::complex<T>>(v1_[n]);
    for (int n = 0; n < ntot; ++n)
      vb_[n] += alpha * static_cast<std::complex<T>>(v2_[n]);
    // v-
    vb_ = (vb + (nc0 + nc) * nwalk);
    // v- = -a*i*(v[Q]-v[-Q])
    auto ialpha(alpha * std::complex<T>(0.0, 1.0));
    for (int n = 0; n < ntot; ++n)
      vb_[n] -= ialpha * static_cast<std::complex<T>>(v1_[n]);
    for (int n = 0; n < ntot; ++n)
      vb_[n] += ialpha * static_cast<std::complex<T>>(v2_[n]);
  }
}

// for n in [0,N), y[incy*n] = beta * y[incy*n] + alpha sum_m^{0,M} opA(A)[n,m] * opB(B)[n,m]
template<typename T, typename Q>
void batched_dot(char TA,
                 char TB,
                 int N,
                 int M,
                 std::complex<T> const alpha,
                 std::complex<Q> const* A,
                 int lda,
                 std::complex<Q> const* B,
                 int ldb,
                 std::complex<T> const beta,
                 std::complex<T>* y,
                 int incy)
{
  bool cA(TA == 'H' || TA == 'C');
  bool cB(TB == 'H' || TB == 'C');
  bool tA(TA == 'H' || TA == 'T');
  bool tB(TB == 'H' || TB == 'T');
  if (not tA && not TB)
  {
    for (int n = 0; n < N; n++)
    {
      std::complex<T> r(0.0, 0.0);
      auto an(A + n * lda);
      auto bn(B + n * ldb);
      if (cA && cB)
      {
        for (int m = 0; m < M; m++, an++, bn++)
          r += std::conj(*an) * std::conj(*bn);
      }
      else if (cA && not cB)
      {
        for (int m = 0; m < M; m++, an++, bn++)
          r += std::conj(*an) * (*bn);
      }
      else if (not cA && cB)
      {
        for (int m = 0; m < M; m++, an++, bn++)
          r += (*an) * std::conj(*bn);
      }
      else
      {
        for (int m = 0; m < M; m++, an++, bn++)
          r += (*an) * (*bn);
      }
      y[incy * n] = beta * y[incy * n] + alpha * r;
    }
  }
  else if (tA && not TB)
  {
    for (int n = 0; n < N; n++)
    {
      std::complex<T> r(0.0, 0.0);
      auto an(A + n);
      auto bn(B + n * ldb);
      if (cA && cB)
      {
        for (int m = 0; m < M; m++, an += lda, bn++)
          r += std::conj(*an) * std::conj(*bn);
      }
      else if (cA && not cB)
      {
        for (int m = 0; m < M; m++, an += lda, bn++)
          r += std::conj(*an) * (*bn);
      }
      else if (not cA && cB)
      {
        for (int m = 0; m < M; m++, an += lda, bn++)
          r += (*an) * std::conj(*bn);
      }
      else
      {
        for (int m = 0; m < M; m++, an += lda, bn++)
          r += (*an) * (*bn);
      }
      y[incy * n] = beta * y[incy * n] + alpha * r;
    }
  }
  else if (not tA && TB)
  {
    for (int n = 0; n < N; n++)
    {
      std::complex<T> r(0.0, 0.0);
      auto an(A + n * lda);
      auto bn(B + n);
      if (cA && cB)
      {
        for (int m = 0; m < M; m++, an++, bn += ldb)
          r += std::conj(*an) * std::conj(*bn);
      }
      else if (cA && not cB)
      {
        for (int m = 0; m < M; m++, an++, bn += ldb)
          r += std::conj(*an) * (*bn);
      }
      else if (not cA && cB)
      {
        for (int m = 0; m < M; m++, an++, bn += ldb)
          r += (*an) * std::conj(*bn);
      }
      else
      {
        for (int m = 0; m < M; m++, an++, bn += ldb)
          r += (*an) * (*bn);
      }
      y[incy * n] = beta * y[incy * n] + alpha * r;
    }
  }
  else
  { // special case, tA && tB
    for (int n = 0; n < N; n++)
      y[incy * n] *= beta;
    for (int m = 0; m < M; m++)
    {
      auto am(A + m * lda);
      auto bm(B + m * ldb);
      if (cA && cB)
      {
        for (int n = 0; n < N; n++, am++, bm++)
          y[incy * n] += alpha * std::conj(*am) * std::conj(*bm);
      }
      else if (cA && not cB)
      {
        for (int n = 0; n < N; n++, am++, bm++)
          y[incy * n] += alpha * std::conj(*am) * (*bm);
      }
      else if (not cA && cB)
      {
        for (int n = 0; n < N; n++, am++, bm++)
          y[incy * n] += alpha * (*am) * std::conj(*bm);
      }
      else
      {
        for (int n = 0; n < N; n++, am++, bm++)
          y[incy * n] += alpha * (*am) * (*bm);
      }
    }
  }
}

// y[s] = y[s] + sum_ab A[s][a][b] * B[s][b][a]
// shapes of arrays are in packed form in n array
template<typename T, typename Q>
void batched_ab_ba(int* n,
                   std::complex<Q>* const* A,
                   int lda,
                   std::complex<Q>* const* B,
                   int ldb,
                   std::complex<T> alpha,
                   std::complex<T>** y,
                   int batchSize)
{
  // not optimal, blocked algorithm is faster
  for (int b = 0; b < batchSize; b++)
  {
    int n1(n[2 * b]);
    int n2(n[2 * b + 1]);
    std::complex<Q> const* A_(A[b]);
    std::complex<Q> const* B_(B[b]);
    std::complex<T> y_(0.0);
    for (int i = 0; i < n1; i++)
      for (int j = 0; j < n2; j++)
        y_ += static_cast<std::complex<T>>(A_[i * lda + j] * B_[j * ldb + i]);
    *(y[b]) += alpha * y_;
  }
}

template<typename T, typename Q>
void batched_diagonal_sum(int* n,
                          std::complex<Q>* const* A,
                          int lda,
                          std::complex<T> alpha,
                          std::complex<T>** y,
                          int batchSize)
{
  for (int b = 0; b < batchSize; b++)
  {
    std::complex<Q> const* A_(A[b]);
    std::complex<T> y_(0.0);
    int n_(n[b]);
    for (int i = 0; i < n_; i++)
      y_ += static_cast<std::complex<T>>(A_[i * lda + i]);
    *(y[b]) = alpha * y_;
  }
}

// C[u][w] = a * sum_n A[u][w][n] * B[u][n]
// should be called: Aijk_Bik_Cij
template<typename T>
void Auwn_Bun_Cuw(int nu, int nw, int na, T alpha, T const* A, T const* B, T* C)
{
  for (int u = 0; u < nu; ++u, B += na)
    for (int iw = 0; iw < nw; ++iw, ++C, A += na)
      (*C) = alpha * ma::dot(na, A, 1, B, 1);
}

// C[u][w] = alpha * sum_i A[w][i][u] * B[i][u]
template<typename T, typename Q>
// Aijk_Bjk_Cki
void Awiu_Biu_Cuw(int nu, int nw, int ni, T alpha, T const* A, Q const* B, int ldb, T* C, int ldc)
{
  for (int w = 0; w < nw; ++w)
  {
    for (int i = 0; i < ni; ++i)
    {
      auto Ci(C + w);
      auto Bi(B + i * ldb);
      for (int u = 0; u < nu; ++u, ++A, ++Bi, Ci += ldc)
        *Ci += alpha * (*A) * static_cast<T>(*Bi);
    }
  }
}

// C[i][k] = sum_i A[i][j][k] * B[k][j]
template<typename T, typename T1>
void Aijk_Bkj_Cik(int ni, int nj, int nk, T const* A, int lda, int stride, T1 const* B, int ldb, T* C, int ldc)
{
  for (int i = 0; i < ni; ++i)
  {
    auto Ci(C + i * ldc);
    auto Ai_(A + i * stride);
    for (int k = 0; k < nk; ++k, ++Ci)
    {
      auto Ak(Ai_ + k);
      auto Bk(B + k * ldb);
      for (int j = 0; j < nj; ++j, Ak += lda, ++Bk)
        *Ci += (*Ak) * static_cast<T>(*Bk);
    }
  }
}

// A[w][i][j] = B[i][w][j]
template<typename T, typename T1>
void viwj_vwij(int nw, int ni, int i0, int iN, T const* B, T1* A)
{
  for (int w = 0; w < nw; ++w)
  {
    for (int i = i0; i < iN; ++i)
    {
      auto A_(A + (w * ni + i) * ni);
      auto B_(B + (i * nw + w) * ni);
      for (int j = 0; j < ni; ++j, ++A_, ++B_)
        *A_ = static_cast<T1>(*B_);
    }
  }
}

// Ckij = transA(Aij) * Bjk
//        conj(Aij)?
template<typename T>
void element_wise_Aij_Bjk_Ckij(char transA,
                               int ni,
                               int nj,
                               int nk,
                               T const* A,
                               int lda,
                               T const* B,
                               int ldb,
                               T* C,
                               int ldc1,
                               int ldc2)
{
  if (transA == 'N')
  {
    for (int k = 0; k < nk; ++k)
    {
      for (int i = 0; i < ni; ++i)
      {
        auto A_(A + i * lda);
        auto B_(B + k);
        auto C_(C + k * ldc1 * ldc2 + i * ldc2);
        for (int j = 0; j < nj; ++j, ++C_, ++A_, B_ += ldb)
          (*C_) = (*A_) * (*B_);
      }
    }
  }
  else if (transA == 'C')
  {
    for (int k = 0; k < nk; ++k)
    {
      for (int i = 0; i < ni; ++i)
      {
        auto A_(A + i * lda);
        auto B_(B + k);
        auto C_(C + k * ldc1 * ldc2 + i * ldc2);
        for (int j = 0; j < nj; ++j, ++C_, ++A_, B_ += ldb)
          (*C_) = ma::conj(*A_) * (*B_);
      }
    }
  }
  else
    throw std::runtime_error(" Invalid parameter in element_wise_Aij_Bjk_Ckij. ");
}

template<typename T1, typename T2>
void element_wise_Aij_Bjk_Ckji(int ni,
                               int nj,
                               int nk,
                               T1 const* A,
                               int lda,
                               T2 const* B,
                               int ldb,
                               T2* C,
                               int ldc,
                               int stride)
{
  for (int k = 0; k < nk; k++)
  {
    for (int j = 0; j < nj; j++)
    {
      auto A_(A + j);
      auto B_(*(B + j * ldb + k));
      auto C_(C + k * stride + j * ldc);
      for (int i = 0; i < ni; i++, A_ += lda, ++C_)
        (*C_) = (*A_) * B_;
    }
  }
}

// A[n][i][j] * = B[i][j]
template<typename T, typename T1>
void inplace_product(int nbatch, int n, int m, T1 const* B, int ldb, std::complex<T>* A, int lda)
{
  for (int nb = 0; nb < nbatch; nb++)
    for (int i = 0; i < n; ++i)
    {
      auto A_(A + (nb * n + i) * lda);
      auto B_(B + i * ldb);
      for (int j = 0; j < m; ++j, ++A_, ++B_)
        *A_ *= static_cast<std::complex<T>>(*B_);
    }
}


} // namespace ma

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
namespace device
{
template<typename T, typename Q>
void batched_Tab_to_Klr(int nterms,
                        int nwalk,
                        int nocc,
                        int nchol_max,
                        int nchol_tot,
                        int ncholQ,
                        int ncholQ0,
                        device_pointer<int> kdiag,
                        device_pointer<Q> Tab,
                        device_pointer<T> Kl,
                        device_pointer<T> Kr)
{
  kernels::batched_Tab_to_Klr(nterms, nwalk, nocc, nchol_max, nchol_tot, ncholQ, ncholQ0, to_address(kdiag),
                              to_address(Tab), to_address(Kl), to_address(Kr));
}

// Not used.
template<typename T, typename Q>
void batched_Tanb_to_Klr(int nterms,
                         int nwalk,
                         int nocc,
                         int nchol_max,
                         int nchol_tot,
                         int ncholQ,
                         int ncholQ0,
                         device_pointer<int> kdiag,
                         device_pointer<Q> Tab,
                         device_pointer<T> Kl,
                         device_pointer<T> Kr)
{
  kernels::batched_Tanb_to_Klr(nterms, nwalk, nocc, nchol_max, nchol_tot, ncholQ, ncholQ0, to_address(kdiag),
                               to_address(Tab), to_address(Kl), to_address(Kr));
}

template<typename T, typename Q>
void Tab_to_Kl(int nwalk, int nocc, int nchol, device_pointer<Q> Tab, device_pointer<T> Kl)
{
  kernels::Tab_to_Kl(nwalk, nocc, nchol, to_address(Tab), to_address(Kl));
}

template<typename T, typename Q>
void Tanb_to_Kl(int nwalk, int nocc, int nchol, int nchol_tot, device_pointer<Q> Tab, device_pointer<T> Kl)
{
  kernels::Tanb_to_Kl(nwalk, nocc, nchol, nchol_tot, to_address(Tab), to_address(Kl));
}

template<typename T, typename Q, typename R>
void batched_dot_wabn_wban(int nbatch,
                           int nwalk,
                           int nocc,
                           int nchol,
                           device_pointer<R> alpha,
                           device_pointer<Q> Tab,
                           T* y,
                           int incy)
{
  kernels::batched_dot_wabn_wban(nbatch, nwalk, nocc, nchol, to_address(alpha), to_address(Tab), y, incy);
}

template<typename T, typename Q, typename R>
void batched_dot_wanb_wbna(int nbatch,
                           int nwalk,
                           int nocc,
                           int nchol,
                           device_pointer<R> alpha,
                           device_pointer<Q> Tab,
                           T* y,
                           int incy)
{
  kernels::batched_dot_wanb_wbna(nbatch, nwalk, nocc, nchol, to_address(alpha), to_address(Tab), y, incy);
}

template<typename T, typename Q, typename R>
void dot_wabn(int nwalk, int nocc, int nchol, R alpha, device_pointer<Q> Tab, T* y, int incy)
{
  kernels::dot_wabn(nwalk, nocc, nchol, alpha, to_address(Tab), y, incy);
}

// Not used.
template<typename T, typename Q, typename R>
void dot_wpan_waqn_Fwpq(int nwalk, int nocc, int nchol, R alpha, device_pointer<Q> Tab, device_pointer<T> F)
{
  kernels::dot_wpan_waqn_Fwpq(nwalk, nocc, nchol, alpha, to_address(Tab), to_address(F));
}

template<typename T, typename Q, typename R>
void dot_wanb(int nwalk, int nocc, int nchol, R alpha, device_pointer<Q> Tab, T* y, int incy)
{
  kernels::dot_wanb(nwalk, nocc, nchol, alpha, to_address(Tab), y, incy);
}

template<typename T, typename Q, typename R>
void vbias_from_v1(int nwalk,
                   int nkpts,
                   int nchol_max,
                   device_pointer<int> Qsym,
                   device_pointer<int> kminus,
                   device_pointer<int> ncholpQ,
                   device_pointer<int> ncholpQ0,
                   R alpha,
                   device_pointer<Q> v1,
                   T* vb)
{
  kernels::vbias_from_v1(nwalk, nkpts, nchol_max, to_address(Qsym), to_address(kminus), to_address(ncholpQ),
                         to_address(ncholpQ0), alpha, to_address(v1), vb);
}

template<typename T1, typename T2, typename T3, typename T4>
void Auwn_Bun_Cuw(int nu, int nw, int na, T1 alpha, device_pointer<T2> A, device_pointer<T3> B, device_pointer<T4> C)
{
  kernels::Auwn_Bun_Cuw(nu, nw, na, alpha, to_address(A), to_address(B), to_address(C));
}

template<typename T1, typename T2, typename T3, typename T4>
void Awiu_Biu_Cuw(int nu,
                  int nw,
                  int ni,
                  T1 alpha,
                  device_pointer<T2> A,
                  device_pointer<T3> B,
                  int ldb,
                  device_pointer<T4> C,
                  int ldc)
{
  kernels::Awiu_Biu_Cuw(nu, nw, ni, alpha, to_address(A), to_address(B), ldb, to_address(C), ldc);
}

template<typename T1, typename T2, typename T3>
void Aijk_Bkj_Cik(int ni,
                  int nj,
                  int nk,
                  device_pointer<T1> A,
                  int lda,
                  int stride,
                  device_pointer<T2> B,
                  int ldb,
                  device_pointer<T3> C,
                  int ldc)
{
  kernels::Aijk_Bkj_Cik(ni, nj, nk, to_address(A), lda, stride, to_address(B), ldb, to_address(C), ldc);
}

// A[w][i][j] = B[i][w][j]
template<typename T, typename T1>
void viwj_vwij(int nw, int ni, int i0, int iN, device_pointer<T> B, device_pointer<T1> A)
{
  kernels::viwj_vwij(nw, ni, i0, iN, to_address(B), to_address(A));
}

template<typename T1, typename T2, typename T3>
void element_wise_Aij_Bjk_Ckij(char transA,
                               int ni,
                               int nj,
                               int nk,
                               device_pointer<T1> A,
                               int lda,
                               device_pointer<T2> B,
                               int ldb,
                               device_pointer<T3> C,
                               int ldc1,
                               int ldc2)
{
  kernels::element_wise_Aij_Bjk_Ckij(transA, ni, nj, nk, to_address(A), lda, to_address(B), ldb, to_address(C), ldc1,
                                     ldc2);
}

template<typename T1, typename T2, typename T3>
void element_wise_Aij_Bjk_Ckji(int ni,
                               int nj,
                               int nk,
                               device_pointer<T1> A,
                               int lda,
                               device_pointer<T2> B,
                               int ldb,
                               device_pointer<T3> C,
                               int ldc,
                               int stride)
{
  kernels::element_wise_Aij_Bjk_Ckji(ni, nj, nk, to_address(A), lda, to_address(B), ldb, to_address(C), ldc, stride);
}

template<typename T, typename T1>
void inplace_product(int nbatch, int n, int m, device_pointer<T1> B, int ldb, device_pointer<T> A, int lda)
{
  kernels::inplace_product(nbatch, n, m, to_address(B), ldb, to_address(A), lda);
}

template<typename T, typename Q>
void batched_dot(char TA,
                 char TB,
                 int N,
                 int M,
                 T alpha,
                 device_pointer<Q> A,
                 int lda,
                 device_pointer<Q> B,
                 int ldb,
                 T beta,
                 device_pointer<T> y,
                 int incy)
{
  //  kernels::batched_dot(TA,TB,N,M,alpha,to_address(A),lda,to_address(B),ldb,beta,to_address(y),incy);
  APP_ABORT(" Error: batched_dot not yet available in gpu.\n");
}

template<typename I, typename T, typename Q, typename T1>
void batched_ab_ba(device_pointer<I> n,
                   device_pointer<Q>* A,
                   int lda,
                   device_pointer<Q>* B,
                   int ldb,
                   T1 alpha,
                   device_pointer<T>* y,
                   int batchSize)
{
  //  kernels::batched_dot(TA,TB,N,M,alpha,to_address(A),lda,to_address(B),ldb,beta,to_address(y),incy);
  APP_ABORT(" Error: batched_ab_ba not yet available in gpu.\n");
}

template<typename I, typename T, typename Q>
void batched_diagonal_sum(device_pointer<I> n,
                          device_pointer<Q>* A,
                          int lda,
                          T alpha,
                          device_pointer<T>* y,
                          int batchSize)
{
  //  kernels::batched_dot(TA,TB,N,M,alpha,to_address(A),lda,to_address(B),ldb,beta,to_address(y),incy);
  APP_ABORT(" Error: batched_diagonal_sum not yet available in gpu.\n");
}

} // namespace device
#endif


#endif
