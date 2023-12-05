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

#ifndef AFQMC_NUMERICS_HELPERS_TENSOR_TRANSPOSITION_HPP
#define AFQMC_NUMERICS_HELPERS_TENSOR_TRANSPOSITION_HPP

#include <cassert>
#include "AFQMC/Numerics/detail/utilities.hpp"
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
#include "AFQMC/Memory/custom_pointers.hpp"
#include "AFQMC/Numerics/device_kernels.hpp"
#endif

namespace ma
{
// Need OpenMP!!!

template<typename Q, typename T>
void KaKjw_to_KKwaj(int nwalk,
                    int nkpts,
                    int npol,
                    int nmo_max,
                    int nmo_tot,
                    int nocc_max,
                    int* nopk,
                    int* nopk0,
                    int* nelpk,
                    int* nelpk0,
                    Q const* A,
                    T* B)
{
  // OpenMP: Combine Ka,Kj loops into single loop and call parallel for
  int napj = nocc_max * npol * nmo_max;
  int na0  = 0;
  for (int Ka = 0; Ka < nkpts; Ka++)
  {
    int na  = nelpk[Ka];
    int nj0 = 0;
    for (int Kj = 0; Kj < nkpts; Kj++)
    {
      int nj = nopk[Kj];
      //auto G_(to_address(GKK[0][Ka][Kj].origin()));
      auto G_(B + (Ka * nkpts + Kj) * nwalk * nocc_max * npol * nmo_max);
      for (int a = 0; a < na; a++)
      {
        //auto Gc_( to_address(Gca[na0+a][p][nj0].origin()) );
        int apj = a * npol * nmo_max;
        for (int p = 0; p < npol; p++)
        {
          auto Gc_(A + ((na0 + a) * npol + p) * nmo_tot * nwalk + nj0 * nwalk);
          for (int j = 0; j < nj; j++, apj++)
          {
            for (int w = 0, wapj = 0; w < nwalk; w++, ++Gc_, wapj += napj)
              G_[wapj + apj] = static_cast<T>(*Gc_);
          }
        }
      }
      nj0 += nj;
    }
    na0 += na;
  }
}

template<typename T, typename T1>
void KaKjw_to_QKajw(int nwalk,
                    int nkpts,
                    int npol,
                    int nmo_max,
                    int nmo_tot,
                    int nocc_max,
                    int* nmo,
                    int* nmo0,
                    int* nocc,
                    int* nocc0,
                    int* QKtok2,
                    T1 const* A,
                    T* B)
{
  // OpenMP: Combine Q,K loops into single loop and call parallel for
  for (int Q = 0; Q < nkpts; Q++)
  {
    for (int K = 0; K < nkpts; K++)
    {
      int Ka  = K;
      int Kj  = QKtok2[Q * nkpts + Ka];
      int na  = nocc[Ka];
      int nj  = nmo[Kj];
      int na0 = nocc0[Ka];
      int nj0 = nmo0[Kj];
      //auto G_(to_address(GKK[Q][K].origin()));
      auto G_(B + (Q * nkpts + K) * nwalk * nocc_max * npol * nmo_max);
      for (int a = 0, a0 = 0; a < na; a++)
      {
        for (int p = 0; p < npol; p++, a0 += nmo_max * nwalk)
        {
          //auto Gc_( to_address(Gca[na0+a][p][nj0].origin()) );
          auto Gc_(A + ((na0 + a) * npol + p) * nmo_tot * nwalk + nj0 * nwalk);
          for (int j = 0, apj = a0; j < nj; j++, apj += nwalk)
          {
            for (int w = 0; w < nwalk; w++, ++Gc_)
            {
              G_[apj + w] = static_cast<T>(*Gc_);
            }
          }
        }
      }
    }
  }
}

template<typename T, typename Q>
void vKKwij_to_vwKiKj(int nwalk, int nkpts, int nmo_max, int nmo_tot, int* kk, int* nopk, int* nopk0, Q const* A, T* B)
{
  for (int w = 0; w < nwalk; w++)
  {
    for (int Ki = 0; Ki < nkpts; Ki++)
    {
      for (int Kj = 0; Kj < nkpts; Kj++)
      {
        int ni  = nopk[Ki];
        int nj  = nopk[Kj];
        int ni0 = nopk0[Ki];
        int nj0 = nopk0[Kj];
        // setup copy/transpose tags
        // 1: copy from [Ki][Kj] without rho^+ term
        // 2: transpose from [Ki][Kj] without rho^+ term
        // -P: copy from [Ki][Kj] and transpose from [nkpts+P-1][Kj]
        int key = kk[Ki * nkpts + Kj];
        if (key == 3)
          continue;
        if (key == 2)
        { // transpose
          auto vb_(B + w * nmo_tot * nmo_tot + ni0 * nmo_tot + nj0);
          auto v_(A + ((Ki * nkpts + Kj) * nwalk + w) * nmo_max * nmo_max);
          for (int i = 0; i < ni; i++)
            for (int j = 0; j < nj; j++)
              vb_[i * nmo_tot + j] += static_cast<T>(v_[j * nmo_max + i]);
        }
        else if ((key == 1) || (key < 0))
        { // copy
          for (int i = 0; i < ni; i++)
          {
            auto vb_(B + w * nmo_tot * nmo_tot + (ni0 + i) * nmo_tot + nj0);
            auto v_(A + (((Ki * nkpts + Kj) * nwalk + w) * nmo_max + i) * nmo_max);
            for (int j = 0; j < nj; j++)
              vb_[j] += static_cast<T>(v_[j]);
          }
        }
        else
        {
          APP_ABORT(" Error: Programming error. \n");
        }
        if (key < 0)
        { // transpose
          key = (-key) - 1;
          auto vb_(B + w * nmo_tot * nmo_tot + nj0 * nmo_tot + ni0);
          auto v_(A + (((nkpts + key) * nkpts + Kj) * nwalk + w) * nmo_max * nmo_max);
          for (int i = 0; i < ni; i++)
            for (int j = 0; j < nj; j++)
              vb_[j * nmo_tot + i] += static_cast<T>(v_[i * nmo_max + j]);
        }
      }
    }
  }
}

/*
 * Performs the generic operation: (limited to matrices for now)
 * A[i,j] = A[i,j] op x[...], 
 *   where op is {+,-,*,/} and x[...] depends on dim (0:i, 1:j, ...}
 */
template<typename T, typename T2>
void term_by_term_matrix_vector(TENSOR_OPERATIONS op, int dim, int nrow, int ncol, T* A, int lda, T2 const* x, int incx)
{
  assert(dim == 0 || dim == 1);
  if (op == TOp_PLUS)
  {
    if (dim == 0)
    {
      // A[i,j] += x[i]
      for (int i = 0; i < nrow; i++, A += lda, x += incx)
        for (int j = 0; j < ncol; j++)
          A[j] += *x;
    }
    else if (dim == 1)
    {
      // A[i,j] += x[j]
      for (int i = 0; i < nrow; i++, A += lda)
        for (int j = 0; j < ncol; j++)
          A[j] += x[j * incx];
    }
  }
  else if (op == TOp_MINUS)
  {
    if (dim == 0)
    {
      // A[i,j] += x[i]
      for (int i = 0; i < nrow; i++, A += lda, x += incx)
        for (int j = 0; j < ncol; j++)
          A[j] -= *x;
    }
    else if (dim == 1)
    {
      // A[i,j] += x[j]
      for (int i = 0; i < nrow; i++, A += lda)
        for (int j = 0; j < ncol; j++)
          A[j] -= x[j * incx];
    }
  }
  else if (op == TOp_MUL)
  {
    if (dim == 0)
    {
      // A[i,j] += x[i]
      for (int i = 0; i < nrow; i++, A += lda, x += incx)
        for (int j = 0; j < ncol; j++)
          A[j] *= *x;
    }
    else if (dim == 1)
    {
      // A[i,j] += x[j]
      for (int i = 0; i < nrow; i++, A += lda)
        for (int j = 0; j < ncol; j++)
          A[j] *= x[j * incx];
    }
  }
  else if (op == TOp_DIV)
  {
    if (dim == 0)
    {
      // A[i,j] += x[i]
      for (int i = 0; i < nrow; i++, A += lda, x += incx)
        for (int j = 0; j < ncol; j++)
          A[j] /= *x;
    }
    else if (dim == 1)
    {
      // A[i,j] += x[j]
      for (int i = 0; i < nrow; i++, A += lda)
        for (int j = 0; j < ncol; j++)
          A[j] /= x[j * incx];
    }
  }
  else
  {
    APP_ABORT(" Error: Unknown operation in term_by_term_matrix_vector. \n");
  }
}

template<typename T, typename Q>
void transpose_wabn_to_wban(int nwalk, int na, int nb, int nchol, T const* Tab, Q* Tba)
{
  for (int w = 0; w < nwalk; w++)
  {
    for (int b = 0; b < nb; b++)
    {
      for (int a = 0; a < na; a++)
      {
        T const* Tn(Tab + ((w * na + a) * nb + b) * nchol);
        for (int n = 0; n < nchol; n++, ++Tba, ++Tn)
          *Tba = static_cast<Q>(*Tn);
      }
    }
  }
}

} //namespace ma

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
namespace device
{
template<typename T, typename Q>
void KaKjw_to_KKwaj(int nwalk,
                    int nkpts,
                    int npol,
                    int nmo_max,
                    int nmo_tot,
                    int nocc_max,
                    device_pointer<int> nmo,
                    device_pointer<int> nmo0,
                    device_pointer<int> nocc,
                    device_pointer<int> nocc0,
                    device_pointer<Q> A,
                    device_pointer<T> B)
{
  kernels::KaKjw_to_KKwaj(nwalk, nkpts, npol, nmo_max, nmo_tot, nocc_max, to_address(nmo), to_address(nmo0),
                          to_address(nocc), to_address(nocc0), to_address(A), to_address(B));
}

template<typename T, typename Q>
void KaKjw_to_QKajw(int nwalk,
                    int nkpts,
                    int npol,
                    int nmo_max,
                    int nmo_tot,
                    int nocc_max,
                    device_pointer<int> nmo,
                    device_pointer<int> nmo0,
                    device_pointer<int> nocc,
                    device_pointer<int> nocc0,
                    device_pointer<int> QKtok2,
                    device_pointer<Q> A,
                    device_pointer<T> B)
{
  kernels::KaKjw_to_QKajw(nwalk, nkpts, npol, nmo_max, nmo_tot, nocc_max, to_address(nmo), to_address(nmo0),
                          to_address(nocc), to_address(nocc0), to_address(QKtok2), to_address(A), to_address(B));
}

template<typename T, typename Q>
void vKKwij_to_vwKiKj(int nwalk,
                      int nkpts,
                      int nmo_max,
                      int nmo_tot,
                      device_pointer<int> kk,
                      device_pointer<int> nmo,
                      device_pointer<int> nmo0,
                      device_pointer<Q> A,
                      device_pointer<T> B)
{
  kernels::vKKwij_to_vwKiKj(nwalk, nkpts, nmo_max, nmo_tot, to_address(kk), to_address(nmo), to_address(nmo0),
                            to_address(A), to_address(B));
}

/*
 * Performs the generic operation: (limited to matrices for now)
 * A[i,j] = A[i,j] op x[...], 
 *   where op is {+,-,*,/} and x[...] depends on dim (0:i, 1:j, ...}
 */
template<typename T, typename T2>
void term_by_term_matrix_vector(ma::TENSOR_OPERATIONS op,
                                int dim,
                                int nrow,
                                int ncol,
                                device_pointer<T> A,
                                int lda,
                                device_pointer<T2> x,
                                int incx)
{
  assert(dim == 0 || dim == 1);
  if (op == ma::TOp_PLUS)
    kernels::term_by_term_mat_vec_plus(dim, nrow, ncol, to_address(A), lda, to_address(x), incx);
  else if (op == ma::TOp_MINUS)
    kernels::term_by_term_mat_vec_minus(dim, nrow, ncol, to_address(A), lda, to_address(x), incx);
  else if (op == ma::TOp_MUL)
    kernels::term_by_term_mat_vec_mult(dim, nrow, ncol, to_address(A), lda, to_address(x), incx);
  else if (op == ma::TOp_DIV)
    kernels::term_by_term_mat_vec_div(dim, nrow, ncol, to_address(A), lda, to_address(x), incx);
  else
    APP_ABORT(" Error: Unknown operation in term_by_term_matrix_vector. \n");
}

template<typename T, typename Q>
void transpose_wabn_to_wban(int nwalk, int na, int nb, int nchol, device_pointer<T> Tab, device_pointer<Q> Tba)
{
  kernels::transpose_wabn_to_wban(nwalk, na, nb, nchol, to_address(Tab), to_address(Tba));
}

} // namespace device
#endif

#endif
