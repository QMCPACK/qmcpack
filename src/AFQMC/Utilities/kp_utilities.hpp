////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by:
// Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_KP_UTILITIES_HPP
#define AFQMC_KP_UTILITIES_HPP

#include <algorithm>
#include "AFQMC/Numerics/ma_operations.hpp"
#include "multi/array.hpp"
#include "multi/array_ref.hpp"

/*
 * Return true if PsiT is in block diagonal form, otherwise return false;
 * If it is block diagonal form, it also returns the number of states in each k-point block. 
 */
template<class Vector, class CSR, class Array>
bool get_nocc_per_kp(Vector const& nmo_per_kp, CSR const& PsiT, Array&& nocc_per_kp, bool noncolin = false)
{
  int nkpts = nmo_per_kp.size();
  int N     = PsiT.size(0);
  int M     = PsiT.size(1);
  int npol  = noncolin ? 2 : 1;
  assert(M % npol == 0);
  assert(nocc_per_kp.size() == nkpts);

  std::fill_n(to_address(nocc_per_kp.origin()), nkpts, 0);
  std::vector<int> bounds(npol * nkpts + 1);
  bounds[0] = 0;
  for (int k = 0; k < npol * nkpts; k++)
    bounds[k + 1] = bounds[k] + nmo_per_kp[k % nkpts];
  int Q = 0;
  for (int i = 0; i < N; i++)
  {
    auto nt = PsiT.num_non_zero_elements(i);
    if (nt == 0)
    {
      std::fill_n(to_address(nocc_per_kp.origin()), nkpts, 0);
      return false;
    }
    auto col = PsiT.non_zero_indices2_data(i);
    // check the kp index of the first non-zero column. Mut be either >= Q
    auto it = std::lower_bound(bounds.begin(), bounds.end(), *col + 1) - 1;
    assert(it != bounds.end());
    int Q_   = std::distance(bounds.begin(), it) % nkpts;
    int pol_ = std::distance(bounds.begin(), it) / nkpts;
    assert(Q_ >= 0 && Q_ < nkpts);
    assert(pol_ == 0 || pol_ == 1);
    if (Q_ < Q)
    {
      std::fill_n(to_address(nocc_per_kp.origin()), nkpts, 0);
      return false;
    }
    Q = Q_;
    for (int ni = 0; ni < nt; ++ni, ++col)
    {
      if ((*col < bounds[Q] || *col >= bounds[Q + 1]) &&
          (*col < bounds[(npol - 1) * nkpts + Q] || *col >= bounds[(npol - 1) * nkpts + Q + 1]))
      {
        std::fill_n(to_address(nocc_per_kp.origin()), nkpts, 0);
        return false;
      }
    }
    ++nocc_per_kp[Q];
  }
  return true;
}

template<class Array, class Vector, class CSR>
Array get_PsiK(Vector const& nmo_per_kp, CSR const& PsiT, int K, bool noncolin = false)
{
  int nkpts = nmo_per_kp.size();
  int N     = PsiT.size(0);
  int M     = PsiT.size(1);
  int npol  = noncolin ? 2 : 1;
  assert(M % npol == 0);

  int nel = 0;
  std::vector<int> bounds(npol * nkpts + 1);
  bounds[0] = 0;
  for (int k = 0; k < npol * nkpts; k++)
    bounds[k + 1] = bounds[k] + nmo_per_kp[k % nkpts];
  int Q = 0;
  for (int i = 0; i < N; i++)
  {
    auto nt = PsiT.num_non_zero_elements(i);
    if (nt == 0)
      APP_ABORT("Error: PsiT not in block-diagonal form in get_PsiK.\n");
    auto col = PsiT.non_zero_indices2_data(i);
    // check the kp index of the first non-zero column. Mut be either >= Q
    auto it = std::lower_bound(bounds.begin(), bounds.end(), *col + 1) - 1;
    assert(it != bounds.end());
    int Q_   = std::distance(bounds.begin(), it) % nkpts;
    int pol_ = std::distance(bounds.begin(), it) / nkpts;
    assert(Q_ >= 0 && Q_ < nkpts);
    assert(pol_ == 0 || pol_ == 1);
    if (Q_ < Q)
      APP_ABORT("Error: PsiT not in block-diagonal form in get_PsiK.\n");
    Q = Q_;
    for (int ni = 0; ni < nt; ++ni, ++col)
      if ((*col < bounds[Q] || *col >= bounds[Q + 1]) &&
          (*col < bounds[(npol - 1) * nkpts + Q] || *col >= bounds[(npol - 1) * nkpts + Q + 1]))
        APP_ABORT("Error: PsiT not in block-diagonal form in get_PsiK.\n");
    if (Q == K)
      nel++;
  }
  using element = typename std::decay<Array>::type::element;
  Array A({nel, npol * nmo_per_kp[K]});
  using std::fill_n;
  fill_n(A.origin(), A.num_elements(), element(0));
  nel = 0;
  for (int i = 0; i < N; i++)
  {
    auto nt  = PsiT.num_non_zero_elements(i);
    auto col = PsiT.non_zero_indices2_data(i);
    auto val = PsiT.non_zero_values_data(i);
    // check the kp index of the first non-zero column. Mut be either >= Q
    auto it = std::lower_bound(bounds.begin(), bounds.end(), *col + 1) - 1;
    int Q   = std::distance(bounds.begin(), it) % nkpts;
    if (Q == K)
    {
      for (int ni = 0; ni < nt; ++ni, ++col, ++val)
      {
        if (*col < bounds[K + 1]) // alpha
          A[nel][*col - bounds[K]] = static_cast<element>(*val);
        else // beta
          A[nel][*col - bounds[nkpts + K] + nmo_per_kp[K]] = static_cast<element>(*val);
      }
      nel++;
    }
    if (Q > K)
      break;
  }
  return A;
}

/*
 * Checks if the 2-electron integrals remain invariant under the symmetry assumption: L_IJ_n = conj(L_JI_n)
 *  when Q==kminus[Q]
 */
template<class T1, class shmIMatrix, class IVector>
bool check_cholesky_symmetry(T1 const& LQKikn,
                             IVector const& nmo_per_kp,
                             IVector const& nchol_per_kp,
                             shmIMatrix const& QKtok2,
                             IVector const& kminus)
{
  std::cout << " Checking kpoint cholesky symmetry. \n";
  using ComplexType = typename T1::value_type::element;
  int nkpts         = nmo_per_kp.size();
  int nmo_max       = *std::max_element(nmo_per_kp.begin(), nmo_per_kp.end());
  int nchol_max     = *std::max_element(nchol_per_kp.begin(), nchol_per_kp.end());
  int Q0            = -1; // stores the index of the Q=(0,0,0) Q-point
                          // this must always exist
  for (int Q = 0; Q < nkpts; Q++)
  {
    if (kminus[Q] == Q)
    {
      bool found = true;
      for (int KI = 0; KI < nkpts; KI++)
        if (KI != QKtok2[Q][KI])
        {
          found = false;
          break;
        }
      if (found)
      {
        if (Q0 > 0)
          APP_ABORT(" Error: @ Q-points satisfy Q=0 condition.\n");
        Q0 = Q;
      }
      else
      {
        if (nkpts % 2 != 0)
          APP_ABORT(" Error: Unexpected situation: Q==(-Q)!=Q0 and odd Nk. \n");
      }
    }
  }
  if (Q0 < 0)
    APP_ABORT(" Error: Can not find Q=0 Q-point.\n");

  boost::multi::array<int, 2> KK2Q({nkpts, nkpts});
  for (int KI = 0; KI < nkpts; KI++)
    for (int KK = 0; KK < nkpts; KK++)
    {
      KK2Q[KI][KK] = -1;
      for (int Q = 0; Q < nkpts; Q++)
        if (QKtok2[Q][KI] == KK)
        {
          KK2Q[KI][KK] = Q;
          break;
        }
      assert(KK2Q[KI][KK] >= 0);
    }

  boost::multi::array<ComplexType, 2> IJKL_({nmo_max * nmo_max, nmo_max * nmo_max});
  boost::multi::array<ComplexType, 2> IJKL2_({nmo_max * nmo_max, nmo_max * nmo_max});

  double mxx(0.0);
  bool res = true;
  for (int Q = 0; Q < nkpts; ++Q)
  {
    if (Q != kminus[Q] || Q == Q0)
      continue;

    for (int KI = 0; KI < nkpts; KI++)
      for (int KL = 0; KL < nkpts; KL++)
      {
        int KK = QKtok2[Q][KI];
        int KJ = QKtok2[Q][KL];
        if (KI > KK || KL < KJ)
          continue;
        int ni = nmo_per_kp[KI];
        int nj = nmo_per_kp[KJ];
        int nk = nmo_per_kp[KK];
        int nl = nmo_per_kp[KL];
        boost::multi::array_ref<ComplexType, 2, const ComplexType*> LKI(std::addressof(*LQKikn[Q][KI].origin()),
                                                                        {ni * nk, nchol_per_kp[Q]});
        boost::multi::array_ref<ComplexType, 2, const ComplexType*> LKL(std::addressof(*LQKikn[Q][KL].origin()),
                                                                        {nl * nj, nchol_per_kp[Q]});
        boost::multi::array_ref<ComplexType, 2> IJKL(std::addressof(*IJKL_.origin()), {ni * nk, nl * nj});
        boost::multi::array_ref<ComplexType, 4> IJKL_4D(std::addressof(*IJKL_.origin()), {ni, nk, nl, nj});
        ma::product(LKI, ma::H(LKL), IJKL);

        boost::multi::array_ref<ComplexType, 2, const ComplexType*> LKJ(std::addressof(*LQKikn[Q][KJ].origin()),
                                                                        {nj * nl, nchol_per_kp[Q]});
        boost::multi::array_ref<ComplexType, 2> IJKL2(std::addressof(*IJKL2_.origin()), {ni * nk, nj * nl});
        boost::multi::array_ref<ComplexType, 4> IJKL2_4D(std::addressof(*IJKL2_.origin()), {ni, nk, nj, nl});
        ma::product(LKI, ma::T(LKJ), IJKL2);
        double mx(0.0);
        for (int i = 0; i < ni; i++)
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int l = 0; l < nl; l++)
                mx = std::max(mx, std::abs(IJKL_4D[i][k][l][j] - IJKL2_4D[i][k][j][l]));
        if (mx > 1e-12)
        {
          std::cout << " WARNING: Matrix elements are not symmetric in check_cholesky_symmetry: " << Q << " " << KI
                    << " " << KK << " " << KJ << " " << KL << " " << mx << std::endl;
          res = false;
        }
        mxx = std::max(mx, mxx);
      }
  }
  std::cout << " Largest ERI difference found: " << mxx << std::endl;
  return res;
}

#endif
