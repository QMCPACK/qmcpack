//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_KP3INDEXIO_HPP
#define QMCPLUSPLUS_AFQMC_KP3INDEXIO_HPP

#include <fstream>

#include "hdf/hdf_archive.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/kp_utilities.hpp"
#include "AFQMC/SlaterDeterminantOperations/rotate.hpp"

#include "AFQMC/HamiltonianOperations/KP3IndexFactorization.hpp"

namespace qmcplusplus
{
namespace afqmc
{
inline HamiltonianOperations loadKP3IndexFactorization(hdf_archive& dump,
                                                       WALKER_TYPES type,
                                                       int NMO,
                                                       int NAEA,
                                                       int NAEB,
                                                       std::vector<PsiT_Matrix>& PsiT,
                                                       TaskGroup_& TGprop,
                                                       TaskGroup_& TGwfn,
                                                       RealType cutvn,
                                                       RealType cutv2)
{
  using shmIMatrix    = boost::multi::array<int, 2, shared_allocator<int>>;
  using shmCVector    = boost::multi::array<ComplexType, 1, shared_allocator<ComplexType>>;
  using shmCMatrix    = boost::multi::array<ComplexType, 2, shared_allocator<ComplexType>>;
  using shmCTensor    = boost::multi::array<ComplexType, 3, shared_allocator<ComplexType>>;
  using shmSpMatrix   = boost::multi::array<SPComplexType, 2, shared_allocator<SPComplexType>>;
  using shmSpTensor   = boost::multi::array<SPComplexType, 3, shared_allocator<SPComplexType>>;
  using SpMatrix      = boost::multi::array<SPComplexType, 2>;
  using SpMatrix_ref  = boost::multi::array_ref<SPComplexType, 2>;
  using Sp3Tensor_ref = boost::multi::array_ref<SPComplexType, 3>;

  if (TGprop.getNGroupsPerTG() > 1)
    APP_ABORT(" Error: Distributed KP3INdex not yet implemented.\n");
  if (TGwfn.getNGroupsPerTG() > 1)
    APP_ABORT(" Error: Distributed KP3INdex not yet implemented.\n");

  std::vector<int> dims(7);
  int ndet   = (type == COLLINEAR ? PsiT.size() / 2 : PsiT.size());
  int nspins = ((type != COLLINEAR) ? 1 : 2);
  if (type == COLLINEAR)
    assert(PsiT.size() % 2 == 0);
  ValueType E0;

  // single reader for now
  if (TGwfn.Global().root())
  {
    dump.push("HamiltonianOperations", false);
    dump.push("KP3IndexFactorization", false);

    if (!dump.read(dims, "dims"))
    {
      app_error() << " Error in loadKP3IndexFactorization: Problems reading dataset. \n";
      APP_ABORT("");
    }
    if (dims[0] != NMO)
    {
      app_error() << " Error in loadKP3IndexFactorization: Inconsistent data in file: NMO. \n";
      APP_ABORT("");
    }
    if (dims[1] != NAEA)
    {
      app_error() << " Error in loadKP3IndexFactorization: Inconsistent data in file: NAEA. \n";
      APP_ABORT("");
    }
    if (dims[2] != NAEB)
    {
      app_error() << " Error in loadKP3IndexFactorization: Inconsistent data in file: NAEB. \n";
      APP_ABORT("");
    }
    if (type == CLOSED && dims[4] != 1)
    {
      app_error() << " Error in loadKP3IndexFactorization: Inconsistent data in file: walker_type. \n";
      APP_ABORT("");
    }
    if (type == COLLINEAR && dims[4] != 2)
    {
      app_error() << " Error in loadKP3IndexFactorization: Inconsistent data in file: walker_type. \n";
      APP_ABORT("");
    }
    if (type == NONCOLLINEAR && dims[4] != 3)
    {
      app_error() << " Error in loadKP3IndexFactorization: Inconsistent data in file: walker_type. \n";
      APP_ABORT("");
    }
    std::vector<ValueType> et;
    if (!dump.read(et, "E0"))
    {
      app_error() << " Error in loadKP3IndexFactorization: Problems reading dataset. \n";
      APP_ABORT("");
    }
    E0 = et[0];
  }
  TGwfn.Global().broadcast_n(dims.data(), 7);
  TGwfn.Global().broadcast_value(E0);
  int nkpts    = dims[3];
  int nsampleQ = dims[5];

  std::vector<int> nmo_per_kp(nkpts);
  std::vector<int> nchol_per_kp(nkpts);
  std::vector<int> kminus(nkpts);
  std::vector<RealType> gQ;
  shmIMatrix QKtok2({nkpts, nkpts}, shared_allocator<int>{TGwfn.Node()});

  if (TGwfn.Global().root())
  {
    if (!dump.read(nmo_per_kp, "NMOPerKP"))
    {
      app_error() << " Error in loadKP3IndexFactorization: Problems reading NMOPerKP dataset. \n";
      APP_ABORT("");
    }
    if (!dump.read(nchol_per_kp, "NCholPerKP"))
    {
      app_error() << " Error in loadKP3IndexFactorization: Problems reading NCholPerKP dataset. \n";
      APP_ABORT("");
    }
    if (!dump.read(kminus, "MinusK"))
    {
      app_error() << " Error in loadKP3IndexFactorization: Problems reading MinusK dataset. \n";
      APP_ABORT("");
    }
    if (!dump.read(QKtok2, "QKTok2"))
    {
      app_error() << " Error in loadKP3IndexFactorization: Problems reading QKTok2 dataset. \n";
      APP_ABORT("");
    }
    if (!dump.read(gQ, "gQ"))
    {
      app_error() << " Error in loadKP3IndexFactorization: Problems reading gQ dataset. \n";
      APP_ABORT("");
    }
  }
  TGwfn.Global().broadcast_n(nmo_per_kp.data(), nmo_per_kp.size(), 0);
  TGwfn.Global().broadcast_n(nchol_per_kp.data(), nchol_per_kp.size(), 0);
  TGwfn.Global().broadcast_n(kminus.data(), kminus.size(), 0);
  TGwfn.Global().broadcast_n(gQ.data(), gQ.size(), 0);
  if (TGwfn.Node().root())
    TGwfn.Cores().broadcast_n(std::addressof(*QKtok2.origin()), QKtok2.num_elements(), 0);

  int Q0 = -1; // stores the index of the Q=(0,0,0) Q-point
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

  int nmo_max   = *std::max_element(nmo_per_kp.begin(), nmo_per_kp.end());
  int nchol_max = *std::max_element(nchol_per_kp.begin(), nchol_per_kp.end());
  shmCTensor H1({nkpts, nmo_max, nmo_max}, shared_allocator<ComplexType>{TGwfn.Node()});
  shmCTensor vn0({nkpts, nmo_max, nmo_max}, shared_allocator<ComplexType>{TGwfn.Node()});
  std::vector<shmSpMatrix> LQKikn;
  LQKikn.reserve(nkpts);
  for (int Q = 0; Q < nkpts; Q++)
    if (Q == Q0)
      LQKikn.emplace_back(
          shmSpMatrix({nkpts, nmo_max * nmo_max * nchol_per_kp[Q]}, shared_allocator<SPComplexType>{TGwfn.Node()}));
    else if (kminus[Q] == Q) // only storing half of K points and using symmetry
      LQKikn.emplace_back(
          shmSpMatrix({nkpts / 2, nmo_max * nmo_max * nchol_per_kp[Q]}, shared_allocator<SPComplexType>{TGwfn.Node()}));
    else if (Q < kminus[Q])
      LQKikn.emplace_back(
          shmSpMatrix({nkpts, nmo_max * nmo_max * nchol_per_kp[Q]}, shared_allocator<SPComplexType>{TGwfn.Node()}));
    else // Q > kminus[Q]
      LQKikn.emplace_back(shmSpMatrix({1, 1}, shared_allocator<SPComplexType>{TGwfn.Node()}));

  if (TGwfn.Global().root())
  {
    {
      boost::multi::array_ref<ComplexType, 2> H1_(std::addressof(*H1.origin()),
                                                  {H1.shape()[0], H1.shape()[1] * H1.shape()[2]});
      if (!dump.read(H1_, "KPH1"))
      {
        app_error() << " Error in loadKP3IndexFactorization: Problems reading dataset. \n";
        APP_ABORT("");
      }
    }
    {
      boost::multi::array_ref<ComplexType, 2> vn0_(std::addressof(*vn0.origin()),
                                                   {vn0.shape()[0], vn0.shape()[1] * vn0.shape()[2]});
      if (!dump.read(vn0_, "KPv0"))
      {
        app_error() << " Error in loadKP3IndexFactorization: Problems reading dataset. \n";
        APP_ABORT("");
      }
    }
    for (int Q = 0; Q < nkpts; Q++)
    {
      if (Q > kminus[Q])
        continue;
      // to avoid reallocation, will abort if sizes don't match
      boost::multi::array_ref<SPComplexType, 2> LQ(std::addressof(*LQKikn[Q].origin()),
                                                   {LQKikn[Q].shape()[0], LQKikn[Q].shape()[1]});
      if (!dump.read(LQ, std::string("L") + std::to_string(Q)))
      {
        app_error() << " Error in loadKP3IndexFactorization: Problems reading dataset. \n";
        APP_ABORT("");
      }
    }
  }
  if (TGwfn.Node().root())
  {
    TGwfn.Cores().broadcast_n(std::addressof(*H1.origin()), H1.num_elements(), 0);
    TGwfn.Cores().broadcast_n(std::addressof(*vn0.origin()), vn0.num_elements(), 0);
    for (int Q = 0; Q < nkpts; Q++)
    {
      if (Q > kminus[Q])
        continue;
      TGwfn.Cores().broadcast_n(std::addressof(*LQKikn[Q].origin()), LQKikn[Q].num_elements(), 0);
    }
  }


  shmIMatrix nocc_per_kp({ndet, nspins * nkpts}, shared_allocator<int>{TGwfn.Node()});
  if (TGwfn.Node().root())
  {
    if (type == COLLINEAR)
    {
      for (int i = 0; i < ndet; i++)
      {
        if (not get_nocc_per_kp(nmo_per_kp, PsiT[2 * i], nocc_per_kp[i]({0, nkpts})))
        {
          app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                      << " Only wavefunctions in block-diagonal form are accepted. \n";
          APP_ABORT("");
        }
        if (not get_nocc_per_kp(nmo_per_kp, PsiT[2 * i + 1], nocc_per_kp[i]({nkpts, 2 * nkpts})))
        {
          app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                      << " Only wavefunctions in block-diagonal form are accepted. \n";
          APP_ABORT("");
        }
      }
    }
    else
    {
      for (int i = 0; i < ndet; i++)
        if (not get_nocc_per_kp(nmo_per_kp, PsiT[i], nocc_per_kp[i]))
        {
          app_error() << " Error in KPFactorizedHamiltonian::getHamiltonianOperations():"
                      << " Only wavefunctions in block-diagonal form are accepted. \n";
          APP_ABORT("");
        }
    }
  }
  TGwfn.Node().barrier();
  int nocc_max = *std::max_element(std::addressof(*nocc_per_kp.origin()),
                                   std::addressof(*nocc_per_kp.origin()) + nocc_per_kp.num_elements());


  std::vector<shmSpMatrix> LQKank;
  LQKank.reserve(ndet * nspins * (nkpts + 1)); // storing 2 components for Q=0, since it is not assumed symmetric
  shmCMatrix haj({ndet * nkpts, (type == COLLINEAR ? 2 : 1) * nocc_max * nmo_max},
                 shared_allocator<ComplexType>{TGwfn.Node()});
  if (TGwfn.Node().root())
    std::fill_n(haj.origin(), haj.num_elements(), ComplexType(0.0));
  int ank_max = nocc_max * nchol_max * nmo_max;
  for (int nd = 0; nd < ndet; nd++)
  {
    for (int Q = 0; Q < (nkpts + 1); Q++)
    {
      LQKank.emplace_back(shmSpMatrix({nkpts, ank_max}, shared_allocator<ComplexType>{TGwfn.Node()}));
    }
    if (type == COLLINEAR)
    {
      for (int Q = 0; Q < (nkpts + 1); Q++)
      {
        LQKank.emplace_back(shmSpMatrix({nkpts, ank_max}, shared_allocator<ComplexType>{TGwfn.Node()}));
      }
    }
  }
  for (int nd = 0, nt = 0, nq0 = 0; nd < ndet; nd++, nq0 += (nkpts + 1) * nspins)
  {
    for (int Q = 0; Q < (nkpts + 1); Q++)
    {
      for (int K = 0; K < nkpts; K++, nt++)
      {
        if (nt % TGwfn.Node().size() == TGwfn.Node().rank())
        {
          std::fill_n(std::addressof(*LQKank[nq0 + Q][K].origin()), LQKank[nq0 + Q][K].num_elements(),
                      SPComplexType(0.0));
          if (type == COLLINEAR)
          {
            std::fill_n(std::addressof(*LQKank[nq0 + nkpts + 1 + Q][K].origin()),
                        LQKank[nq0 + nkpts + 1 + Q][K].num_elements(), SPComplexType(0.0));
          }
        }
      }
    }
  }
  TGwfn.Node().barrier();
  boost::multi::array<SPComplexType, 2> buff({nmo_max, nchol_max});
  for (int nd = 0, nt = 0, nq0 = 0; nd < ndet; nd++, nq0 += (nkpts + 1) * nspins)
  {
    for (int Q = 0; Q < nkpts; Q++)
    {
      for (int K = 0; K < nkpts; K++, nt++)
      {
        if (nt % TGwfn.Global().size() == TGwfn.Global().rank())
        {
          // haj and add half-transformed right-handed rotation for Q=0
          int Qm    = kminus[Q];
          int QK    = QKtok2[Q][K];
          int na    = nocc_per_kp[nd][K];
          int nb    = nocc_per_kp[nd][nkpts + QK];
          int ni    = nmo_per_kp[K];
          int nk    = nmo_per_kp[QK];
          int nchol = nchol_per_kp[Q];
          if (Q == 0)
          {
            if (type == COLLINEAR)
            {
              { // Alpha
                auto Psi = get_PsiK<boost::multi::array<ComplexType, 2>>(nmo_per_kp, PsiT[2 * nd], K);
                assert(Psi.shape()[0] == na);
                boost::multi::array_ref<ComplexType, 2> haj_r(std::addressof(*haj[nd * nkpts + K].origin()), {na, ni});
                ma::product(Psi, H1[K]({0, ni}, {0, ni}), haj_r);
              }
              { // Beta
                auto Psi = get_PsiK<boost::multi::array<ComplexType, 2>>(nmo_per_kp, PsiT[2 * nd + 1], K);
                assert(Psi.shape()[0] == nb);
                boost::multi::array_ref<ComplexType, 2> haj_r(std::addressof(*haj[nd * nkpts + K].origin()) + na * ni,
                                                              {nb, ni});
                ma::product(Psi, H1[K]({0, ni}, {0, ni}), haj_r);
              }
            }
            else
            {
              auto Psi = get_PsiK<boost::multi::array<ComplexType, 2>>(nmo_per_kp, PsiT[nd], K);
              assert(Psi.shape()[0] == na);
              boost::multi::array_ref<ComplexType, 2> haj_r(std::addressof(*haj[nd * nkpts + K].origin()), {na, ni});
              ma::product(ComplexType(2.0), Psi, H1[K]({0, ni}, {0, ni}), ComplexType(0.0), haj_r);
            }
          }
          if (type == COLLINEAR)
          {
            { // Alpha
              // change get_PsiK to cast to the value_type of the result
              auto Psi = get_PsiK<boost::multi::array<SPComplexType, 2>>(nmo_per_kp, PsiT[2 * nd], K);
              assert(Psi.shape()[0] == nocc_per_kp[nd][K]);
              if (Q < Qm || Q == Q0 || ((Q == Qm) && (K < QK)))
              {
                int kpos = K;
                if (Q == Qm && Q != Q0)
                { //find position in symmetric list
                  kpos = 0;
                  for (int K_ = 0; K_ < K; K_++)
                    if (K_ < QKtok2[Q][K_])
                      kpos++;
                }
                Sp3Tensor_ref Likn(std::addressof(*LQKikn[Q][kpos].origin()), {ni, nk, nchol});
                Sp3Tensor_ref Lank(std::addressof(*LQKank[nq0 + Q][K].origin()), {na, nchol, nk});
                ma_rotate::getLank(Psi, Likn, Lank, buff);
                if (Q == Q0)
                {
                  assert(K == QK);
                  Sp3Tensor_ref Lank(std::addressof(*LQKank[nq0 + nkpts][K].origin()), {na, nchol, nk});
                  ma_rotate::getLank_from_Lkin(Psi, Likn, Lank, buff);
                }
              }
              else
              {
                int kpos = QK;
                if (Q == Qm)
                { //find position in symmetric list
                  kpos = 0;
                  for (int K_ = 0; K_ < QK; K_++)
                    if (K_ < QKtok2[Q][K_])
                      kpos++;
                }
                Sp3Tensor_ref Lkin(std::addressof(*LQKikn[Qm][QK].origin()), {nk, ni, nchol});
                Sp3Tensor_ref Lank(std::addressof(*LQKank[nq0 + Q][K].origin()), {na, nchol, nk});
                ma_rotate::getLank_from_Lkin(Psi, Lkin, Lank, buff);
              }
            }
            { // Beta
              // change get_PsiK to cast to the value_type of the result
              auto Psi = get_PsiK<boost::multi::array<SPComplexType, 2>>(nmo_per_kp, PsiT[2 * nd + 1], K);
              assert(Psi.shape()[0] == nb);
              if (Q < Qm || Q == Q0 || ((Q == Qm) && (K < QK)))
              {
                int kpos = K;
                if (Q == Qm && Q != Q0)
                { //find position in symmetric list
                  kpos = 0;
                  for (int K_ = 0; K_ < K; K_++)
                    if (K_ < QKtok2[Q][K_])
                      kpos++;
                }
                Sp3Tensor_ref Likn(std::addressof(*LQKikn[Q][kpos].origin()), {ni, nk, nchol});
                Sp3Tensor_ref Lank(std::addressof(*LQKank[nq0 + nkpts + 1 + Q][K].origin()), {nb, nchol, nk});
                ma_rotate::getLank(Psi, Likn, Lank, buff);
                if (Q == Q0)
                {
                  assert(K == QK);
                  Sp3Tensor_ref Lank(std::addressof(*LQKank[nq0 + nkpts + 1 + nkpts][K].origin()), {nb, nchol, nk});
                  ma_rotate::getLank_from_Lkin(Psi, Likn, Lank, buff);
                }
              }
              else
              {
                int kpos = QK;
                if (Q == Qm)
                { //find position in symmetric list
                  kpos = 0;
                  for (int K_ = 0; K_ < QK; K_++)
                    if (K_ < QKtok2[Q][K_])
                      kpos++;
                }
                Sp3Tensor_ref Lkin(std::addressof(*LQKikn[Qm][QK].origin()), {nk, ni, nchol});
                Sp3Tensor_ref Lank(std::addressof(*LQKank[nq0 + nkpts + 1 + Q][K].origin()), {nb, nchol, nk});
                ma_rotate::getLank_from_Lkin(Psi, Lkin, Lank, buff);
              }
            }
          }
          else
          {
            // change get_PsiK to cast to the value_type of the result
            auto Psi = get_PsiK<SpMatrix>(nmo_per_kp, PsiT[nd], K);
            assert(Psi.shape()[0] == na);
            if (Q < Qm || Q == Q0 || ((Q == Qm) && (K < QK)))
            {
              int kpos = K;
              if (Q == Qm && Q != Q0)
              { //find position in symmetric list
                kpos = 0;
                for (int K_ = 0; K_ < K; K_++)
                  if (K_ < QKtok2[Q][K_])
                    kpos++;
              }
              Sp3Tensor_ref Likn(std::addressof(*LQKikn[Q][kpos].origin()), {ni, nk, nchol});
              Sp3Tensor_ref Lank(std::addressof(*LQKank[nq0 + Q][K].origin()), {na, nchol, nk});
              ma_rotate::getLank(Psi, Likn, Lank, buff);
              if (Q == Q0)
              {
                assert(K == QK);
                Sp3Tensor_ref Lank(std::addressof(*LQKank[nq0 + nkpts][K].origin()), {na, nchol, nk});
                ma_rotate::getLank_from_Lkin(Psi, Likn, Lank, buff);
              }
            }
            else
            {
              int kpos = QK;
              if (Q == Qm)
              { //find position in symmetric list
                kpos = 0;
                for (int K_ = 0; K_ < QK; K_++)
                  if (K_ < QKtok2[Q][K_])
                    kpos++;
              }
              Sp3Tensor_ref Lkin(std::addressof(*LQKikn[Qm][kpos].origin()), {nk, ni, nchol});
              Sp3Tensor_ref Lank(std::addressof(*LQKank[nq0 + Q][K].origin()), {na, nchol, nk});
              ma_rotate::getLank_from_Lkin(Psi, Lkin, Lank, buff);
            }
          }
        }
      }
    }
  }
  TGwfn.Global().barrier();
  if (TGwfn.Node().root())
  {
    TGwfn.Cores().all_reduce_in_place_n(std::addressof(*haj.origin()), haj.num_elements(), std::plus<>());
    for (int Q = 0; Q < LQKank.size(); Q++)
      TGwfn.Cores().all_reduce_in_place_n(std::addressof(*LQKank[Q].origin()), LQKank[Q].num_elements(), std::plus<>());
  }
  TGwfn.Node().barrier();

  int global_ncvecs = std::accumulate(nchol_per_kp.begin(), nchol_per_kp.end(), 0);

  if (TGwfn.Global().root())
  {
    dump.pop();
    dump.pop();
  }

  return HamiltonianOperations(
      KP3IndexFactorization(TGwfn.TG_local(), type, std::move(nmo_per_kp), std::move(nchol_per_kp), std::move(kminus),
                            std::move(nocc_per_kp), std::move(QKtok2), std::move(H1), std::move(haj), std::move(LQKikn),
                            std::move(LQKank), std::move(vn0), std::move(gQ), nsampleQ, E0, global_ncvecs));
}

// single writer right now
template<class shmIMatrix, class shmC3Tensor, class shmSpMatrix>
inline void writeKP3IndexFactorization(hdf_archive& dump,
                                       WALKER_TYPES type,
                                       int NMO,
                                       int NAEA,
                                       int NAEB,
                                       TaskGroup_& TGprop,
                                       TaskGroup_& TGwfn,
                                       std::vector<int>& nopk,
                                       std::vector<int>& ncholpQ,
                                       std::vector<int>& kminus,
                                       shmIMatrix& QKToK2,
                                       shmC3Tensor& H1,
                                       std::vector<shmSpMatrix>& Lik,
                                       shmC3Tensor& vn0,
                                       int nsampleQ,
                                       std::vector<RealType>& gQ,
                                       ValueType E0,
                                       int gncv)
{
  if (TGprop.getNGroupsPerTG() > 1)
    APP_ABORT(" Error: Distributed KP3INdex not yet implemented.\n");
  if (TGwfn.getNGroupsPerTG() > 1)
    APP_ABORT(" Error: Distributed KP3INdex not yet implemented.\n");

  if (TGwfn.Global().root())
  {
    dump.push("HamiltonianOperations");
    dump.push("KP3IndexFactorization");

    std::vector<int> dims{NMO, NAEA, NAEB, int(nopk.size()), type, nsampleQ, gncv};
    dump.write(dims, "dims");
    std::vector<ValueType> et{E0};
    dump.write(nopk, "NMOPerKP");
    dump.write(ncholpQ, "NCholPerKP");
    dump.write(kminus, "MinusK");
    dump.write(QKToK2, "QKTok2");
    dump.write(gQ, "gQ");
    {
      boost::multi::array_ref<ComplexType, 2> H1_(std::addressof(*H1.origin()),
                                                  {H1.shape()[0], H1.shape()[1] * H1.shape()[2]});
      dump.write(H1_, "KPH1");
    }
    {
      boost::multi::array_ref<ComplexType, 2> vn0_(std::addressof(*vn0.origin()),
                                                   {vn0.shape()[0], vn0.shape()[1] * vn0.shape()[2]});
      dump.write(vn0_, "KPv0");
    }
    int nkpts = Lik.size();
    // use hyperslabs on the roots of TGnumber=0 in distributed case
    for (int Q = 0; Q < nkpts; Q++)
      dump.write(Lik[Q], std::string("L") + std::to_string(Q));
    dump.write(et, "E0");
  }

  if (TGwfn.Global().root())
  {
    dump.pop();
    dump.pop();
  }
  TGwfn.Global().barrier();
}

} // namespace afqmc
} // namespace qmcplusplus

#endif
