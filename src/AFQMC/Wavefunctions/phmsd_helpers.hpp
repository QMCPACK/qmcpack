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

#ifndef QMCPLUSPLUS_AFQMC_PHMSD_HELPERS_HPP
#define QMCPLUSPLUS_AFQMC_PHMSD_HELPERS_HPP

#include "AFQMC/config.h"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Numerics/ma_small_mat_ops.hpp"

namespace qmcplusplus
{
namespace afqmc
{
// for GPU implementation, dispatch Kernels from here!
// it is probably easier to have 2 types of kernels, one that loads the
// appropriate terms in the matrices and a second one that just computes determinants

// using simple round-robin scheme for parallelization within TG_local
// assumes that reference determinant is already on [0]
template<class Array1D, class MatA, class MatB, class PH_EXCT>
inline void calculate_overlaps(int rank, int ngrp, int spin, PH_EXCT const& abij, MatA&& T, MatB&& Qwork, Array1D&& ov)
{
  std::vector<int> IWORK(abij.maximum_excitation_number()[spin]);
  for (int nex = 1, nd = 1; nex < abij.maximum_excitation_number()[spin]; nex++)
  {
    // expanding some of them by hand for efficiency
    if (nex == 1)
    {
      for (auto it = abij.unique_begin(1)[spin]; it < abij.unique_end(1)[spin]; ++it, ++nd)
        if (nd % ngrp == rank)
        {
          ov[nd] = T[*((*it) + 1)][*(*it)];
        }
    }
    else if (nex == 2)
    {
      for (auto it = abij.unique_begin(2)[spin]; it < abij.unique_end(2)[spin]; ++it, ++nd)
        if (nd % ngrp == rank)
        {
          auto e = *it;
          ov[nd] = ma::D2x2(T[e[2]][e[0]], T[e[2]][e[1]], T[e[3]][e[0]], T[e[3]][e[1]]);
        }
    }
    else if (nex == 3)
    {
      for (auto it = abij.unique_begin(3)[spin]; it < abij.unique_end(3)[spin]; ++it, ++nd)
        if (nd % ngrp == rank)
        {
          auto e = *it;
          ov[nd] = ma::D3x3(T[e[3]][e[0]], T[e[3]][e[1]], T[e[3]][e[2]], T[e[4]][e[0]], T[e[4]][e[1]], T[e[4]][e[2]],
                            T[e[5]][e[0]], T[e[5]][e[1]], T[e[5]][e[2]]);
        }
    }
    else if (nex == 4)
    {
      for (auto it = abij.unique_begin(4)[spin]; it < abij.unique_end(4)[spin]; ++it, ++nd)
        if (nd % ngrp == rank)
        {
          auto e = *it;
          ov[nd] = ma::D4x4(T[e[4]][e[0]], T[e[4]][e[1]], T[e[4]][e[2]], T[e[4]][e[3]], T[e[5]][e[0]], T[e[5]][e[1]],
                            T[e[5]][e[2]], T[e[5]][e[3]], T[e[6]][e[0]], T[e[6]][e[1]], T[e[6]][e[2]], T[e[6]][e[3]],
                            T[e[7]][e[0]], T[e[7]][e[1]], T[e[7]][e[2]], T[e[7]][e[3]]);
        }
    }
    else if (nex == 5)
    {
      for (auto it = abij.unique_begin(5)[spin]; it < abij.unique_end(5)[spin]; ++it, ++nd)
        if (nd % ngrp == rank)
        {
          auto e = *it;
          ov[nd] = ma::D5x5(T[e[5]][e[0]], T[e[5]][e[1]], T[e[5]][e[2]], T[e[5]][e[3]], T[e[5]][e[4]], T[e[6]][e[0]],
                            T[e[6]][e[1]], T[e[6]][e[2]], T[e[6]][e[3]], T[e[6]][e[4]], T[e[7]][e[0]], T[e[7]][e[1]],
                            T[e[7]][e[2]], T[e[7]][e[3]], T[e[7]][e[4]], T[e[8]][e[0]], T[e[8]][e[1]], T[e[8]][e[2]],
                            T[e[8]][e[3]], T[e[8]][e[4]], T[e[9]][e[0]], T[e[9]][e[1]], T[e[9]][e[2]], T[e[9]][e[3]],
                            T[e[9]][e[4]]);
        }
    }
    else
    {
      boost::multi::array_ref<ComplexType, 2> Qwork_(Qwork.origin(), {nex, nex});
      boost::multi::array_ref<ComplexType, 1> Qwork2_(Qwork.origin() + Qwork_.num_elements(),
                                                      iextensions<1u>{nex * nex});
      for (auto it = abij.unique_begin(nex)[spin]; it < abij.unique_end(nex)[spin]; ++it, ++nd)
        if (nd % ngrp == rank)
        {
          auto exct = *it;
          for (int p = 0; p < nex; p++)
            for (int q = 0; q < nex; q++)
              Qwork_[p][q] = T[exct[p + nex]][exct[q]];
          ov[nd] = ma::determinant<ComplexType>(Qwork_, IWORK, Qwork2_, 0.0);
        }
    }
  }
}

// using simple round-robin scheme for parallelization within TG_local
// assumes that reference determinant is already on [0]
template<class Array1D, class MatA, class MatB, class MatC, class PH_EXCT, class index_aos>
inline void calculate_R(int rank,
                        int ngrp,
                        int spin,
                        PH_EXCT const& abij,
                        index_aos& couplings,
                        MatA&& T,
                        MatB&& Qwork,
                        Array1D&& ov,
                        ComplexType ov0,
                        MatC& R)
{
  using ma::conj;
  using std::get;
  std::vector<int> IWORK(abij.maximum_excitation_number()[spin]);
  std::vector<ComplexType> WORK(abij.maximum_excitation_number()[spin] * abij.maximum_excitation_number()[spin]);
  auto confgs = abij.configurations_begin();
  auto refc   = abij.reference_configuration(spin);
  for (int i = 0; i < std::get<0>(R.sizes()); i++)
    std::fill_n(R[i].origin(), std::get<1>(R.sizes()), ComplexType(0));
  int NEL = std::get<1>(T.sizes());
  std::vector<int> orbs(NEL);
  ComplexType ov_a;
  // add reference contribution!!!
  if (rank == 0)
  {
    ComplexType w(0.0);
    auto it  = to_address(couplings.values()) + (*couplings.pointers_begin(0));
    auto ite = to_address(couplings.values()) + (*couplings.pointers_end(0));
    if (spin == 0)
      for (; it < ite; ++it)
        w += ma::conj(get<2>(*(confgs + (*it)))) * ov[get<1>(*(confgs + (*it)))];
    else
      for (; it < ite; ++it)
        w += ma::conj(get<2>(*(confgs + (*it)))) * ov[get<0>(*(confgs + (*it)))];
    w *= ov0;
    // Wrong if NAEB < NAEA!!! FIX FIX FIX
    for (int i = 0; i < NEL; ++i)
      R[i][refc[i]] += w;
  }
  for (int nex = 1, nd = 1; nex < abij.maximum_excitation_number()[spin]; nex++)
  {
    boost::multi::array_ref<ComplexType, 2> Q(Qwork.origin(), {nex, nex});
    for (auto it = abij.unique_begin(nex)[spin]; it < abij.unique_end(nex)[spin]; ++it, ++nd)
    {
      if (nd % ngrp == rank)
      {
        auto e = *it;
        abij.get_configuration(spin, nd, orbs);
        if (nex == 1)
        {
          ov_a    = T[*((*it) + 1)][*(*it)];
          Q[0][0] = 1.0 / ov_a;
        }
        else if (nex == 2)
        {
          ov_a = ma::I2x2(T[e[2]][e[0]], T[e[2]][e[1]], T[e[3]][e[0]], T[e[3]][e[1]], Q);
        }
        else if (nex == 3)
        {
          ov_a = ma::I3x3(T[e[3]][e[0]], T[e[3]][e[1]], T[e[3]][e[2]], T[e[4]][e[0]], T[e[4]][e[1]], T[e[4]][e[2]],
                          T[e[5]][e[0]], T[e[5]][e[1]], T[e[5]][e[2]], Q);
        }
        else
        {
          for (int p = 0; p < nex; p++)
            for (int q = 0; q < nex; q++)
              Q[p][q] = T[e[p + nex]][e[q]];
          ov_a = ma::invert<ComplexType>(Q, IWORK, WORK, 0.0);
        }
        ComplexType w(0.0);
        auto it  = to_address(couplings.values()) + (*couplings.pointers_begin(nd));
        auto ite = to_address(couplings.values()) + (*couplings.pointers_end(nd));
        if (spin == 0)
          for (; it < ite; ++it)
            w += ma::conj(get<2>(*(confgs + (*it)))) * ov[get<1>(*(confgs + (*it)))];
        else
          for (; it < ite; ++it)
            w += ma::conj(get<2>(*(confgs + (*it)))) * ov[get<0>(*(confgs + (*it)))];
        w *= ov_a * ov0;
        if (std::abs(w) > 1e-10)
        {
          // add term coming from identity
          for (int i = 0; i < NEL; ++i)
            R[i][orbs[i]] += w;
          for (int p = 0; p < nex; ++p)
          {
            auto Rp = R[e[p]];
            auto Ip = Q[p];
            for (int q = 0; q < nex; ++q)
            {
              auto Ipq = Ip[q];
              auto Tq  = T[e[q + nex]];
              for (int i = 0; i < NEL; ++i)
                Rp[orbs[i]] -= w * Ipq * Tq[i];
              Rp[orbs[e[q]]] += w * Ipq;
            }
          }
        }
      }
    }
  }
}

template<class MatE, class MatO, class MatQ, class MatB, class MatT, class MatP, class index_aos>
void calculate_ph_energies_v1(int spin,
                              int rank,
                              int size,
                              MatE&& E,
                              MatO&& Ov,
                              MatQ const& QQ0,
                              MatB&& Qwork,
                              MatT const& Scu,
                              MatT const& Jcb,
                              MatT const& Xcb,
                              MatP const& cPua,
                              MatT const& G0,
                              MatT&& Guv,
                              MatT const& Muv,
                              MatT&& T1,
                              MatT&& T2,
                              MatT&& T3,
                              MatT&& xn,
                              ph_excitations<int, ComplexType> const& abij,
                              std::array<index_aos, 2> const& det_couplings)
{
  using ma::conj;
  using std::get;
  /*
  std::vector<int> IWORK(abij.maximum_excitation_number()[spin]);
  std::vector<ComplexType> WORK(abij.maximum_excitation_number()[spin]*abij.maximum_excitation_number()[spin]);
  auto confgs = abij.configurations_begin();
  auto refc = abij.reference_configuration(spin);
  for(int i=0; i<R.size(0); i++)
    std::fill_n(R[i].origin(),R.size(1),ComplexType(0));
  int NEL = T.size(1);
  std::vector<int> orbs(NEL);
  ComplexType ov_a;
  // add reference contribution!!!
  if(rank==0){
    ComplexType w(0.0);
    auto it = to_address(couplings.values()) + (*couplings.pointers_begin(0));
    auto ite = to_address(couplings.values()) + (*couplings.pointers_end(0));
    if(spin==0)
      for(; it<ite; ++it)
        w += ma::conj(get<2>(*(confgs+(*it)))) * ov[get<1>(*(confgs+(*it)))];
    else
      for(; it<ite; ++it)
        w += ma::conj(get<2>(*(confgs+(*it)))) * ov[get<0>(*(confgs+(*it)))];
    w *= ov0;
// Wrong if NAEB < NAEA!!! FIX FIX FIX
    for(int i=0; i<NEL; ++i)
      R[i][refc[i]] += w;
  }
  for(int nex = 1, nd=1 ; nex<abij.maximum_excitation_number()[spin]; nex++) {
    boost::multi::array_ref<ComplexType,2> Q(Qwork.origin(),{nex,nex});
    for(auto it = abij.unique_begin(nex)[spin]; it<abij.unique_end(nex)[spin]; ++it, ++nd) {
      if(nd%ngrp==rank) {
        auto e = *it;
        abij.get_configuration(spin,nd,orbs);
        if(nex==1) {
          ov_a = T[*((*it)+1)][*(*it)];
          Q[0][0] = 1.0/ov_a;
        } else if(nex==2) {
          ov_a = ma::I2x2( T[e[2]][e[0]], T[e[2]][e[1]],
                           T[e[3]][e[0]], T[e[3]][e[1]], Q);
        } else if(nex==3) {
          ov_a = ma::I3x3( T[e[3]][e[0]], T[e[3]][e[1]],T[e[3]][e[2]],
                           T[e[4]][e[0]], T[e[4]][e[1]],T[e[4]][e[2]],
                           T[e[5]][e[0]], T[e[5]][e[1]],T[e[5]][e[2]],
                          Q);
        } else {
          for(int p=0; p<nex; p++)
            for(int q=0; q<nex; q++)
              Q[p][q] = T[e[p+nex]][e[q]];
          ov_a = ma::invert(Q,IWORK,WORK,0.0);
        }
        ComplexType w(0.0);
        auto it = to_address(couplings.values()) + (*couplings.pointers_begin(nd));
        auto ite = to_address(couplings.values()) + (*couplings.pointers_end(nd));
        if(spin==0)
          for(; it<ite; ++it)
            w += ma::conj(get<2>(*(confgs+(*it)))) * ov[get<1>(*(confgs+(*it)))];
        else
          for(; it<ite; ++it)
            w += ma::conj(get<2>(*(confgs+(*it)))) * ov[get<0>(*(confgs+(*it)))];
        w *= ov_a*ov0;
        if(std::abs(w) > 1e-10) {
          // add term coming from identity
          for(int i=0; i<NEL; ++i)
            R[i][orbs[i]] += w;
          for(int p=0; p<nex; ++p) {
            auto Rp = R[e[p]];
            auto Ip = Q[p];
            for(int q=0; q<nex; ++q) {
              auto Ipq = Ip[q];
              auto Tq = T[e[q+nex]];
              for(int i=0; i<NEL; ++i)
               Rp[orbs[i]] -= w*Ipq*Tq[i];
              Rp[orbs[e[q]]] += w*Ipq;
            }
          }
        }
      }
    }
  }
*/
}


} // namespace afqmc
} // namespace qmcplusplus


#endif
