//////////////////////////////////////////////////////////////////////////////
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

#ifndef QMCPLUSPLUS_AFQMC_ROTATEHAMILTONIAN_HELPER2_H
#define QMCPLUSPLUS_AFQMC_ROTATEHAMILTONIAN_HELPER2_H

#include <vector>
#include <tuple>
#include <mpi.h>

#include "Utilities/FairDivide.h"
#include "AFQMC/config.h"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Matrix/csr_matrix.hpp"
#include "AFQMC/Utilities/afqmc_TTI.hpp"

namespace qmcplusplus
{
namespace afqmc
{
// <a,b | k,l> = Ta(ka,lb) - Tb(la,kb)
// Ta(ka,lb) = Q(k,a,:)*R(l,b,:)
template<class MatQk, class MatRl, class MatTa>
inline void count_Qk_x_Rl(WALKER_TYPES walker_type,
                          SPComplexType EJX,
                          TaskGroup_& TG,
                          std::vector<std::size_t>& sz,
                          int k0,
                          int kN,
                          int l0,
                          int lN,
                          int NMO,
                          int NAEA,
                          int NAEB,
                          MatQk const& Qk,
                          MatRl const& Rl,
                          MatTa&& Ta,
                          const SPRealType cut)
{
  using Type = typename std::decay<MatTa>::type::element;
  assert(std::get<0>(Qk.sizes()) == std::get<0>(Ta.sizes()));
  assert(std::get<1>(Qk.sizes()) == std::get<0>(Rl.sizes()));
  assert(std::get<1>(Rl.sizes()) == std::get<1>(Rl.sizes()));
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();

  bool amIAlpha = true;
  if (l0 >= NMO && lN >= NMO)
    amIAlpha = false;

  int bl0 = -1, blN = -1;
  int nwork = std::min(int(std::get<1>(Rl.sizes())), ncores);
  if (coreid < nwork)
    std::tie(bl0, blN) = FairDivideBoundary(coreid, int(std::get<1>(Rl.sizes())), nwork);
  int ka0 = -1, kaN = -1;
  nwork = std::min(int(std::get<0>(Qk.sizes())), ncores);
  if (coreid < nwork)
    std::tie(ka0, kaN) = FairDivideBoundary(coreid, int(std::get<0>(Qk.sizes())), nwork);

  Type four(4.0);
  Type two(2.0);
  if (walker_type == CLOSED)
  {
    // <a,b | k,l> = Ta(ka,lb) - Tb(la,kb)
    // Ta(ka,lb) = Q(k,a,:)*R(l,b,:)
    ma::product(Qk, Rl(Rl.extension(0), {bl0, blN}), Ta(Ta.extension(0), {bl0, blN}));
    TG.node_barrier();
    for (int ka = ka0; ka < kaN; ka++)
    {                         // ka = local range index
      int k = ka / NAEA + k0; // global index
      int a = ka % NAEA;      // global index
      for (int b = a; b < NAEA; b++)
      {
        for (int l = l0; l < lN; l++)
        {
          int la     = (l - l0) * NAEA + a;
          int lb     = (l - l0) * NAEA + b;
          int kb     = (k - k0) * NAEA + b;
          Type qkalb = Ta[ka][lb]; // Ta(ka,lb)
          Type qlakb = Ta[kb][la]; // Ta(kb,la)
          if (std::abs(EJX * four * qkalb - two * qlakb) > cut)
            if (a != b || k <= l)
              ++sz[a * NMO + k];
        }
      }
    }
  }
  else if (walker_type == COLLINEAR)
  {
    if (k0 < NMO && (kN - 1) < NMO)
    {
      // k is alpha

      if (amIAlpha)
      {
        // <a,b | k,l> = Ta(ka,lb) - Tb(la,kb)
        // Ta(ka,lb) = Q(k,a,:)*R(l,b,:)
        ma::product(Qk, Rl(Rl.extension(0), {bl0, blN}), Ta(Ta.extension(0), {bl0, blN}));
        TG.node_barrier();
        for (int k = k0, ka = 0; k < kN; k++)
        {
          for (int a = 0; a < NAEA; a++, ka++)
          {
            for (int lb = bl0; lb < blN; lb++)
            { // lb = (l-l0)*NAEA+b
              int b = lb % NAEA;
              if (b < a)
                continue;
              int l      = lb / NAEA + l0;
              int la     = (l - l0) * NAEA + a;
              int kb     = (k - k0) * NAEA + b;
              Type qkalb = Ta[ka][lb]; // Ta(ka,lb)
              Type qlakb = Ta[kb][la]; // Ta(kb,la)
              if (std::abs(EJX * qkalb - qlakb) > cut)
                if (a != b || k <= l)
                  ++sz[a * NMO + k];
            }
          }
        }
      }
      else if (std::abs(EJX) > 1e-8)
      {
        // <a,b | k,l> = Ta(ka,lb) = Q(k,a,:)*R(l,b,:)
        ma::product(Qk, Rl(Rl.extension(0), {bl0, blN}), Ta(Ta.extension(0), {bl0, blN}));
        TG.node_barrier();
        for (int k = k0, ka = 0; k < kN; k++)
        {
          for (int a = 0; a < NAEA; a++, ka++)
          {
            for (int lb = bl0; lb < blN; lb++)
            {                          // lb = (l-l0)*NAEB+b
              Type qkalb = Ta[ka][lb]; // Ta(ka,lb)
              if (std::abs(EJX * qkalb) > cut)
                ++sz[a * NMO + k];
            }
          }
        }
      }
    }
    else if (k0 >= NMO && kN >= NMO)
    {
      // k is beta
      if (!amIAlpha)
      {
        // <a,b | k,l> = Ta(ka,lb) - Tb(la,kb)
        // Ta(ka,lb) = Q(k,a,:)*R(l,b,:)
        ma::product(Qk, Rl(Rl.extension(0), {bl0, blN}), Ta(Ta.extension(0), {bl0, blN}));
        TG.node_barrier();
        for (int k = k0, ka = 0; k < kN; k++)
        {
          for (int a = 0; a < NAEB; a++, ka++)
          {
            for (int lb = bl0; lb < blN; lb++)
            { // lb = (l-l0)*NAEB+b
              int b = lb % NAEB;
              if (b < a)
                continue;
              int l      = lb / NAEB + l0;
              int la     = (l - l0) * NAEB + a;
              int kb     = (k - k0) * NAEB + b;
              Type qkalb = Ta[ka][lb]; // Ta(ka,lb)
              Type qlakb = Ta[kb][la]; // Ta(kb,la)
              if (std::abs(EJX * qkalb - qlakb) > cut)
                if (a != b || k <= l)
                  ++sz[NAEA * NMO + a * NMO + k - NMO];
            }
          }
        }
      }
    }
    else
    {
      APP_ABORT(" Error: This should not happen. \n\n\n");
    }
  }
  else if (walker_type == NONCOLLINEAR)
  {
    APP_ABORT(" Error in count_Qk_x_Rl: GHF not implemented. \n\n\n");
  }
}

// <a,b | k,l> = Ta(ka,lb) - Tb(la,kb)
// Ta(ka,lb) = Q(k,a,:)*R(l,b,:)
template<class MatQk, class MatRl, class MatTa, class Container>
inline void Qk_x_Rl(WALKER_TYPES walker_type,
                    SPComplexType EJX,
                    TaskGroup_& TG,
                    int k0,
                    int kN,
                    int l0,
                    int lN,
                    int NMO,
                    int NAEA,
                    int NAEB,
                    MatQk const& Qk,
                    MatRl const& Rl,
                    MatTa&& Ta,
                    Container& Vijkl,
                    const SPRealType cut)
{
  using Type = typename std::decay<MatTa>::type::element;
  assert(std::get<0>(Qk.sizes()) == std::get<0>(Ta.sizes()));
  assert(std::get<1>(Qk.sizes()) == std::get<0>(Rl.sizes()));
  assert(std::get<1>(Rl.sizes()) == std::get<1>(Rl.sizes()));
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();

  bool amIAlpha = true;
  if (l0 >= NMO && lN >= NMO)
    amIAlpha = false;

  int bl0 = -1, blN = -1;
  int ka0 = -1, kaN = -1;
  int nwork = std::min(int(std::get<1>(Rl.sizes())), ncores);
  if (coreid < nwork)
    std::tie(bl0, blN) = FairDivideBoundary(coreid, int(std::get<1>(Rl.sizes())), nwork);
  nwork = std::min(int(std::get<0>(Qk.sizes())), ncores);
  if (coreid < nwork)
    std::tie(ka0, kaN) = FairDivideBoundary(coreid, int(std::get<0>(Qk.sizes())), nwork);

  Type four(4.0);
  Type two(2.0);
  if (walker_type == CLOSED)
  {
    // <a,b | k,l> = Ta(ka,lb) - Tb(la,kb)
    // Ta(ka,lb) = Q(k,a,:)*R(l,b,:)
    ma::product(Qk, Rl(Rl.extension(0), {bl0, blN}), Ta(Ta.extension(0), {bl0, blN}));
    TG.node_barrier();
    for (int ka = ka0; ka < kaN; ka++)
    {                         // ka = local range index
      int k = ka / NAEA + k0; // global index
      int a = ka % NAEA;      // global index
      for (int b = a; b < NAEA; b++)
      {
        for (int l = l0; l < lN; l++)
        {
          int la     = (l - l0) * NAEA + a;
          int lb     = (l - l0) * NAEA + b;
          int kb     = (k - k0) * NAEA + b;
          Type qkalb = Ta[ka][lb]; // Ta(ka,lb)
          Type qlakb = Ta[kb][la]; // Ta(kb,la)
          if (std::abs(EJX * four * qkalb - two * qlakb) > cut)
          {
            if (a != b || k < l)
            {
              emplace(Vijkl, std::forward_as_tuple(a * NMO + k, b * NMO + l, two * (EJX * four * qkalb - two * qlakb)));
            }
            else if (k == l)
            {
              emplace(Vijkl, std::forward_as_tuple(a * NMO + k, b * NMO + l, (EJX * four * qkalb - two * qlakb)));
            }
          }
        }
      }
    }
  }
  else if (walker_type == COLLINEAR)
  {
    if (k0 < NMO && (kN - 1) < NMO)
    {
      // k is alpha

      if (amIAlpha)
      {
        // <a,b | k,l> = Ta(ka,lb) - Tb(la,kb)
        // Ta(ka,lb) = Q(k,a,:)*R(l,b,:)
        ma::product(Qk, Rl(Rl.extension(0), {bl0, blN}), Ta(Ta.extension(0), {bl0, blN}));
        TG.node_barrier();
        for (int ka = ka0; ka < kaN; ka++)
        {                         // ka = local range index
          int k = ka / NAEA + k0; // global index
          int a = ka % NAEA;      // global index
          for (int b = a; b < NAEA; b++)
          {
            for (int l = l0; l < lN; l++)
            {
              int la     = (l - l0) * NAEA + a;
              int lb     = (l - l0) * NAEA + b;
              int kb     = (k - k0) * NAEA + b;
              Type qkalb = Ta[ka][lb]; // Ta(ka,lb)
              Type qlakb = Ta[kb][la]; // Ta(kb,la)
              if (std::abs(EJX * qkalb - qlakb) > cut)
              {
                if (a != b || k < l)
                {
                  emplace(Vijkl, std::forward_as_tuple(a * NMO + k, b * NMO + l, two * (EJX * qkalb - qlakb)));
                }
                else if (k == l)
                {
                  emplace(Vijkl, std::forward_as_tuple(a * NMO + k, b * NMO + l, EJX * qkalb - qlakb));
                }
              }
            }
          }
        }
      }
      else if (std::abs(EJX) > 1e-8)
      {
        // <a,b | k,l> = Ta(ka,lb) = Q(k,a,:)*R(l,b,:)
        ma::product(Qk, Rl(Rl.extension(0), {bl0, blN}), Ta(Ta.extension(0), {bl0, blN}));
        TG.node_barrier();
        for (int ka = ka0; ka < kaN; ka++)
        {                         // ka = local range index
          int k = ka / NAEA + k0; // global index
          int a = ka % NAEA;      // global index
          for (int b = 0; b < NAEB; b++)
          {
            for (int l = l0; l < lN; l++)
            {
              int lb     = (l - l0) * NAEB + b;
              Type qkalb = Ta[ka][lb]; // Ta(ka,lb)
              if (std::abs(EJX * qkalb) > cut)
                emplace(Vijkl, std::forward_as_tuple(a * NMO + k, NMO * NAEA + b * NMO + l - NMO, two * (EJX * qkalb)));
            }
          }
        }
      }
    }
    else if (k0 >= NMO && kN >= NMO)
    {
      // k is beta
      if (!amIAlpha)
      {
        // <a,b | k,l> = Ta(ka,lb) - Tb(la,kb)
        // Ta(ka,lb) = Q(k,a,:)*R(l,b,:)
        ma::product(Qk, Rl(Rl.extension(0), {bl0, blN}), Ta(Ta.extension(0), {bl0, blN}));
        TG.node_barrier();
        for (int ka = ka0; ka < kaN; ka++)
        {                         // ka = local range index
          int k = ka / NAEB + k0; // global index
          int a = ka % NAEB;      // global index
          for (int b = a; b < NAEB; b++)
          {
            for (int l = l0; l < lN; l++)
            {
              int la     = (l - l0) * NAEB + a;
              int lb     = (l - l0) * NAEB + b;
              int kb     = (k - k0) * NAEB + b;
              Type qkalb = Ta[ka][lb]; // Ta(ka,lb)
              Type qlakb = Ta[kb][la]; // Ta(kb,la)
              if (std::abs(EJX * qkalb - qlakb) > cut)
              {
                if (a != b || k < l)
                {
                  emplace(Vijkl,
                          std::forward_as_tuple(NMO * NAEA + a * NMO + k - NMO, NMO * NAEA + b * NMO + l - NMO,
                                                two * (EJX * qkalb - qlakb)));
                }
                else if (k == l)
                {
                  emplace(Vijkl,
                          std::forward_as_tuple(NMO * NAEA + a * NMO + k - NMO, NMO * NAEA + b * NMO + l - NMO,
                                                EJX * qkalb - qlakb));
                }
              }
            }
          }
        }
      }
    }
    else
    {
      APP_ABORT(" Error: This should not happen. \n\n\n");
    }
  }
  else if (walker_type == NONCOLLINEAR)
  {
    APP_ABORT(" Error in createHamiltonianForGeneralDeterminant: GHF not implemented. \n\n\n");
  }
}

} // namespace afqmc

} // namespace qmcplusplus

#endif
