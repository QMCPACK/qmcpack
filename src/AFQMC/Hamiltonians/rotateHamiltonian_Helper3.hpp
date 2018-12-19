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

#ifndef QMCPLUSPLUS_AFQMC_ROTATEHAMILTONIAN_HELPER3_H
#define QMCPLUSPLUS_AFQMC_ROTATEHAMILTONIAN_HELPER3_H

#include<vector>
#include<tuple>
#include<mpi.h>

#include "Numerics/OhmmsBlas.h"
#include <Utilities/FairDivide.h>
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
// Ta(ka,lb) = Q(k,a,:)*Q(l,b,:)
template<class MattQk,
         class MatQk,
         class MatTa>
inline void count_Qk_x_Qk_symmetric(WALKER_TYPES walker_type, TaskGroup_& TG, std::vector<std::size_t>& sz, int k0, int kN, int l0, int lN, int NMO, int NAEA, int NAEB, MattQk const& tQk, MatQk const& Qk, MatTa && Ta, const RealType cut)
{
  using Type = typename std::decay<MatTa>::type::element;
  assert(tQk.shape()[0] == Ta.shape()[0]);
  assert(tQk.shape()[1] == Qk.shape()[1]);
  assert(Qk.shape()[0] == Ta.shape()[1]);
  int nvec = tQk.shape()[1];
  int nk = kN-k0;
  int norb = lN-l0;
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();

  bool amIAlpha = true;
  if( l0 >= NMO && lN >= NMO )
    amIAlpha = false;
  if(amIAlpha && (k0 > NMO || kN > NMO))
    APP_ABORT(" Error: Block error, this should not happen. \n\n\n");
  if(not amIAlpha && (k0 < NMO || kN < NMO))
    APP_ABORT(" Error: Beta block error, this should not happen. \n\n\n");

  int bl0=-1, blN=-1;
  int nwork = std::min(int(Qk.shape()[0]),ncores);
  if(coreid <nwork)
    std::tie(bl0, blN) = FairDivideBoundary(coreid,int(Qk.shape()[0]),nwork);

  using ma::T;
  if(walker_type == CLOSED) {

    // <a,b | k,l> = Ta(ka,lb) - Tb(la,kb)
    // Ta(ka,lb) = Q(k,a,:)*Q(l,b,:)
    ma::product(tQk,T(Qk[indices[range_t(bl0,blN)][range_t()]]),
                   Ta[indices[range_t()][range_t(bl0,blN)]]);
    TG.node_barrier();
    Type four = Type(4.0);
    Type two = Type(2.0);
    for(int k=k0, ka=0; k<kN; k++) {
      for(int a=0; a<NAEA; a++, ka++) {
        for(int lb=bl0; lb<blN; lb++) { // lb = (l-l0)*NAEA+b
          int b = lb%NAEA;
          if(b<a) continue;
          int l = lb/NAEA+l0;
          int la = (l-l0)*NAEA+a;
          int kb = (k-k0)*NAEA+b;
          Type qkalb = Ta[ka][lb]; // Ta(ka,lb)
          Type qlakb = Ta[kb][la]; // Ta(kb,la)
          //if(std::abs( four*qkalb - two*qlakb ) > cut)
          if(std::abs( two*qlakb ) > cut)
            if(a!=b || k<=l)
              ++sz[a*NMO+k];
        }
      }
    }
  } else if(walker_type == COLLINEAR) {

    if(amIAlpha) {

      // <a,b | k,l> = Ta(ka,lb) - Tb(la,kb)
      // Ta(ka,lb) = Q(k,a,:)*R(l,b,:)
      ma::product(tQk,T(Qk[indices[range_t(bl0,blN)][range_t()]]),
                 Ta[indices[range_t()][range_t(bl0,blN)]]);
      TG.node_barrier();
      for(int k=k0, ka=0; k<kN; k++) {
        for(int a=0; a<NAEA; a++, ka++) {
          for(int lb=bl0; lb<blN; lb++) { // lb = (l-l0)*NAEA+b
            int b = lb%NAEA;
            if(b<a) continue;
            int l = lb/NAEA+l0;
            int la = (l-l0)*NAEA+a;
            int kb = (k-k0)*NAEA+b;
            Type qkalb = Ta[ka][lb];  // Ta(ka,lb)
            Type qlakb = Ta[kb][la];  // Ta(kb,la)
            //if(std::abs( qkalb - qlakb ) > cut)
            if(std::abs( qlakb ) > cut)
              if(a!=b || k<=l)
                ++sz[a*NMO+k];
          }
        }
      }

    } else {

      // k is beta
      // <a,b | k,l> = Ta(ka,lb) - Tb(la,kb)
      // Ta(ka,lb) = Q(k,a,:)*R(l,b,:)
      ma::product(tQk,T(Qk[indices[range_t(bl0,blN)][range_t()]]),
                 Ta[indices[range_t()][range_t(bl0,blN)]]);
      TG.node_barrier();
      for(int k=k0, ka=0; k<kN; k++) {
        for(int a=0; a<NAEB; a++, ka++) {
          for(int lb=bl0; lb<blN; lb++) { // lb = (l-l0)*NAEB+b
            int b = lb%NAEB;
            if(b<a) continue;
            int l = lb/NAEB+l0;
            int la = (l-l0)*NAEB+a;
            int kb = (k-k0)*NAEB+b;
            Type qkalb = Ta[ka][lb]; // Ta(ka,lb)
            Type qlakb = Ta[kb][la]; // Ta(kb,la)
            //if(std::abs( qkalb - qlakb ) > cut)
            if(std::abs( qlakb ) > cut)
              if(a!=b || k<=l)
                ++sz[NAEA*NMO+a*NMO+k-NMO];
          }
        }
      }
    }

  } else if(walker_type == NONCOLLINEAR) {
    APP_ABORT(" Error in count_Qk_x_Qk: GHF not implemented. \n\n\n");
  }

}

// <a,b | k,l> = Ta(ka,lb) - Tb(la,kb)
// Ta(ka,lb) = Q(k,a,:)*Q(l,b,:)
template<class MattQk,
         class MatQk,
         class MatTa,
         class Container>
inline void Qk_x_Qk_symmetric(WALKER_TYPES walker_type, TaskGroup_& TG, int k0, int kN, int l0, int lN, int NMO, int NAEA, int NAEB, MattQk const& tQk, MatQk const& Qk, MatTa && Ta, Container& Vijkl, const RealType cut)
{
  using Type = typename std::decay<MatTa>::type::element;
  assert(tQk.shape()[0] == Ta.shape()[0]);
  assert(tQk.shape()[1] == Qk.shape()[1]);
  assert(Qk.shape()[0] == Ta.shape()[1]);
  int nvec = tQk.shape()[1];
  int nk = kN-k0;
  int norb = lN-l0;
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();

  bool amIAlpha = true;
  if( l0 >= NMO && lN >= NMO )
    amIAlpha = false;
  if(amIAlpha && (k0 > NMO || kN > NMO))
    APP_ABORT(" Error: Block error, this should not happen. \n\n\n");
  if(not amIAlpha && (k0 < NMO || kN < NMO))
    APP_ABORT(" Error: Beta block error, this should not happen. \n\n\n");

  int bl0=-1, blN=-1;
  int ka0=-1, kaN=-1;
  int nwork = std::min(int(Qk.shape()[0]),ncores);
  if(coreid <nwork)
    std::tie(bl0, blN) = FairDivideBoundary(coreid,int(Qk.shape()[0]),nwork);
  nwork = std::min(int(tQk.shape()[0]),ncores);
  if(coreid <nwork)
    std::tie(ka0, kaN) = FairDivideBoundary(coreid,int(tQk.shape()[0]),nwork);

  using ma::T;
  if(walker_type == CLOSED) {

    // <a,b | k,l> = Ta(ka,lb) - Tb(la,kb)
    // Ta(ka,lb) = Q(k,a,:)*Q(l,b,:)
    ma::product(tQk,T(Qk[indices[range_t(bl0,blN)][range_t()]]),
                   Ta[indices[range_t()][range_t(bl0,blN)]]);
    TG.node_barrier();
    Type four = Type(4.0);
    Type two = Type(2.0);
    for(int ka=ka0; ka<kaN; ka++) {  // ka = local range index
      int k = ka/NAEA+k0;  // global index
      int a = ka%NAEA;     // global index
      for(int b=a; b<NAEA; b++) {
        for(int l=l0; l<lN; l++) {
          int la = (l-l0)*NAEA+a;
          int lb = (l-l0)*NAEA+b;
          int kb = (k-k0)*NAEA+b;
          Type qkalb = Ta[ka][lb];  // Ta(ka,lb)
          Type qlakb = Ta[kb][la];  // Ta(kb,la)
          //if(std::abs( four*qkalb - two*qlakb ) > cut)
          if(std::abs( two*qlakb ) > cut)
            if(a!=b || k<l)
              emplace(Vijkl,std::forward_as_tuple(a*NMO+k, b*NMO+l, 2.0*(-two*qlakb)));
            else if(k==l)
              emplace(Vijkl,std::forward_as_tuple(a*NMO+k, b*NMO+l, (-two*qlakb)));
        }
      }
    }
  } else if(walker_type == COLLINEAR) {

    if(amIAlpha) {

      // <a,b | k,l> = Ta(ka,lb) - Tb(la,kb)
      // Ta(ka,lb) = Q(k,a,:)*R(l,b,:)
      ma::product(tQk,T(Qk[indices[range_t(bl0,blN)][range_t()]]),
                 Ta[indices[range_t()][range_t(bl0,blN)]]);
      TG.node_barrier();
      for(int ka=ka0; ka<kaN; ka++) {  // ka = local range index
        int k = ka/NAEA+k0;  // global index
        int a = ka%NAEA;     // global index
        for(int b=a; b<NAEA; b++) {
          for(int l=l0; l<lN; l++) {
            int la = (l-l0)*NAEA+a;
            int lb = (l-l0)*NAEA+b;
            int kb = (k-k0)*NAEA+b;
            Type qkalb = Ta[ka][lb];  // Ta(ka,lb)
            Type qlakb = Ta[kb][la];  // Ta(kb,la)
            //if(std::abs( qkalb - qlakb ) > cut)
            if(std::abs( qlakb ) > cut)
              if(a!=b || k<l)
                emplace(Vijkl, std::forward_as_tuple(a*NMO+k, b*NMO+l, 2.0*(-qlakb)) );
              else if(k==l)
                emplace(Vijkl, std::forward_as_tuple(a*NMO+k, b*NMO+l, (-qlakb)) );
          }
        }
      }

    } else {

      // k is beta
      // <a,b | k,l> = Ta(ka,lb) - Tb(la,kb)
      // Ta(ka,lb) = Q(k,a,:)*R(l,b,:)
      ma::product(tQk,T(Qk[indices[range_t(bl0,blN)][range_t()]]),
                 Ta[indices[range_t()][range_t(bl0,blN)]]);
      TG.node_barrier();
      for(int ka=ka0; ka<kaN; ka++) {  // ka = local range index
        int k = ka/NAEB+k0;  // global index
        int a = ka%NAEB;     // global index
        for(int b=a; b<NAEB; b++) {
          for(int l=l0; l<lN; l++) {
            int la = (l-l0)*NAEB+a;
            int lb = (l-l0)*NAEB+b;
            int kb = (k-k0)*NAEB+b;
            Type qkalb = Ta[ka][lb];  // Ta(ka,lb)
            Type qlakb = Ta[kb][la];  // Ta(kb,la)
            //if(std::abs( qkalb - qlakb ) > cut)
            if(std::abs( qlakb ) > cut)
              if(a!=b || k<l)
                emplace(Vijkl, std::forward_as_tuple(NMO*NAEA+a*NMO+k-NMO, NMO*NAEA+b*NMO+l-NMO, 2.0*(-qlakb)) );
              else if(k==l)
                emplace(Vijkl, std::forward_as_tuple(NMO*NAEA+a*NMO+k-NMO, NMO*NAEA+b*NMO+l-NMO, (-qlakb)) );
          }
        }
      }
    }

  } else if(walker_type == NONCOLLINEAR) {
    APP_ABORT(" Error in count_Qk_x_Qk: GHF not implemented. \n\n\n");
  }

}

}

}

#endif


