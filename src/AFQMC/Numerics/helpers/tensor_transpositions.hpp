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

#include<cassert>
#if defined(QMC_CUDA)
#include "AFQMC/Memory/CUDA/cuda_gpu_pointer.hpp"
#include "AFQMC/Kernels/KaKjw_to_KKwaj.cuh"
#include "AFQMC/Kernels/KaKjw_to_QKajw.cuh"
#include "AFQMC/Kernels/vKKwij_to_vwKiKj.cuh"
#endif

namespace ma
{

// Need OpenMP!!!

template<typename Q, typename T>
void KaKjw_to_KKwaj( int nwalk, int nkpts, int nmo_max, int nmo_tot,
                     int nocc_max, int* nopk, int* nopk0,
                     int* nelpk, int* nelpk0,
                     Q const* A, T*  B) 
{
// OpenMP: Combine Ka,Kj loops into single loop and call parallel for
  int naj=nocc_max*nmo_max;
  int na0 = 0;
  for(int Ka=0; Ka<nkpts; Ka++) {
    int na = nelpk[Ka];
    int nj0=0;
    for(int Kj=0; Kj<nkpts; Kj++) {
      int nj = nopk[Kj];
      //auto G_(to_address(GKK[0][Ka][Kj].origin()));
      auto G_( B + (Ka*nkpts+Kj)*nwalk*nocc_max*nmo_max);
      for(int a=0; a<na; a++) {
        //auto Gc_( to_address(Gca[na0+a][nj0].origin()) );
        auto Gc_( A + (na0+a)*nmo_tot*nwalk + nj0*nwalk); 
        int aj = a*nmo_max;
        for(int j=0; j<nj; j++, aj++) {
          for(int w=0, waj=0; w<nwalk; w++, ++Gc_, waj+=naj)
            G_[waj+aj] = static_cast<T>(*Gc_);
        }
      }
      nj0 += nj;
    }
    na0 += na;
  }
}

template<typename T, typename T1>
void KaKjw_to_QKajw( int nwalk, int nkpts, int nmo_max, int nmo_tot,
                     int nocc_max, int* nmo, int* nmo0,
                     int* nocc, int* nocc0, int* QKtok2,
                     T1 const* A, T*  B) 
{
// OpenMP: Combine Q,K loops into single loop and call parallel for
  for(int Q=0; Q<nkpts; Q++) {
    for(int K=0; K<nkpts; K++) {
      int Ka = K;
      int Kj = QKtok2[Q*nkpts+Ka];
      int na = nocc[Ka];
      int nj = nmo[Kj];
      int na0=nocc0[Ka];
      int nj0=nmo0[Kj];
      //auto G_(to_address(GKK[Q][K].origin()));
      auto G_( B + (Q*nkpts+K)*nwalk*nocc_max*nmo_max);
      for(int a=0, a0=0; a<na; a++, a0+=nmo_max*nwalk) {
        //auto Gc_( to_address(Gca[na0+a][nj0].origin()) );
        auto Gc_( A + (na0+a)*nmo_tot*nwalk + nj0*nwalk);
        for(int j=0, aj=a0; j<nj; j++, aj+=nwalk) {
          for(int w=0; w<nwalk; w++, ++Gc_)
            G_[aj + w] = static_cast<T>(*Gc_);
        }
      }
    }
  }
}

template<typename T, typename Q>
void vKKwij_to_vwKiKj( int nwalk, int nkpts, int nmo_max, int nmo_tot,
                     bool* kk, int* nopk, int* nopk0,
                     Q const* A, T*  B) 
{
  for(int w=0; w<nwalk; w++) {
    for(int Ki=0; Ki<(nkpts+1); Ki++) {
      for(int Kj=0; Kj<nkpts; Kj++) {
        int Ki_ = (Ki==nkpts?Kj:Ki); // Ki==nkpts is the second term of Q0
        int ni = nopk[Ki_];
        int nj = nopk[Kj];
        int ni0 = nopk0[Ki_];
        int nj0 = nopk0[Kj];
        if(kk[Ki*nkpts+Kj]) { // transpose
          auto vb_( B + w*nmo_tot*nmo_tot + ni0*nmo_tot + nj0);
          auto v_( A + ((Ki*nkpts+Kj)*nwalk + w )*nmo_max*nmo_max);
          for(int i=0; i<ni; i++) 
            for(int j=0; j<nj; j++)
              vb_[i*nmo_tot + j] += static_cast<T>(v_[j*nmo_max+i]);
        } else { // copy
          for(int i=0; i<ni; i++) {
            auto vb_( B + w*nmo_tot*nmo_tot + (ni0+i)*nmo_tot + nj0);
            auto v_( A + (((Ki*nkpts+Kj)*nwalk + w )*nmo_max + i)*nmo_max);
            for(int j=0; j<nj; j++) 
              vb_[j] += static_cast<T>(v_[j]);
          }
        }
      }  
    }
  }
}

#ifdef QMC_CUDA

template<typename T, typename Q>
void KaKjw_to_KKwaj( int nwalk, int nkpts, int nmo_max, int nmo_tot,
                     int nocc_max, cuda_gpu_ptr<int> nmo, cuda_gpu_ptr<int> nmo0,
                     cuda_gpu_ptr<int> nocc, cuda_gpu_ptr<int> nocc0,
                     cuda_gpu_ptr<Q> A, cuda_gpu_ptr<T>  B) {
  kernels::KaKjw_to_KKwaj(nwalk,nkpts,nmo_max,nmo_tot,nocc_max,to_address(nmo),to_address(nmo0),
                          to_address(nocc),to_address(nocc0),to_address(A),to_address(B));
}

template<typename T, typename Q>
void KaKjw_to_QKajw( int nwalk, int nkpts, int nmo_max, int nmo_tot,
                     int nocc_max, cuda_gpu_ptr<int> nmo, cuda_gpu_ptr<int> nmo0,
                     cuda_gpu_ptr<int> nocc, cuda_gpu_ptr<int> nocc0, cuda_gpu_ptr<int> QKtok2,
                     cuda_gpu_ptr<Q> A, cuda_gpu_ptr<T>  B) {
  kernels::KaKjw_to_QKajw(nwalk,nkpts,nmo_max,nmo_tot,nocc_max,to_address(nmo),to_address(nmo0),
              to_address(nocc),to_address(nocc0),to_address(QKtok2),to_address(A),to_address(B));
}

template<typename T, typename Q>
void vKKwij_to_vwKiKj( int nwalk, int nkpts, int nmo_max, int nmo_tot,
                     cuda_gpu_ptr<bool> kk, cuda_gpu_ptr<int> nmo, cuda_gpu_ptr<int> nmo0,
                     cuda_gpu_ptr<Q> A, cuda_gpu_ptr<T>  B) {
  kernels::vKKwij_to_vwKiKj(nwalk,nkpts,nmo_max,nmo_tot,to_address(kk),to_address(nmo),
                            to_address(nmo0),to_address(A),to_address(B));
}

#endif

}


#endif
