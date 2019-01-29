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

#include<cassert>
#if defined(QMC_CUDA)
#include "AFQMC/Memory/CUDA/cuda_gpu_pointer.hpp"
#include "AFQMC/Kernels/batched_dot_wabn_wban.cuh"
#include "AFQMC/Kernels/batched_Tab_to_Klr.cuh"
#include "AFQMC/Kernels/vbias_from_v1.cuh"
#endif

namespace ma
{

// Tab [nbatch][nwalk][nocc][nocc][nchol]
template<typename T, typename Q>
void batched_dot_wabn_wban( int nbatch, int nwalk, int nocc, int nchol,
                    std::complex<Q> const* alpha, std::complex<Q> const* Tab,
                    std::complex<T>* y, int incy)
{
  int nocc2nc = nocc*nocc*nchol;
  for(int batch=0; batch<nbatch; ++batch) {
    for(int w=0; w<nwalk; ++w) {
      std::complex<Q> E_(0.0);
      auto A_(Tab + (2*batch*nwalk+w)*nocc2nc); 
      auto B_(Tab + ((2*batch+1)*nwalk+w)*nocc2nc);
      using ma::dot;
      for(int a=0; a<nocc; ++a)
        for(int b=0; b<nocc; ++b)
          E_ += ma::dot(nchol,A_+(a*nocc+b)*nchol,1,B_+(b*nocc+a)*nchol,1);
      y[w*incy] += static_cast<std::complex<T>>(alpha[batch]*E_);
    }
  }
}


template<typename T, typename Q>
void batched_Tab_to_Klr(int nterms, int nwalk, int nocc, int nchol_max,
                    int nchol_tot, int ncholQ, int ncholQ0, int* kdiag,
                    Q const* Tab, T*  Kl, T* Kr)
{
  for(int w=0; w<nwalk; ++w) {
    for( int k=0; k<nterms; k++) {
      int batch = kdiag[k];
      for(int a=0; a<nocc; a++) {
        auto Tba_(Tab + batch*nwalk*nocc*nocc*nchol_max
                                           + ((w*nocc+a)*nocc + a)*nchol_max);
        auto Kr_(Kr + w*nchol_tot + ncholQ0);
        for(int c=0; c<ncholQ; ++c)
          Kr_[c] += static_cast<T>(Tba_[c]);
      }
    }
    for( int k=0; k<nterms; k++) {
      int batch = kdiag[k];
      for(int a=0; a<nocc; a++) {
        auto Tab_(Tab + (batch+1)*nwalk*nocc*nocc*nchol_max
                                      + ((w*nocc+a)*nocc + a)*nchol_max);
        auto Kl_(Kl + w*nchol_tot + ncholQ0);
        for(int c=0; c<ncholQ; ++c)
          Kl_[c] += static_cast<T>(Tab_[c]);
      }
    }
  }
}

template<typename T, typename T1>
void vbias_from_v1( int nwalk, int nkpts, int nchol_max, int Q0, int* kminus,
                    int* ncholpQ, int* ncholpQ0, std::complex<T> const alpha,
                    std::complex<T1> const* v1, std::complex<T>* vb)
{
  for(int Q=0; Q<nkpts; Q++) {
    int Qm = kminus[Q];
    int nc0 = ncholpQ0[Q];
    int nc = ncholpQ[Q];
    int ncm = ncholpQ[Qm];
    int Qm_ = Qm;
    if(Q==Q0) Qm_=nkpts;

    // v+
    auto vb_(vb + nc0*nwalk);
    auto v1_(v1 + Q*nchol_max*nwalk );
    auto v2_(v1 + Qm_*nchol_max*nwalk );
    // v+ = a*(v[Q]+v[-Q]) 
    for(int n=0, nw=0; n<nc; ++n)
      for(int w=0; w<nwalk; ++w, nw+=nwalk)
        vb_[ nw+w ] += alpha*static_cast<std::complex<T>>(v1_[ nw+w ]);
    for(int n=0, nw=0; n<ncm; ++n, nw+=nwalk)
      for(int w=0; w<nwalk; ++w)
        vb_[ nw+w ] += alpha*static_cast<std::complex<T>>(v2_[ nw+w ]);
    // v-
    vb_ = (vb + (nc0+nc)*nwalk);
    // v- = -a*i*(v[Q]-v[-Q]) 
    auto ialpha(alpha*std::complex<double>(0.0,1.0));
    for(int n=0, nw=0; n<nc; n++)
      for(int w=0; w<nwalk; w++, nw+=nwalk)
        vb_[ nw+w ] -= ialpha*static_cast<std::complex<T>>(v1_[ nw+w ]);
    for(int n=0, nw=0; n<ncm; n++)
      for(int w=0; w<nwalk; w++, nw+=nwalk)
        vb_[ nw+w ] += ialpha*static_cast<std::complex<T>>(v2_[ nw+w ]);
  }
}

#ifdef QMC_CUDA

template<typename T, typename Q>
void batched_Tab_to_Klr(int nterms, int nwalk, int nocc, int nchol_max,
                    int nchol_tot, int ncholQ, int ncholQ0, cuda_gpu_ptr<int> kdiag,
                    cuda_gpu_ptr<Q> Tab, cuda_gpu_ptr<T>  Kl,
                    cuda_gpu_ptr<T> Kr)
{
  kernels::batched_Tab_to_Klr(nterms,nwalk,nocc,nchol_max,nchol_tot,ncholQ,ncholQ0,
                             to_address(kdiag), to_address(Tab),
                             to_address(Kl), to_address(Kr));
}

template<typename T, typename Q, typename R>
void batched_dot_wabn_wban( int nbatch, int nwalk, int nocc, int nchol,
                    cuda_gpu_ptr<R> alpha, cuda_gpu_ptr<Q> Tab,
                    cuda_gpu_ptr<T> y , int incy)
{
  kernels::batched_dot_wabn_wban(nbatch,nwalk,nocc,nchol,to_address(alpha),to_address(Tab),
                                 to_address(y),incy);

}

template<typename T, typename Q, typename R>
void vbias_from_v1( int nwalk, int nkpts, int nchol_max, int Q0, cuda_gpu_ptr<int> kminus,
                    cuda_gpu_ptr<int> ncholpQ, cuda_gpu_ptr<int> ncholpQ0, R alpha,
                    cuda_gpu_ptr<Q> v1, cuda_gpu_ptr<T> vb)
{
  kernels::vbias_from_v1(nwalk,nkpts,nchol_max,Q0,to_address(kminus),to_address(ncholpQ),
            to_address(ncholpQ0),alpha,to_address(v1),to_address(vb));
}

#endif

}


#endif
