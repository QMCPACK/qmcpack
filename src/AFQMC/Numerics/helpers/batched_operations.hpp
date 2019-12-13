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
#if defined(ENABLE_CUDA)
#include "AFQMC/Memory/CUDA/cuda_gpu_pointer.hpp"
#include "AFQMC/Numerics/detail/CUDA/Kernels/batched_dot_wabn_wban.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/dot_wabn.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/batched_Tab_to_Klr.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/Tab_to_Kl.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/vbias_from_v1.cuh"
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
void batched_dot_wanb_wbna( int nbatch, int nwalk, int nocc, int nchol,
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
          E_ += ma::dot(nchol,A_+a*nocc*nchol+b,nocc,B_+b*nocc*nchol+a,nocc);
      y[w*incy] += static_cast<std::complex<T>>(alpha[batch]*E_);
    }
  }
}

template<typename T, typename Q>
void dot_wabn( int nwalk, int nocc, int nchol,
                    std::complex<Q> alpha, std::complex<Q> const* Tab,
                    std::complex<T>* y, int incy)
{
  int nocc2nc = nocc*nocc*nchol;
  for(int w=0; w<nwalk; ++w) {
    std::complex<Q> E_(0.0);
    auto A_(Tab + w*nocc2nc);
    using ma::dot;
    for(int a=0; a<nocc; ++a)
      for(int b=0; b<nocc; ++b)
        E_ += ma::dot(nchol,A_+(a*nocc+b)*nchol,1,A_+(b*nocc+a)*nchol,1);
    y[w*incy] += static_cast<std::complex<T>>(alpha*E_);
  }
}

template<typename T, typename Q>
void dot_wanb( int nwalk, int nocc, int nchol,
                    std::complex<Q> alpha, std::complex<Q> const* Tab,
                    std::complex<T>* y, int incy)
{
  int nocc2nc = nocc*nchol*nocc;
  for(int w=0; w<nwalk; ++w) {
    std::complex<Q> E_(0.0);
    auto A_(Tab + w*nocc2nc);
    using ma::dot;
    for(int a=0; a<nocc; ++a)
      for(int b=0; b<nocc; ++b)
        E_ += ma::dot(nchol,A_+a*nocc*nchol+b,nocc,A_+b*nocc*nchol+a,nocc);
    y[w*incy] += static_cast<std::complex<T>>(alpha*E_);
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

template<typename T, typename Q>
void batched_Tanb_to_Klr(int nterms, int nwalk, int nocc, int nchol_max,
                    int nchol_tot, int ncholQ, int ncholQ0, int* kdiag,
                    Q const* Tab, T*  Kl, T* Kr)
{
  for(int w=0; w<nwalk; ++w) {
    for( int k=0; k<nterms; k++) {
      int batch = kdiag[k];
      for(int a=0; a<nocc; a++) {
        auto Tba_(Tab + batch*nwalk*nocc*nocc*nchol_max
                                           + ((w*nocc+a)*nocc)*nchol_max+a);
        auto Kr_(Kr + w*nchol_tot + ncholQ0);
        for(int c=0; c<ncholQ; ++c)
          Kr_[c] += static_cast<T>(Tba_[c*nocc]);
      }
    }
    for( int k=0; k<nterms; k++) {
      int batch = kdiag[k];
      for(int a=0; a<nocc; a++) {
        auto Tab_(Tab + (batch+1)*nwalk*nocc*nocc*nchol_max
                                      + ((w*nocc+a)*nocc)*nchol_max+a);
        auto Kl_(Kl + w*nchol_tot + ncholQ0);
        for(int c=0; c<ncholQ; ++c)
          Kl_[c] += static_cast<T>(Tab_[c*nocc]);
      }
    }
  }
}

template<typename T, typename Q>
void Tab_to_Kl(int nwalk, int nocc, int nchol, Q const* Tab, T*  Kl) 
{
  for(int w=0; w<nwalk; ++w) {
    for(int a=0; a<nocc; a++) {
      auto Tab_(Tab + ((w*nocc+a)*nocc + a)*nchol);
      auto Kl_(Kl + w*nchol);
      for(int c=0; c<nchol; ++c)
        Kl_[c] += static_cast<T>(Tab_[c]);
    }
  }
}

template<typename T, typename Q>
void Tanb_to_Kl(int nwalk, int nocc, int nchol, int nchol_tot, Q const* Tab, T*  Kl)
{
  for(int w=0; w<nwalk; ++w) {
    for(int a=0; a<nocc; a++) {
      auto Tab_(Tab + ((w*nocc+a)*nocc)*nchol+a);
      auto Kl_(Kl + w*nchol_tot);
      for(int c=0; c<nchol; ++c)
        Kl_[c] += static_cast<T>(Tab_[c*nocc]);
    }
  }
}

template<typename T, typename T1>
void vbias_from_v1( int nwalk, int nkpts, int nchol_max, int* Qsym, int* kminus,
                    int* ncholpQ, int* ncholpQ0, std::complex<T> const alpha,
                    std::complex<T1> const* v1, std::complex<T>* vb)
{
  for(int Q=0; Q<nkpts; Q++) {
    if( Qsym[Q] < 0 ) return;
    int Qm = kminus[Q];
    int nc0 = ncholpQ0[Q];
    int nc = ncholpQ[Q];
    int ncm = ncholpQ[Qm];
    int Qm_ = Qm;
    int ntot = nc*nwalk;
    if( Qsym[Q] > 0 ) Qm_ = nkpts+Qsym[Q]-1;

    // v+
    auto vb_(vb + nc0*nwalk);
    auto v1_(v1 + Q*nchol_max*nwalk );
    auto v2_(v1 + Qm_*nchol_max*nwalk );
    // v+ = a*(v[Q]+v[-Q]) 
    for(int n=0; n<ntot; ++n)
      vb_[ n ] += alpha*static_cast<std::complex<T>>(v1_[ n ]);
    for(int n=0; n<ntot; ++n)
      vb_[ n ] += alpha*static_cast<std::complex<T>>(v2_[ n ]);
    // v-
    vb_ = (vb + (nc0+nc)*nwalk);
    // v- = -a*i*(v[Q]-v[-Q]) 
    auto ialpha(alpha*std::complex<double>(0.0,1.0));
    for(int n=0; n<ntot; ++n)
      vb_[ n ] -= ialpha*static_cast<std::complex<T>>(v1_[ n ]);
    for(int n=0; n<ntot; ++n)
      vb_[ n ] += ialpha*static_cast<std::complex<T>>(v2_[ n ]);
  }
}

// for n in [0,N), y[incy*n] = beta * y[incy*n] + alpha sum_m^{0,M} opA(A)[n,m] * opB(B)[n,m]  
template<typename T, typename Q>
void batched_dot( char TA, char TB, int N, int M, std::complex<T> const alpha, 
                  std::complex<Q> const* A, int lda, std::complex<Q> const* B, int ldb,
                  std::complex<T> const beta, std::complex<T>* y, int incy)
{
  bool cA( TA == 'H' || TA == 'C');  
  bool cB( TB == 'H' || TB == 'C');  
  bool tA( TA == 'H' || TA == 'T');  
  bool tB( TB == 'H' || TB == 'T');  
  if(not tA && not TB) {
    for(int n=0; n<N; n++) {
      std::complex<T> r(0.0,0.0);
      auto an(A+n*lda); 
      auto bn(B+n*ldb); 
      if(cA && cB) {
        for(int m=0; m<M; m++, an++, bn++) r += std::conj(*an) * std::conj(*bn);
      } else if(cA && not cB) {  
        for(int m=0; m<M; m++, an++, bn++) r += std::conj(*an) * (*bn);
      } else if(not cA && cB) {  
        for(int m=0; m<M; m++, an++, bn++) r += (*an) * std::conj(*bn);
      } else { 
        for(int m=0; m<M; m++, an++, bn++) r += (*an) * (*bn);
      }  
      y[incy*n] = beta * y[incy*n] + alpha * r;
    }  
  } else if(tA && not TB) {
    for(int n=0; n<N; n++) {  
      std::complex<T> r(0.0,0.0);
      auto an(A+n); 
      auto bn(B+n*ldb);
      if(cA && cB) {
        for(int m=0; m<M; m++, an+=lda, bn++) r += std::conj(*an) * std::conj(*bn);
      } else if(cA && not cB) {
        for(int m=0; m<M; m++, an+=lda, bn++) r += std::conj(*an) * (*bn);
      } else if(not cA && cB) {
        for(int m=0; m<M; m++, an+=lda, bn++) r += (*an) * std::conj(*bn);
      } else { 
        for(int m=0; m<M; m++, an+=lda, bn++) r += (*an) * (*bn);
      }
      y[incy*n] = beta * y[incy*n] + alpha * r;
    }
  } else if(not tA && TB) {
    for(int n=0; n<N; n++) {
      std::complex<T> r(0.0,0.0);
      auto an(A+n*lda);    
      auto bn(B+n);
      if(cA && cB) {
        for(int m=0; m<M; m++, an++, bn+=ldb) r += std::conj(*an) * std::conj(*bn);
      } else if(cA && not cB) {
        for(int m=0; m<M; m++, an++, bn+=ldb) r += std::conj(*an) * (*bn);
      } else if(not cA && cB) {
        for(int m=0; m<M; m++, an++, bn+=ldb) r += (*an) * std::conj(*bn);
      } else { 
        for(int m=0; m<M; m++, an++, bn+=ldb) r += (*an) * (*bn);
      }
      y[incy*n] = beta * y[incy*n] + alpha * r;
    } 
  } else {  // special case, tA && tB
    for(int n=0; n<N; n++) y[incy*n] *= beta; 
    for(int m=0; m<M; m++) { 
      auto am(A+m*lda);
      auto bm(B+m*ldb);
      if(cA && cB) {
        for(int n=0; n<N; n++, am++, bm++) y[incy*n] += alpha * std::conj(*am) * std::conj(*bm);      
      } else if(cA && not cB) {
        for(int n=0; n<N; n++, am++, bm++) y[incy*n] += alpha * std::conj(*am) * (*bm);      
      } else if(not cA && cB) {
        for(int n=0; n<N; n++, am++, bm++) y[incy*n] += alpha * (*am) * std::conj(*bm);      
      } else { 
        for(int n=0; n<N; n++, am++, bm++) y[incy*n] += alpha * (*am) * (*bm);      
      }
    }
  } 
}


} // namespace ma

#ifdef ENABLE_CUDA
namespace qmc_cuda{

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

template<typename T, typename Q>
void batched_Tanb_to_Klr(int nterms, int nwalk, int nocc, int nchol_max,
                    int nchol_tot, int ncholQ, int ncholQ0, cuda_gpu_ptr<int> kdiag,
                    cuda_gpu_ptr<Q> Tab, cuda_gpu_ptr<T>  Kl,
                    cuda_gpu_ptr<T> Kr)
{
  kernels::batched_Tanb_to_Klr(nterms,nwalk,nocc,nchol_max,nchol_tot,ncholQ,ncholQ0,
                             to_address(kdiag), to_address(Tab),
                             to_address(Kl), to_address(Kr));
}

template<typename T, typename Q>
void Tab_to_Kl(int nwalk, int nocc, int nchol,
                    cuda_gpu_ptr<Q> Tab, cuda_gpu_ptr<T>  Kl)
{
  kernels::Tab_to_Kl(nwalk,nocc,nchol,to_address(Tab),to_address(Kl));
}

template<typename T, typename Q>
void Tanb_to_Kl(int nwalk, int nocc, int nchol, int nchol_tot,
                    cuda_gpu_ptr<Q> Tab, cuda_gpu_ptr<T>  Kl)
{
  kernels::Tanb_to_Kl(nwalk,nocc,nchol,nchol_tot,to_address(Tab),to_address(Kl));
}

template<typename T, typename Q, typename R>
void batched_dot_wabn_wban( int nbatch, int nwalk, int nocc, int nchol,
                    cuda_gpu_ptr<R> alpha, cuda_gpu_ptr<Q> Tab,
                    T* y , int incy)
{
  kernels::batched_dot_wabn_wban(nbatch,nwalk,nocc,nchol,to_address(alpha),to_address(Tab),
                                 y,incy);
}

template<typename T, typename Q, typename R>
void batched_dot_wanb_wbna( int nbatch, int nwalk, int nocc, int nchol,
                    cuda_gpu_ptr<R> alpha, cuda_gpu_ptr<Q> Tab,
                    T* y , int incy)
{
  kernels::batched_dot_wanb_wbna(nbatch,nwalk,nocc,nchol,to_address(alpha),to_address(Tab),
                                 y,incy);
}

template<typename T, typename Q, typename R>
void dot_wabn( int nwalk, int nocc, int nchol, R alpha, cuda_gpu_ptr<Q> Tab,
                    T* y , int incy)
{
  kernels::dot_wabn(nwalk,nocc,nchol,alpha,to_address(Tab),y,incy);
}

template<typename T, typename Q, typename R>
void dot_wanb( int nwalk, int nocc, int nchol, R alpha, cuda_gpu_ptr<Q> Tab,
                    T* y , int incy)
{
  kernels::dot_wanb(nwalk,nocc,nchol,alpha,to_address(Tab),y,incy);
}

template<typename T, typename Q, typename R>
void vbias_from_v1( int nwalk, int nkpts, int nchol_max, cuda_gpu_ptr<int> Qsym, cuda_gpu_ptr<int> kminus,
                    cuda_gpu_ptr<int> ncholpQ, cuda_gpu_ptr<int> ncholpQ0, R alpha,
                    cuda_gpu_ptr<Q> v1, T* vb)
{
  kernels::vbias_from_v1(nwalk,nkpts,nchol_max,to_address(Qsym),to_address(kminus),
            to_address(ncholpQ),to_address(ncholpQ0),alpha,to_address(v1),vb);
}

template<typename T, typename Q>
void batched_dot( char TA, char TB, int N, int M, T alpha,
                  cuda_gpu_ptr<Q> A, int lda, cuda_gpu_ptr<Q> B, int ldb,
                  T beta, cuda_gpu_ptr<T> y, int incy)
{
//  kernels::batched_dot(TA,TB,N,M,alpha,to_address(A),lda,to_address(B),ldb,beta,to_address(y),incy);
    APP_ABORT(" Error: batched_dot not yet available in gpu.\n");
}

}
#endif



#endif
