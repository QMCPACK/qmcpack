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

#ifndef AFQMC_LAPACK_GPU_HPP
#define AFQMC_LAPACK_GPU_HPP

#include<cassert>
#include "AFQMC/Utilities/type_conversion.hpp"
#include "AFQMC/Memory/custom_pointers.hpp"
#include "AFQMC/Numerics/detail/lapack_cpu.hpp"
#include "AFQMC/Numerics/detail/CUDA/cublas_wrapper.hpp"
#include "AFQMC/Numerics/detail/CUDA/cusolver_wrapper.hpp"
#include "AFQMC/Kernels/setIdentity.cuh"

namespace qmc_cuda 
{
  using qmcplusplus::afqmc::remove_complex;

  // hevr
  template<typename T,
           class ptr,
           class ptrR,
           class ptrI,
           typename = typename std::enable_if_t< (ptr::memory_type != GPU_MEMORY_POINTER_TYPE) or 
                                                 (ptrR::memory_type != GPU_MEMORY_POINTER_TYPE) or 
                                                 (ptrI::memory_type != GPU_MEMORY_POINTER_TYPE) 
                                               >
          >
  inline static void hevr (char JOBZ, char RANGE, char UPLO, int N, 
                         ptr A, int LDA, T VL, T VU,int IL, int IU, T ABSTOL, int &M, 
                         ptrR W, ptr Z, int LDZ, ptrI ISUPPZ, 
                         ptr WORK, int &LWORK, 
                         ptrR RWORK, int &LRWORK, 
                         ptrI IWORK, int &LIWORK, int& INFO)
  {
    using ma::hevr;
    hevr (JOBZ,RANGE,UPLO,N,to_address(A),LDA,VL,VU,IL,IU,ABSTOL,M,to_address(W),to_address(Z),LDZ,to_address(ISUPPZ),
           to_address(WORK),LWORK,to_address(RWORK),LRWORK,to_address(IWORK),LIWORK,INFO);
  }

  template<typename T,
           class ptr,
           class ptrR,
           class ptrI,
           typename = typename std::enable_if_t< (ptr::memory_type == GPU_MEMORY_POINTER_TYPE) and
                                                 (ptrR::memory_type == GPU_MEMORY_POINTER_TYPE) and 
                                                 (ptrI::memory_type == GPU_MEMORY_POINTER_TYPE)
                                               >,
           typename = void 
          >
  inline static void hevr (char JOBZ, char RANGE, char UPLO, int N,    
                         ptr A, int LDA, T VL, T VU,int IL, int IU, T ABSTOL, int &M,                                                      
                         ptrR W, ptr Z, int LDZ, ptrI ISUPPZ,                  
                         ptr WORK, int &LWORK,      
                         ptrR RWORK, int &LRWORK,               
                         ptrI IWORK, int &LIWORK, int& INFO)
  {
    throw std::runtime_error("Error: hevr not implemented in gpu."); 
  }

  // getrf_bufferSize
  template<class ptr,
           typename = typename std::enable_if_t< (ptr::memory_type != GPU_MEMORY_POINTER_TYPE) > 
          >
  inline static void getrf_bufferSize (const int n, const int m, ptr a, int lda, int& lwork) 
  {
    using ma::getrf_bufferSize;
    getrf_bufferSize(n, m, to_address(a), lda, &lwork);
  }

  template<class ptr,
           typename = typename std::enable_if_t< (ptr::memory_type == GPU_MEMORY_POINTER_TYPE) >,
           typename = void 
          >
  inline static void getrf_bufferSize (const int n, const int m, ptr a, int lda, int& lwork)
  {
    cusolver::cusolver_getrf_bufferSize(*a.handles.cusolverDn_handle,n, m, to_address(a), lda, &lwork);
  }

  // getrf
  template<class ptr,
           class ptrW,  
           class ptrI,
           typename = typename std::enable_if_t< (ptr::memory_type != GPU_MEMORY_POINTER_TYPE) or
                                                 (ptrI::memory_type != GPU_MEMORY_POINTER_TYPE)
            >
          >
  inline static void getrf (const int n, const int m, ptr const a, int lda, ptrI piv, int &st, ptrW work) 
  {
    using ma::getrf;
    getrf(n, m, to_address(a), lda, to_address(piv), st);
  }

  template<class ptr,
           class ptrW,
           class ptrI,
           typename = typename std::enable_if_t< (ptr::memory_type == GPU_MEMORY_POINTER_TYPE) and
                                                 (ptrI::memory_type == GPU_MEMORY_POINTER_TYPE)
                                               >,
           typename = void 
          >
  inline static void getrf (const int n, const int m, ptr && a, int lda, ptrI && piv, int &st, ptrW work) 
  {
    cusolverStatus_t status = cusolver::cusolver_getrf(*a.handles.cusolverDn_handle, n, m,
                                       to_address(a), lda, to_address(work), to_address(piv), to_address(piv)+n);
    if(CUSOLVER_STATUS_SUCCESS != status) { 
      int st;
      cudaMemcpy(&st,to_address(piv)+n,sizeof(int),cudaMemcpyDeviceToHost);
      std::cerr<<" cublas_getrf status, info: " <<status <<" " <<st <<std::endl; std::cerr.flush();
      throw std::runtime_error("Error: cublas_getrf returned error code."); 
    }
    cudaMemcpy(&st,to_address(piv)+n,sizeof(int),cudaMemcpyDeviceToHost);

  }

  // getrfBatched
  template<class ptr,
           class ptrI,
           typename = typename std::enable_if_t< (ptr::memory_type != GPU_MEMORY_POINTER_TYPE) and 
                                                 (ptrI::memory_type != GPU_MEMORY_POINTER_TYPE)
            >
          >
  inline static void getrfBatched (const int n, ptr * a, int lda, ptrI piv, ptrI info, int batchSize)
  {
    using Q = typename ptr::value_type;
    Q **A_h;
    A_h = new Q*[batchSize];
    for(int i=0; i<batchSize; i++)
      A_h[i] = to_address(a[i]);
    using ma::getrfBatched;
    getrfBatched(n, A_h, lda, to_address(piv), to_address(info), batchSize);
    delete [] A_h;
  }

  template<class ptr,
           class ptrI,
           typename = typename std::enable_if_t< (ptr::memory_type == GPU_MEMORY_POINTER_TYPE) and
                                                 (ptrI::memory_type == GPU_MEMORY_POINTER_TYPE) >,
           typename = void 
          >
  inline static void getrfBatched (const int n, ptr * a, int lda, ptrI piv, ptrI info, int batchSize)
  {
    using Q = typename ptr::value_type;
    Q **A_d;
    Q **A_h;
    A_h = new Q*[batchSize];
    for(int i=0; i<batchSize; i++) 
      A_h[i] = to_address(a[i]);
    cudaMalloc((void **)&A_d,  batchSize*sizeof(*A_h));
    cudaMemcpy(A_d, A_h, batchSize*sizeof(*A_h), cudaMemcpyHostToDevice);
    cublasStatus_t status = cublas::cublas_getrfBatched(*(a[0]).handles.cublas_handle, n, A_d, lda, 
                                                        to_address(piv), to_address(info), batchSize); 
    if(CUBLAS_STATUS_SUCCESS != status) 
      throw std::runtime_error("Error: cublas_getrf returned error code.");
    cudaFree(A_d);
    delete [] A_h;
  }

  // getri_bufferSize
  template<class ptr,
           typename = typename std::enable_if_t< (ptr::memory_type != GPU_MEMORY_POINTER_TYPE)
                                                >
          >
  inline static void getri_bufferSize (int n, ptr a, int lda, int& lwork)
  {
    using ma::getri_bufferSize;
    getri_bufferSize(n, to_address(a), lda, lwork);
  }

  template<class ptr,
           typename = typename std::enable_if_t< (ptr::memory_type == GPU_MEMORY_POINTER_TYPE)
                                                >,
           typename = void
          >
  inline static void getri_bufferSize (int n, ptr a, int lda, int& lwork)
  {
    // gpu uses getrs to invert matrix, which requires n*n workspace 
    lwork = n*n;
  }

  // getri: will fail if not called correctly, but removing checks on ptrI and ptrW for now
  template<class ptr,
           class ptrI,
           class ptrW,
           typename = typename std::enable_if_t< (ptr::memory_type != GPU_MEMORY_POINTER_TYPE) > 
          >
  inline static void getri(int n, ptr a, int n0, ptrI piv, ptrW work, int n1, int& status)
  {
    using ma::getri;
    getri(n, to_address(a), n0, to_address(piv), to_address(work), n1, status);
  }

  // write separate query function to avoid hack!!!
  template<class ptr,
           class ptrI,
           class ptrW,
           typename = typename std::enable_if_t< (ptr::memory_type == GPU_MEMORY_POINTER_TYPE)>,
           typename = void
          >
  inline static void getri(int n, ptr a, int lda, ptrI piv, ptrW work, int n1, int& status)
  {
    if(n1 < n*n)
      throw std::runtime_error("Error: getri<GPU_MEMORY_POINTER_TYPE> required buffer space of n*n."); 
    if(lda != n)
      throw std::runtime_error("Error: getri<GPU_MEMORY_POINTER_TYPE> required lda = 1."); 
/*
    using detail::to_address;
    std::vector<typename ptr::value_type> buff(n*n);
    std::vector<typename ptr::value_type> w_(n*n);
    std::vector<int> Ibuff(n+1);
    cudaMemcpy(buff.data(),to_address(a),sizeof(typename ptr::value_type)*n*n,cudaMemcpyDeviceToHost);
    cudaMemcpy(Ibuff.data(),to_address(piv),sizeof(int)*(n+1),cudaMemcpyDeviceToHost);
    using ma::getri;
    getri(n, buff.data(), lda, Ibuff.data(), w_.data(), n1, status);
    cudaMemcpy(to_address(a),buff.data(),sizeof(typename ptr::value_type)*n*n,cudaMemcpyHostToDevice);
*/

    kernels::setIdentity(n,to_address(work),n);
    if(CUSOLVER_STATUS_SUCCESS != cusolver::cusolver_getrs(*a.handles.cusolverDn_handle, CUBLAS_OP_N, n, n,
                   to_address(a), lda, to_address(piv), to_address(work), n, to_address(piv)+n))    
      throw std::runtime_error("Error: cusolver_getrs returned error code."); 
    cudaMemcpy(to_address(a),to_address(work),n*n*sizeof(typename ptr::value_type),cudaMemcpyDeviceToDevice);
    cudaMemcpy(&status,to_address(piv)+n,sizeof(int),cudaMemcpyDeviceToHost);

  }

  // getriBatched
  template<class ptr,
           class ptrI,
           typename = typename std::enable_if_t< (ptr::memory_type != GPU_MEMORY_POINTER_TYPE) and
                                                 (ptrI::memory_type != GPU_MEMORY_POINTER_TYPE)
            >
          >
  inline static void getriBatched (int n, ptr * a, int lda, ptrI piv, ptr * c, int lwork, ptrI info, int batchSize)
  {
    using Q = typename ptr::value_type;
    Q **A_h, **C_h;
    A_h = new Q*[batchSize];
    C_h = new Q*[batchSize];
    for(int i=0; i<batchSize; i++) {
      A_h[i] = to_address(a[i]);
      C_h[i] = to_address(c[i]);
    }
    using ma::getriBatched;
    getriBatched(n, to_address(a), lda, to_address(piv), C_h, lwork, to_address(info), batchSize);
    delete [] A_h;
    delete [] C_h;
  }

  template<class ptr,
           class ptrI,
           typename = typename std::enable_if_t< (ptr::memory_type == GPU_MEMORY_POINTER_TYPE) and
                                                 (ptrI::memory_type == GPU_MEMORY_POINTER_TYPE)
            >,
           typename = void 
          >
  inline static void getriBatched (int n, ptr * a, int lda, ptrI piv, ptr * c, int lwork, ptrI info, int batchSize)
  {
    assert(lda == n);
    assert(lwork >= n*n);
    using Q = typename ptr::value_type;
    Q **A_d, **C_d;
    Q **A_h, **C_h;
    A_h = new Q*[batchSize];
    C_h = new Q*[batchSize];
    for(int i=0; i<batchSize; i++) {
      A_h[i] = to_address(a[i]);
      C_h[i] = to_address(c[i]);
    }
    cudaMalloc((void **)&A_d,  batchSize*sizeof(*A_h));
    cudaMalloc((void **)&C_d,  batchSize*sizeof(*C_h));
    cudaMemcpy(A_d, A_h, batchSize*sizeof(*A_h), cudaMemcpyHostToDevice);
    cudaMemcpy(C_d, C_h, batchSize*sizeof(*C_h), cudaMemcpyHostToDevice);
    cublasStatus_t status = cublas::cublas_getriBatched(*(a[0]).handles.cublas_handle, n, A_d, lda,
                                                        to_address(piv), C_d, n, to_address(info), batchSize);
    if(CUBLAS_STATUS_SUCCESS != status)
      throw std::runtime_error("Error: cublas_getri returned error code.");
    for(int i=0; i<batchSize; i++) {
      cudaMemcpy(A_h[i], C_h[i], n*n*sizeof(Q), cudaMemcpyHostToDevice);
    }
    cudaFree(A_d);
    cudaFree(C_d);
    delete [] A_h;
    delete [] C_h;
  }

  // geqrf
  template<class ptr,
           typename = typename std::enable_if_t< (ptr::memory_type != GPU_MEMORY_POINTER_TYPE) >
          >
  inline static void geqrf_bufferSize (int m, int n, ptr a, int lda, int& lwork)
  {
    using ma::geqrf_bufferSize;
    geqrf_bufferSize(m, n, to_address(a), lda, lwork); 
  }


  template<class ptr,
           typename = typename std::enable_if_t< (ptr::memory_type == GPU_MEMORY_POINTER_TYPE) >,
           typename = void
          >
  inline static void geqrf_bufferSize (int m, int n, ptr a, int lda, int& lwork)
  {
    if(CUSOLVER_STATUS_SUCCESS != cusolver::cusolver_geqrf_bufferSize(*a.handles.cusolverDn_handle,
                m, n, to_address(a), lda, &lwork))
      throw std::runtime_error("Error: cusolver_geqrf_bufferSize returned error code.");
  }

  template<class ptr,
           typename = typename std::enable_if_t< (ptr::memory_type != GPU_MEMORY_POINTER_TYPE) > 
          >
  inline static void geqrf(int M, int N, ptr A, const int LDA, ptr TAU, ptr WORK, int LWORK, int& INFO) 
  {
    using ma::geqrf;
    geqrf(M, N, to_address(A), LDA, to_address(TAU), to_address(WORK), LWORK, INFO);
  }

  template<class ptr,
           typename = typename std::enable_if_t< (ptr::memory_type == GPU_MEMORY_POINTER_TYPE) >,
           typename = void
          >
  inline static void geqrf(int M, int N, ptr A, const int LDA, ptr TAU, ptr WORK, int LWORK, int& INFO) 
  {
    // allocating here for now
    int* piv;
    if(cudaSuccess != cudaMalloc ((void**)&piv,sizeof(int))) {
      std::cerr<<" Error geqrf: Error allocating on GPU." <<std::endl;
      throw std::runtime_error("Error: cudaMalloc returned error code.");
    }
    
    cusolverStatus_t status = cusolver::cusolver_geqrf(*A.handles.cusolverDn_handle, M, N,
                   to_address(A), LDA, to_address(TAU), to_address(WORK), LWORK, piv);
    cudaMemcpy(&INFO,piv,sizeof(int),cudaMemcpyDeviceToHost);
    if(CUSOLVER_STATUS_SUCCESS != status) {
      int st;
      std::cerr<<" cublas_getrf status, info: " <<status <<" " <<INFO <<std::endl; std::cerr.flush();
      throw std::runtime_error("Error: cublas_geqrf returned error code.");
    }
    cudaFree(piv);
  }


  // gelqf
  template<class ptr,
           typename = typename std::enable_if_t< (ptr::memory_type != GPU_MEMORY_POINTER_TYPE) >
          >
  inline static void gelqf_bufferSize (int m, int n, ptr a, int lda, int& lwork)
  {
    using ma::gelqf_bufferSize;
    gelqf_bufferSize(m, n, to_address(a), lda, lwork);
  }


  template<class ptr,
           typename = typename std::enable_if_t< (ptr::memory_type == GPU_MEMORY_POINTER_TYPE) >,
           typename = void
          >
  inline static void gelqf_bufferSize (int m, int n, ptr a, int lda, int& lwork)
  {
      throw std::runtime_error("Error: gelqf_bufferSize not implemented in CUDA backend. \n"); 
  }

  template<class ptr,
           typename = typename std::enable_if_t< (ptr::memory_type != GPU_MEMORY_POINTER_TYPE) >
          >
  inline static void gelqf(int M, int N, ptr A, const int LDA, ptr TAU, ptr WORK, int LWORK, int& INFO)
  {
    using ma::gelqf;
    gelqf(M, N, to_address(A), LDA, to_address(TAU), to_address(WORK), LWORK, INFO);
  }

  template<class ptr,
           typename = typename std::enable_if_t< (ptr::memory_type == GPU_MEMORY_POINTER_TYPE) >,
           typename = void
          >
  inline static void gelqf(int M, int N, ptr A, const int LDA, ptr TAU, ptr WORK, int LWORK, int& INFO)
  {
      throw std::runtime_error("Error: gelqf not implemented in CUDA backend. \n"); 
  }

 // gqr
  template<class ptr,
           typename = typename std::enable_if_t< (ptr::memory_type != GPU_MEMORY_POINTER_TYPE) >
          >
  static void gqr_bufferSize (int m, int n, int k, ptr a, int lda, int& lwork)
  {
    using ma::gqr_bufferSize;
    gqr_bufferSize(m, n, k, to_address(a), lda, lwork);
  }

  template<class ptr,
           typename = typename std::enable_if_t< (ptr::memory_type == GPU_MEMORY_POINTER_TYPE) >,
           typename = void
          >
  static void gqr_bufferSize (int m, int n, int k, ptr a, int lda, int& lwork)
  {
    if(CUSOLVER_STATUS_SUCCESS != cusolver::cusolver_gqr_bufferSize(*a.handles.cusolverDn_handle,
                                            m,n,k,to_address(a),lda,&lwork))
      throw std::runtime_error("Error: cusolver_gqr_bufferSize returned error code.");
  }

  template<class ptr,
           typename = typename std::enable_if_t< (ptr::memory_type != GPU_MEMORY_POINTER_TYPE) >
          >
  void static gqr(int M, int N, int K, ptr A, const int LDA, ptr TAU, ptr WORK, int LWORK, int& INFO)
  {
    using ma::gqr;
    gqr(M, N, K, to_address(A), LDA, to_address(TAU), to_address(WORK), LWORK, INFO);
  }

  template<class ptr,
           typename = typename std::enable_if_t< (ptr::memory_type == GPU_MEMORY_POINTER_TYPE) >,
           typename = void
          >
  void static gqr(int M, int N, int K, ptr A, const int LDA, ptr TAU, ptr WORK, int LWORK, int& INFO)
  {
    // allocating here for now
    int* piv;
    if(cudaSuccess != cudaMalloc ((void**)&piv,sizeof(int))) {
      std::cerr<<" Error gqr: Error allocating on GPU." <<std::endl;
      throw std::runtime_error("Error: cudaMalloc returned error code.");
    }

    cusolverStatus_t status = cusolver::cusolver_gqr(*A.handles.cusolverDn_handle, M, N, K,
                   to_address(A), LDA, to_address(TAU), to_address(WORK), LWORK, piv);
    cudaMemcpy(&INFO,piv,sizeof(int),cudaMemcpyDeviceToHost);
    if(CUSOLVER_STATUS_SUCCESS != status) {
      int st;
      std::cerr<<" cublas_getrf status, info: " <<status <<" " <<INFO <<std::endl; std::cerr.flush();
      throw std::runtime_error("Error: cublas_gqr returned error code.");
    }
    cudaFree(piv);
  }

  // glq 
  template<class ptr,
           typename = typename std::enable_if_t< (ptr::memory_type != GPU_MEMORY_POINTER_TYPE) >
          >
  static void glq_bufferSize (int m, int n, int k, ptr a, int lda, int& lwork)
  {
    using ma::glq_bufferSize;
    glq_bufferSize(m, n, k, to_address(a), lda, lwork);
  }

  template<class ptr,
           typename = typename std::enable_if_t< (ptr::memory_type == GPU_MEMORY_POINTER_TYPE) >,
           typename = void
          >
  static void glq_bufferSize (int m, int n, int k, ptr a, int lda, int& lwork)
  {
      throw std::runtime_error("Error: glq not implemented in CUDA backend. \n"); 
  }

  template<class ptr,
           typename = typename std::enable_if_t< (ptr::memory_type != GPU_MEMORY_POINTER_TYPE) >
          >
  void static glq(int M, int N, int K, ptr A, const int LDA, ptr TAU, ptr WORK, int LWORK, int& INFO)
  {
    using ma::glq;
    glq(M, N, K, to_address(A), LDA, to_address(TAU), to_address(WORK), LWORK, INFO);
  }

  template<class ptr,
           typename = typename std::enable_if_t< (ptr::memory_type == GPU_MEMORY_POINTER_TYPE) >,
           typename = void
          >
  void static glq(int M, int N, int K, ptr A, const int LDA, ptr TAU, ptr WORK, int LWORK, int& INFO)
  {
      throw std::runtime_error("Error: glq not implemented in CUDA backend. \n"); 
  }


}

#endif
