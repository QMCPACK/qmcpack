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
#include "AFQMC/Numerics/detail/CUDA/cublas_wrapper.hpp"
#include "AFQMC/Numerics/detail/CUDA/cusolver_wrapper.hpp"
#include "AFQMC/Numerics/detail/CUDA/Kernels/setIdentity.cuh"

namespace qmc_cuda 
{
  using qmcplusplus::afqmc::remove_complex;

  template<typename T, typename R, typename I>
  inline static void hevr (char JOBZ, char RANGE, char UPLO, int N,    
                         cuda_gpu_ptr<T> A, int LDA, T VL, T VU,int IL, int IU, T ABSTOL, int &M,                                                      
                         cuda_gpu_ptr<T> W, cuda_gpu_ptr<T> Z, int LDZ, cuda_gpu_ptr<I> ISUPPZ,                  
                         cuda_gpu_ptr<T> WORK, int &LWORK,      
                         cuda_gpu_ptr<R> RWORK, int &LRWORK,               
                         cuda_gpu_ptr<I> IWORK, int &LIWORK, int& INFO)
  {
    throw std::runtime_error("Error: hevr not implemented in gpu."); 
  }

  // getrf_bufferSize
  template<typename T>
  inline static void getrf_bufferSize (const int n, const int m, cuda_gpu_ptr<T> a, int lda, int& lwork)
  {
    cusolver::cusolver_getrf_bufferSize(*a.handles.cusolverDn_handle,n, m, to_address(a), lda, &lwork);
  }

  template<typename T, typename R, typename I>
  inline static void getrf (const int n, const int m, cuda_gpu_ptr<T> && a, int lda, 
                            cuda_gpu_ptr<I> && piv, int &st, cuda_gpu_ptr<R> work) 
  {
    cusolverStatus_t status = cusolver::cusolver_getrf(*a.handles.cusolverDn_handle, n, m,
                     to_address(a), lda, to_address(work), to_address(piv), to_address(piv)+n);
    cudaMemcpy(&st,to_address(piv)+n,sizeof(int),cudaMemcpyDeviceToHost);
    if(CUSOLVER_STATUS_SUCCESS != status) { 
      std::cerr<<" cublas_getrf status, info: " <<status <<" " <<st <<std::endl; std::cerr.flush();
      throw std::runtime_error("Error: cublas_getrf returned error code."); 
    }
  }

  // getrfBatched
  template<typename T, typename I>
  inline static void getrfBatched (const int n, cuda_gpu_ptr<T> * a, int lda, cuda_gpu_ptr<I> piv, cuda_gpu_ptr<I> info, int batchSize)
  {
    T **A_d;
    T **A_h;
    A_h = new T*[batchSize];
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
  template<typename T>
  inline static void getri_bufferSize (int n, cuda_gpu_ptr<T> a, int lda, int& lwork)
  {
    // gpu uses getrs to invert matrix, which requires n*n workspace 
    lwork = n*n;
  }

  // write separate query function to avoid hack!!!
  template<typename T, typename R, typename I>
  inline static void getri(int n, cuda_gpu_ptr<T> a, int lda, cuda_gpu_ptr<I> piv, cuda_gpu_ptr<R> work, int n1, int& status)
  {
    if(n1 < n*n)
      throw std::runtime_error("Error: getri<GPU_MEMORY_POINTER_TYPE> required buffer space of n*n."); 
    if(lda != n)
      throw std::runtime_error("Error: getri<GPU_MEMORY_POINTER_TYPE> required lda = 1."); 

    int* info;
    if(cudaSuccess != cudaMalloc ((void**)&info,sizeof(int))) {
      std::cerr<<" Error getri: Error allocating on GPU." <<std::endl;
      throw std::runtime_error("Error: cudaMalloc returned error code.");
    }

    kernels::set_identity(n,n,to_address(work),n);
    if(CUSOLVER_STATUS_SUCCESS != cusolver::cusolver_getrs(*a.handles.cusolverDn_handle, CUBLAS_OP_N, n, n,
                   to_address(a), lda, to_address(piv), to_address(work), n, info))    
      throw std::runtime_error("Error: cusolver_getrs returned error code."); 
    cudaMemcpy(to_address(a),to_address(work),n*n*sizeof(T),cudaMemcpyDeviceToDevice);
    cudaMemcpy(&status,info,sizeof(int),cudaMemcpyDeviceToHost);
    cudaFree(info);

  }

  // getriBatched
  template<typename T, typename I>
  inline static void getriBatched (int n, cuda_gpu_ptr<T> * a, int lda, cuda_gpu_ptr<I> piv, cuda_gpu_ptr<T> * ainv, int ldc, cuda_gpu_ptr<I> info, int batchSize)
  {
    T **A_d, **C_d;
    T **A_h, **C_h;
    A_h = new T*[batchSize];
    C_h = new T*[batchSize];
    for(int i=0; i<batchSize; i++) {
      A_h[i] = to_address(a[i]);
      C_h[i] = to_address(ainv[i]);
    }
    cudaMalloc((void **)&A_d,  batchSize*sizeof(*A_h));
    cudaMalloc((void **)&C_d,  batchSize*sizeof(*C_h));
    cudaMemcpy(A_d, A_h, batchSize*sizeof(*A_h), cudaMemcpyHostToDevice);
    cudaMemcpy(C_d, C_h, batchSize*sizeof(*C_h), cudaMemcpyHostToDevice);
    cublasStatus_t status = cublas::cublas_getriBatched(*(a[0]).handles.cublas_handle, n, A_d, lda,
                                                        to_address(piv), C_d, ldc, to_address(info), batchSize);
    if(CUBLAS_STATUS_SUCCESS != status)
      throw std::runtime_error("Error: cublas_getri returned error code.");
    cudaFree(A_d);
    cudaFree(C_d);
    delete [] A_h;
    delete [] C_h;
  }

  // matinveBatched
  template<typename T1, typename T2, typename I>
  inline static void matinvBatched (int n, cuda_gpu_ptr<T1> * a, int lda, cuda_gpu_ptr<T2> * ainv, int lda_inv, cuda_gpu_ptr<I> info, int batchSize)
  {
    T1 **A_d, **A_h;
    T2 **C_d, **C_h;
    A_h = new T1*[batchSize];
    C_h = new T2*[batchSize];
    for(int i=0; i<batchSize; i++) {
      A_h[i] = to_address(a[i]);
      C_h[i] = to_address(ainv[i]);
    }
    cudaMalloc((void **)&A_d,  batchSize*sizeof(*A_h));
    cudaMalloc((void **)&C_d,  batchSize*sizeof(*C_h));
    cudaMemcpy(A_d, A_h, batchSize*sizeof(*A_h), cudaMemcpyHostToDevice);
    cudaMemcpy(C_d, C_h, batchSize*sizeof(*C_h), cudaMemcpyHostToDevice);
    cublasStatus_t status = cublas::cublas_matinvBatched(*(a[0]).handles.cublas_handle, n, 
                    A_d, lda, C_d, lda_inv, to_address(info), batchSize);
    if(CUBLAS_STATUS_SUCCESS != status)
      throw std::runtime_error("Error: cublas_matinv returned error code.");
    cudaFree(A_d);
    cudaFree(C_d);
    delete [] A_h;
    delete [] C_h;
  }

  // geqrf
  template<typename T>
  inline static void geqrf_bufferSize (int m, int n, cuda_gpu_ptr<T> a, int lda, int& lwork)
  {
    if(CUSOLVER_STATUS_SUCCESS != cusolver::cusolver_geqrf_bufferSize(*a.handles.cusolverDn_handle,
                m, n, to_address(a), lda, &lwork))
      throw std::runtime_error("Error: cusolver_geqrf_bufferSize returned error code.");
  }

  template<typename T>
  inline static void geqrf(int M, int N, cuda_gpu_ptr<T> A, const int LDA, cuda_gpu_ptr<T> TAU, cuda_gpu_ptr<T> WORK, int LWORK, int& INFO) 
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
      std::cerr<<" cublas_geqrf status, info: " <<status <<" " <<INFO <<std::endl; std::cerr.flush();
      throw std::runtime_error("Error: cublas_geqrf returned error code.");
    }
    cudaFree(piv);
  }

  // gelqf
  template<typename T>
  inline static void gelqf_bufferSize (int m, int n, cuda_gpu_ptr<T> a, int lda, int& lwork)
  {
      lwork = 0;  
  }

  template<typename T>
  inline static void gelqf(int M, int N, cuda_gpu_ptr<T> A, const int LDA, cuda_gpu_ptr<T> TAU, cuda_gpu_ptr<T> WORK, int LWORK, int& INFO)
  {
      throw std::runtime_error("Error: gelqf not implemented in CUDA backend. \n"); 
  }

 // gqr
  template<typename T>
  static void gqr_bufferSize (int m, int n, int k, cuda_gpu_ptr<T> a, int lda, int& lwork)
  {
    if(CUSOLVER_STATUS_SUCCESS != cusolver::cusolver_gqr_bufferSize(*a.handles.cusolverDn_handle,
                                            m,n,k,to_address(a),lda,&lwork))
      throw std::runtime_error("Error: cusolver_gqr_bufferSize returned error code.");
  }

  template<typename T>
  void static gqr(int M, int N, int K, cuda_gpu_ptr<T> A, const int LDA, cuda_gpu_ptr<T> TAU, cuda_gpu_ptr<T> WORK, int LWORK, int& INFO)
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
      std::cerr<<" cublas_gqr status, info: " <<status <<" " <<INFO <<std::endl; std::cerr.flush();
      throw std::runtime_error("Error: cublas_gqr returned error code.");
    }
    cudaFree(piv);
  }

  template<typename T, typename I>
  void static gqrStrided(int M, int N, int K, cuda_gpu_ptr<T> A, const int LDA, const int Astride, cuda_gpu_ptr<T> TAU, const int Tstride, cuda_gpu_ptr<T> WORK, int LWORK, cuda_gpu_ptr<I> info, int batchSize )
  {
    cusolverStatus_t status = cusolver::cusolver_gqr_strided(*A.handles.cusolverDn_handle, M, N, K,
                    to_address(A), LDA, Astride, to_address(TAU), Tstride, to_address(WORK), LWORK, 
                    to_address(info), batchSize);
    if(CUSOLVER_STATUS_SUCCESS != status) {
      std::cerr<<" cublas_gqr_strided status: " <<status <<std::endl; std::cerr.flush();
      throw std::runtime_error("Error: cublas_gqr_strided returned error code.");
    }
  }

  // glq 
  template<typename T>
  static void glq_bufferSize (int m, int n, int k, cuda_gpu_ptr<T> a, int lda, int& lwork)
  {
      lwork = 0;  
  }

  template<typename T>
  void static glq(int M, int N, int K, cuda_gpu_ptr<T> A, const int LDA, cuda_gpu_ptr<T> TAU, cuda_gpu_ptr<T> WORK, int LWORK, int& INFO)
  {
      throw std::runtime_error("Error: glq not implemented in CUDA backend. \n"); 
  }

  // batched geqrf
  template<typename T, typename I>
  inline static void geqrfBatched(int M, int N, cuda_gpu_ptr<T> *A, const int LDA, cuda_gpu_ptr<T>* TAU, cuda_gpu_ptr<I> info, int batchSize)
  {
    T **B_h = new T*[2*batchSize];
    T **A_h(B_h);
    T **T_h(B_h+batchSize);
    for(int i=0; i<batchSize; i++) 
      A_h[i] = to_address(A[i]);
    for(int i=0; i<batchSize; i++) 
      T_h[i] = to_address(TAU[i]);
    T **B_d;
    std::vector<int> inf(batchSize);
    cudaMalloc((void **)&B_d,  2*batchSize*sizeof(*B_h));
    cudaMemcpy(B_d, B_h, 2*batchSize*sizeof(*B_h), cudaMemcpyHostToDevice);
    T **A_d(B_d);
    T **T_d(B_d+batchSize);
    cublasStatus_t status = cublas::cublas_geqrfBatched(*(A[0]).handles.cublas_handle, M, N, A_d, LDA,
                                                        T_d, to_address(inf.data()), batchSize);
    if(CUBLAS_STATUS_SUCCESS != status)
      throw std::runtime_error("Error: cublas_geqrfBatched returned error code.");
    cudaFree(B_d);
    delete [] B_h;
  }

  template<typename T, typename I>
  inline static void geqrfStrided(int M, int N, cuda_gpu_ptr<T> A, const int LDA, const int Astride, cuda_gpu_ptr<T> TAU, const int Tstride, cuda_gpu_ptr<I> info, int batchSize)
  {
/*
    T **B_h = new T*[2*batchSize];
    T **A_h(B_h);
    T **T_h(B_h+batchSize);
    for(int i=0; i<batchSize; i++)
      A_h[i] = to_address(A)+i*Astride;
    for(int i=0; i<batchSize; i++)
      T_h[i] = to_address(TAU)+i*Tstride;
    T **B_d;
    cudaMalloc((void **)&B_d,  2*batchSize*sizeof(*B_h));
    cudaMemcpy(B_d, B_h, 2*batchSize*sizeof(*B_h), cudaMemcpyHostToDevice);
    T **A_d(B_d);
    T **T_d(B_d+batchSize);
*/

    std::vector<int> inf(batchSize);
    T **A_h = new T*[batchSize];
    T **T_h = new T*[batchSize];
    for(int i=0; i<batchSize; i++)
      A_h[i] = to_address(A)+i*Astride;
    for(int i=0; i<batchSize; i++)
      T_h[i] = to_address(TAU)+i*Tstride;
    T **A_d, **T_d;
    cudaMalloc((void **)&A_d,  batchSize*sizeof(*A_h));
    cudaMemcpy(A_d, A_h, batchSize*sizeof(*A_h), cudaMemcpyHostToDevice);
    cudaMalloc((void **)&T_d,  batchSize*sizeof(*T_h));
    cudaMemcpy(T_d, T_h, batchSize*sizeof(*T_h), cudaMemcpyHostToDevice);
    cublasStatus_t status = cublas::cublas_geqrfBatched(*A.handles.cublas_handle, M, N, A_d, LDA,
                                                        T_d, to_address(inf.data()), batchSize);
    for(int i=0; i<batchSize; i++)
      assert(inf[i]==0);
    if(CUBLAS_STATUS_SUCCESS != status)
      throw std::runtime_error("Error: cublas_geqrfBatched returned error code.");
    cudaFree(A_d);
    delete [] A_h;
    cudaFree(T_d);
    delete [] T_h;
  }

  // gesvd_bufferSize
  template<typename T>
  inline static void gesvd_bufferSize (const int m, const int n, cuda_gpu_ptr<T> a, int& lwork)
  {
    cusolver::cusolver_gesvd_bufferSize(*a.handles.cusolverDn_handle, m, n, to_address(a), &lwork);
  }

  template<typename T, typename R>
  inline static void gesvd (char jobU, char jobVT, const int m, const int n, 
                            cuda_gpu_ptr<T> && A, int lda, cuda_gpu_ptr<R> && S,
                            cuda_gpu_ptr<T> && U, int ldu,
                            cuda_gpu_ptr<T> && VT, int ldvt,
                            cuda_gpu_ptr<T> && W, int lw,
                            int &st)
  {
    int *devSt;
    cudaMalloc((void **)&devSt,  sizeof(int));
    cusolverStatus_t status = cusolver::cusolver_gesvd(*A.handles.cusolverDn_handle, 
                    jobU, jobVT, m, n, to_address(A), lda, to_address(S), 
                    to_address(U), ldu, to_address(VT), ldvt, to_address(W), lw,  
                    devSt);
    cudaMemcpy(&st,devSt,sizeof(int),cudaMemcpyDeviceToHost);
    if(CUSOLVER_STATUS_SUCCESS != status) {
      std::cerr<<" cublas_gesvd status, info: " <<status <<" " <<st <<std::endl; std::cerr.flush();
      throw std::runtime_error("Error: cublas_gesvd returned error code.");
    }
    cudaFree(devSt);
  }

  template<typename T, typename R>
  inline static void gesvd (char jobU, char jobVT, const int m, const int n,
                            cuda_gpu_ptr<T> && A, int lda, cuda_gpu_ptr<R> && S,
                            cuda_gpu_ptr<T> && U, int ldu,
                            cuda_gpu_ptr<T> && VT, int ldvt,
                            cuda_gpu_ptr<T> && W, int lw,
                            cuda_gpu_ptr<R> && RW,
                            int &st)
  {
    int *devSt;
    cudaMalloc((void **)&devSt,  sizeof(int));
    cusolverStatus_t status = cusolver::cusolver_gesvd(*A.handles.cusolverDn_handle,
                    jobU, jobVT, m, n, to_address(A), lda, to_address(S),
                    to_address(U), ldu, to_address(VT), ldvt, to_address(W), lw,
                    devSt);
    cudaMemcpy(&st,devSt,sizeof(int),cudaMemcpyDeviceToHost);
    if(CUSOLVER_STATUS_SUCCESS != status) {
      std::cerr<<" cublas_gesvd status, info: " <<status <<" " <<st <<std::endl; std::cerr.flush();
      throw std::runtime_error("Error: cublas_gesvd returned error code.");
    }
    cudaFree(devSt);
  }
  

}

#endif
