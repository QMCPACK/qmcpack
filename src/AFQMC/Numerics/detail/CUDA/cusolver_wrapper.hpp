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

#ifndef CUSOLVER_FUNCTIONDEFS_H
#define CUSOLVER_FUNCTIONDEFS_H

#include<cassert>
#include <cuda_runtime.h>
#include "cusolverDn.h"
#include "AFQMC/Memory/CUDA/cuda_utilities.hpp"

namespace cusolver {

  using qmc_cuda::cublasOperation;

  inline cusolverStatus_t 
  cusolver_getrf_bufferSize(cusolverDnHandle_t handle,
                      int m,
                      int n,
                      float *A,
                      int lda,
                      int *Lwork ) 
  {
    cusolverStatus_t sucess =
              cusolverDnSgetrf_bufferSize(handle,m,n,A,lda,Lwork); 
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t
  cusolver_getrf_bufferSize(cusolverDnHandle_t handle,
                      int m,
                      int n,
                      std::complex<float> *A,
                      int lda,
                      int *Lwork ) 
  { 
    cusolverStatus_t sucess =
              cusolverDnCgetrf_bufferSize(handle,m,n,
                                          reinterpret_cast<cuComplex *>(A),lda,Lwork);
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t
  cusolver_getrf_bufferSize(cusolverDnHandle_t handle,
                      int m,
                      int n,
                      double *A,
                      int lda,
                      int *Lwork ) 
  { 
    cusolverStatus_t sucess =
              cusolverDnDgetrf_bufferSize(handle,m,n,A,lda,Lwork);
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t
  cusolver_getrf_bufferSize(cusolverDnHandle_t handle,
                      int m,
                      int n,
                      std::complex<double> *A,
                      int lda,
                      int *Lwork ) 
  { 
    cusolverStatus_t sucess =
              cusolverDnZgetrf_bufferSize(handle,m,n,
                                          reinterpret_cast<cuDoubleComplex *>(A),lda,Lwork);
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t cusolver_getrf(cusolverDnHandle_t handle,
                                   int m, 
                                   int n,
                                   float *A,
                                   int lda,
                                   float *Work,
                                   int *devIpiv,
                                   int *devInfo )
  {
    cusolverStatus_t sucess =
              cusolverDnSgetrf(handle,m,n,A,lda,Work,devIpiv,devInfo);
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t cusolver_getrf(cusolverDnHandle_t handle,
                                   int m,
                                   int n,
                                   double *A,
                                   int lda,
                                   double *Work,
                                   int *devIpiv,
                                   int *devInfo )
  {
    cusolverStatus_t sucess =
              cusolverDnDgetrf(handle,m,n,A,lda,Work,devIpiv,devInfo);
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t cusolver_getrf(cusolverDnHandle_t handle,
                                   int m,
                                   int n,
                                   std::complex<float> *A,
                                   int lda,
                                   std::complex<float> *Work,
                                   int *devIpiv,
                                   int *devInfo )
  {
    cusolverStatus_t sucess =
              cusolverDnCgetrf(handle,m,n,
                              reinterpret_cast<cuComplex *>(A),lda,
                              reinterpret_cast<cuComplex *>(Work),devIpiv,devInfo);
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t cusolver_getrf(cusolverDnHandle_t handle,
                                   int m,
                                   int n,
                                   std::complex<double> *A,
                                   int lda,
                                   std::complex<double> *Work,
                                   int *devIpiv,
                                   int *devInfo )
  {
    cusolverStatus_t sucess =
              cusolverDnZgetrf(handle,m,n,
                              reinterpret_cast<cuDoubleComplex *>(A),lda,
                              reinterpret_cast<cuDoubleComplex *>(Work),devIpiv,devInfo);
    cudaDeviceSynchronize ();
    return sucess;
  }


  // getrs
  inline cusolverStatus_t 
  cusolver_getrs(cusolverDnHandle_t handle,
           cublasOperation_t trans,
           int n,
           int nrhs,
           const float *A,
           int lda,
           const int *devIpiv,
           float *B,
           int ldb,
           int *devInfo )
  {
    cusolverStatus_t sucess =
              cusolverDnSgetrs(handle,trans,n,nrhs,A,lda,devIpiv,B,ldb,devInfo);
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t
  cusolver_getrs(cusolverDnHandle_t handle,
           cublasOperation_t trans,
           int n,
           int nrhs,
           const double *A,
           int lda,
           const int *devIpiv,
           double *B,
           int ldb,
           int *devInfo )
  {
    cusolverStatus_t sucess =
              cusolverDnDgetrs(handle,trans,n,nrhs,A,lda,devIpiv,B,ldb,devInfo);
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t
  cusolver_getrs(cusolverDnHandle_t handle,
           cublasOperation_t trans,
           int n,
           int nrhs,
           const std::complex<float> *A,
           int lda,
           const int *devIpiv,
           std::complex<float> *B,
           int ldb,
           int *devInfo )
  {
    cusolverStatus_t sucess =
              cusolverDnCgetrs(handle,trans,n,nrhs,
                               reinterpret_cast<cuComplex const*>(A),lda,devIpiv,
                               reinterpret_cast<cuComplex *>(B),ldb,devInfo);
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t
  cusolver_getrs(cusolverDnHandle_t handle,
           cublasOperation_t trans,
           int n,
           int nrhs,
           const std::complex<double> *A,
           int lda,
           const int *devIpiv,
           std::complex<double> *B,
           int ldb,
           int *devInfo )
  { 
    cusolverStatus_t sucess =
              cusolverDnZgetrs(handle,trans,n,nrhs,
                               reinterpret_cast<cuDoubleComplex const*>(A),lda,devIpiv,
                               reinterpret_cast<cuDoubleComplex *>(B),ldb,devInfo);
    cudaDeviceSynchronize ();
    return sucess;
  }

  //geqrf_bufferSize
  inline cusolverStatus_t 
  cusolver_geqrf_bufferSize(cusolverDnHandle_t handle,
                      int m,
                      int n,
                      float *A,
                      int lda,
                      int *Lwork )
  {
    cusolverStatus_t sucess =
              cusolverDnSgeqrf_bufferSize(handle,m,n,A,lda,Lwork);
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t
  cusolver_geqrf_bufferSize(cusolverDnHandle_t handle,
                      int m,
                      int n,
                      double *A,
                      int lda,
                      int *Lwork )
  {
    cusolverStatus_t sucess =
              cusolverDnDgeqrf_bufferSize(handle,m,n,A,lda,Lwork);
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t
  cusolver_geqrf_bufferSize(cusolverDnHandle_t handle,
                      int m,
                      int n,
                      std::complex<float> *A,
                      int lda,
                      int *Lwork )
  {
    cusolverStatus_t sucess =
              cusolverDnCgeqrf_bufferSize(handle,m,n,reinterpret_cast<cuComplex *>(A),lda,Lwork);
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t
  cusolver_geqrf_bufferSize(cusolverDnHandle_t handle,
                      int m,
                      int n,
                      std::complex<double> *A,
                      int lda,
                      int *Lwork )
  {
    cusolverStatus_t sucess =
              cusolverDnZgeqrf_bufferSize(handle,m,n,reinterpret_cast<cuDoubleComplex *>(A),lda,Lwork);
    cudaDeviceSynchronize ();
    return sucess;
  }

  //geqrf
  inline cusolverStatus_t 
  cusolver_geqrf(cusolverDnHandle_t handle,
           int m,
           int n,
           float *A,
           int lda,
           float *TAU,
           float *Workspace,
           int Lwork,
           int *devInfo )
  {
    cusolverStatus_t sucess =
              cusolverDnSgeqrf(handle,m,n,A,lda,TAU,Workspace,Lwork,devInfo);
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t
  cusolver_geqrf(cusolverDnHandle_t handle,
           int m,
           int n,
           double *A,
           int lda,
           double *TAU,
           double *Workspace,
           int Lwork,
           int *devInfo )
  {
    cusolverStatus_t sucess =
              cusolverDnDgeqrf(handle,m,n,A,lda,TAU,Workspace,Lwork,devInfo);
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t
  cusolver_geqrf(cusolverDnHandle_t handle,
           int m,
           int n,
           std::complex<float> *A,
           int lda,
           std::complex<float> *TAU,
           std::complex<float> *Workspace,
           int Lwork,
           int *devInfo )
  {
    cusolverStatus_t sucess =
              cusolverDnCgeqrf(handle,m,n,reinterpret_cast<cuComplex *>(A),lda,
                               reinterpret_cast<cuComplex *>(TAU),reinterpret_cast<cuComplex *>(Workspace),Lwork,devInfo);
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t
  cusolver_geqrf(cusolverDnHandle_t handle,
           int m,
           int n,
           std::complex<double> *A,
           int lda,
           std::complex<double> *TAU,
           std::complex<double> *Workspace,
           int Lwork,
           int *devInfo )
  {
    cusolverStatus_t sucess =
              cusolverDnZgeqrf(handle,m,n,reinterpret_cast<cuDoubleComplex *>(A),lda,
                               reinterpret_cast<cuDoubleComplex *>(TAU),reinterpret_cast<cuDoubleComplex *>(Workspace),Lwork,devInfo);
    cudaDeviceSynchronize ();
    return sucess;
  }



  //gqr_bufferSize
  inline cusolverStatus_t 
  cusolver_gqr_bufferSize(
    cusolverDnHandle_t handle,
    int m,
    int n,
    int k,
    const float *A,
    int lda,
    const float *tau,
    int *lwork)
  {
    cusolverStatus_t sucess =
              cusolverDnSorgqr_bufferSize(handle,m,n,k,A,lda,tau,lwork);
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t
  cusolver_gqr_bufferSize(
    cusolverDnHandle_t handle,
    int m,
    int n,
    int k,
    const double *A,
    int lda,
    const double *tau,
    int *lwork)
  {
    cusolverStatus_t sucess =
              cusolverDnDorgqr_bufferSize(handle,m,n,k,A,lda,tau,lwork);
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t
  cusolver_gqr_bufferSize(
    cusolverDnHandle_t handle,
    int m,
    int n,
    int k,
    const std::complex<float> *A,
    int lda,
    const std::complex<float> *tau,
    int *lwork)
  {
    cusolverStatus_t sucess =
              cusolverDnCungqr_bufferSize(handle,m,n,k,reinterpret_cast<cuComplex const *>(A),lda,  
                                          reinterpret_cast<cuComplex const *>(tau),lwork);
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t
  cusolver_gqr_bufferSize(
    cusolverDnHandle_t handle,
    int m,
    int n,
    int k,
    const std::complex<double> *A,
    int lda,
    const std::complex<double> *tau,
    int *lwork)
  {
    cusolverStatus_t sucess =
              cusolverDnZungqr_bufferSize(handle,m,n,k,reinterpret_cast<cuDoubleComplex const*>(A),lda,
                                          reinterpret_cast<cuDoubleComplex const *>(tau),lwork);
    cudaDeviceSynchronize ();
    return sucess;
  }

  //gqr
  inline cusolverStatus_t 
  cusolver_gqr(
    cusolverDnHandle_t handle,
    int m,
    int n,
    int k,
    float *A,
    int lda,
    const float *tau,
    float *work,
    int lwork,
    int *devInfo)
  {
    cusolverStatus_t sucess =
              cusolverDnSorgqr(handle,m,n,k,A,lda,tau,work,lwork,devInfo);
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t
  cusolver_gqr(
    cusolverDnHandle_t handle,
    int m,
    int n,
    int k,
    double *A,
    int lda,
    const double *tau,
    double *work,
    int lwork,
    int *devInfo)
  {
    cusolverStatus_t sucess =
              cusolverDnDorgqr(handle,m,n,k,A,lda,tau,work,lwork,devInfo);
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t
  cusolver_gqr(
    cusolverDnHandle_t handle,
    int m,
    int n,
    int k,
    std::complex<float> *A,
    int lda,
    const std::complex<float> *tau,
    std::complex<float> *work,
    int lwork,
    int *devInfo)
  {
    cusolverStatus_t sucess =
              cusolverDnCungqr(handle,m,n,k,reinterpret_cast<cuComplex *>(A),lda,
                               reinterpret_cast<cuComplex const *>(tau),reinterpret_cast<cuComplex *>(work),lwork,devInfo);
    cudaDeviceSynchronize ();
    return sucess;
  }

  inline cusolverStatus_t
  cusolver_gqr(
    cusolverDnHandle_t handle,
    int m,
    int n,
    int k,
    std::complex<double> *A,
    int lda,
    const std::complex<double> *tau,
    std::complex<double> *work,
    int lwork,
    int *devInfo)
  {
    cusolverStatus_t sucess =
              cusolverDnZungqr(handle,m,n,k,reinterpret_cast<cuDoubleComplex *>(A),lda,
                               reinterpret_cast<cuDoubleComplex const *>(tau),reinterpret_cast<cuDoubleComplex *>(work),lwork,devInfo);
    cudaDeviceSynchronize ();
    return sucess;
  }



}

#endif
