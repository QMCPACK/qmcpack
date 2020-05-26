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

#ifndef HIPBLAS_FUNCTIONDEFS_H
#define HIPBLAS_FUNCTIONDEFS_H

#include<cassert>
#include <hip_runtime.h>
#include "hipblas.h"
#include "AFQMC/Memory/HIP/hip_utilities.h"

namespace hipblas {

  using qmc_hip::hipblasOperation;

  // Level-1
  inline hipblasStatus_t hipblas_copy(hipblasHandle_t handle, int n,
                           float *x, int incx,
                           float *y, int incy)
  {
    hipblasStatus_t sucess =
                hipblasScopy(handle,n,x,incx,y,incy);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_copy(hipblasHandle_t handle, int n,
                           double *x, int incx,
                           double *y, int incy)
  {
    hipblasStatus_t sucess =
                hipblasDcopy(handle,n,x,incx,y,incy);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_copy(hipblasHandle_t handle, int n,
                           std::complex<float> *x, int incx,
                           std::complex<float> *y, int incy)
  {
    hipblasStatus_t sucess =
                hipblasCcopy(handle,n,
                        reinterpret_cast<hipComplex *>(x),incx,
                        reinterpret_cast<hipComplex *>(y),incy);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_copy(hipblasHandle_t handle, int n,
                           std::complex<double> *x, int incx,
                           std::complex<double> *y, int incy)
  {
    hipblasStatus_t sucess =
                hipblasZcopy(handle,n,
                        reinterpret_cast<hipDoubleComplex *>(x),incx,
                        reinterpret_cast<hipDoubleComplex *>(y),incy);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_scal(hipblasHandle_t handle, int n,
                           const float alpha, float *x, int incx)
  {
    hipblasStatus_t sucess =
                hipblasSscal(handle,n,&alpha,x,incx);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_scal(hipblasHandle_t handle, int n,
                           const double alpha, double *x, int incx)
  {
    hipblasStatus_t sucess =
                hipblasDscal(handle,n,&alpha,x,incx);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_scal(hipblasHandle_t handle, int n,
                           const std::complex<float> alpha, std::complex<float> *x, int incx)
  {
    hipblasStatus_t sucess =
                hipblasCscal(handle,n,
                                reinterpret_cast<hipComplex const*>(&alpha),
                                reinterpret_cast<hipComplex *>(x),incx);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_scal(hipblasHandle_t handle, int n,
                           const std::complex<double> alpha, std::complex<double> *x, int incx)
  {
    hipblasStatus_t sucess =
                hipblasZscal(handle,n,
                                reinterpret_cast<hipDoubleComplex const*>(&alpha),
                                reinterpret_cast<hipDoubleComplex *>(x),incx);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline float hipblas_dot(hipblasHandle_t handle, int n,
                           const float *x, int incx,
                           const float *y, int incy)
  {
    float result;
    hipblasStatus_t sucess = hipblasSdot(handle,n,x,incx,y,incy,&result);
    hipDeviceSynchronize ();
    if(HIPBLAS_STATUS_SUCCESS != sucess)
      throw std::runtime_error("Error: hipblas_dot returned error code.");
    return result;
  }

  inline double hipblas_dot(hipblasHandle_t handle, int n,
                           const double *x, int incx,
                           const double *y, int incy)
  {
    double result;
    hipblasStatus_t sucess = hipblasDdot(handle,n,x,incx,y,incy,&result);
    hipDeviceSynchronize ();
    if(HIPBLAS_STATUS_SUCCESS != sucess)
      throw std::runtime_error("Error: hipblas_dot returned error code.");
    return result;
  }

  inline std::complex<float> hipblas_dot(hipblasHandle_t handle, int n,
                           const std::complex<float> *x, int incx,
                           const std::complex<float> *y, int incy)
  {
    std::complex<float> result;
    hipblasStatus_t sucess = hipblasCdotu(handle,n,
                                reinterpret_cast<hipComplex const*>(x),incx,
                                reinterpret_cast<hipComplex const*>(y),incy,
                                reinterpret_cast<hipComplex *>(&result));
    hipDeviceSynchronize ();
    if(HIPBLAS_STATUS_SUCCESS != sucess)
      throw std::runtime_error("Error: hipblas_dot returned error code.");
    return result;
  }

  inline std::complex<double> hipblas_dot(hipblasHandle_t handle, int n,
                           const std::complex<double> *x, int incx,
                           const std::complex<double> *y, int incy)
  {
    std::complex<double> result;
    hipblasStatus_t sucess = hipblasZdotu(handle,n,
                                reinterpret_cast<hipDoubleComplex const*>(x),incx,
                                reinterpret_cast<hipDoubleComplex const*>(y),incy,
                                reinterpret_cast<hipDoubleComplex *>(&result));
    hipDeviceSynchronize ();
    if(HIPBLAS_STATUS_SUCCESS != sucess)
      throw std::runtime_error("Error: hipblas_dot returned error code.");
    return result;
  }

  inline std::complex<double> hipblas_dot(hipblasHandle_t handle, int n,
                           const double *x, int incx,
                           const std::complex<double> *y, int incy)
  {
    int incy_ = 2*incy;
    const double* y_ = reinterpret_cast<const double*>(y);
    const double* y1_ = y_+1;
    double resR, resI;
    hipblasStatus_t sucess = hipblasDdot(handle,n,x,incx,y_,incy_,&resR);
    hipDeviceSynchronize ();
    if(HIPBLAS_STATUS_SUCCESS != sucess)
      throw std::runtime_error("Error: hipblas_dot returned error code.");
    sucess = hipblasDdot(handle,n,x,incx,y1_,incy_,&resI);
    hipDeviceSynchronize ();
    if(HIPBLAS_STATUS_SUCCESS != sucess)
      throw std::runtime_error("Error: hipblas_dot returned error code.");
    return std::complex<double>{resR,resI};
  }

  inline std::complex<double> hipblas_dot(hipblasHandle_t handle, int n,
                           const std::complex<double> *x, int incx,
                           const double *y, int incy)
  {
    int incx_ = 2*incx;
    const double* x_ = reinterpret_cast<const double*>(x);
    const double* x1_ = x_+1;
    double resR, resI;
    hipblasStatus_t sucess = hipblasDdot(handle,n,x_,incx_,y,incy,&resR);
    hipDeviceSynchronize ();
    if(HIPBLAS_STATUS_SUCCESS != sucess)
      throw std::runtime_error("Error: hipblas_dot returned error code.");
    sucess = hipblasDdot(handle,n,x1_,incx_,y,incy,&resI);
    hipDeviceSynchronize ();
    if(HIPBLAS_STATUS_SUCCESS != sucess)
      throw std::runtime_error("Error: hipblas_dot returned error code.");
    return std::complex<double>{resR,resI};
  }

  inline hipblasStatus_t hipblas_axpy(hipblasHandle_t handle, int n,
                           const float alpha, const float *x, int incx,
                           float *y, int incy)
  {
    hipblasStatus_t sucess =
                hipblasSaxpy(handle,n,&alpha,x,incx,y,incy);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_axpy(hipblasHandle_t handle, int n,
                           const double alpha, const double *x, int incx,
                           double *y, int incy)
  {
    hipblasStatus_t sucess =
                hipblasDaxpy(handle,n,&alpha,x,incx,y,incy);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_axpy(hipblasHandle_t handle, int n,
                           const std::complex<float> alpha, const std::complex<float> *x, int incx,
                           std::complex<float> *y, int incy)
  {
    hipblasStatus_t sucess =
                hipblasCaxpy(handle,n,
                                reinterpret_cast<hipComplex const*>(&alpha),
                                reinterpret_cast<hipComplex const*>(x),incx,
                                reinterpret_cast<hipComplex *>(y),incy);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_axpy(hipblasHandle_t handle, int n,
                           const std::complex<double> alpha, const std::complex<double> *x, int incx,
                           std::complex<double> *y, int incy)
  {
    hipblasStatus_t sucess =
                hipblasZaxpy(handle,n,
                                reinterpret_cast<hipDoubleComplex const*>(&alpha),
                                reinterpret_cast<hipDoubleComplex const*>(x),incx,
                                reinterpret_cast<hipDoubleComplex *>(y),incy);
    hipDeviceSynchronize ();
    return sucess;
  }

  // Level-2
  inline hipblasStatus_t hipblas_gemv(hipblasHandle_t handle,
                          char Atrans, int M, int N,
                          const float alpha,
                          const float * A, int lda,
                          const float * x, int incx,
                          const float beta,
                          float * y, int incy)
  {
    hipblasStatus_t sucess =
                hipblasSgemv(handle,hipblasOperation(Atrans),
                           M,N,&alpha,A,lda,x,incx,&beta,y,incy);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_gemv(hipblasHandle_t handle,
                          char Atrans, int M, int N,
                          const double alpha,
                          const double * A, int lda,
                          const double * x, int incx,
                          const double beta,
                          double * y, int incy)
  {
    hipblasStatus_t sucess =
                hipblasDgemv(handle,hipblasOperation(Atrans),
                           M,N,&alpha,A,lda,x,incx,&beta,y,incy);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_gemv(hipblasHandle_t handle,
                          char Atrans, int M, int N,
                          const std::complex<float> alpha,
                          const std::complex<float> * A, int lda,
                          const std::complex<float> * x, int incx,
                          const std::complex<float> beta,
                          std::complex<float> *y, int incy)
  {
    hipblasStatus_t sucess =
                hipblasCgemv(handle,hipblasOperation(Atrans),M,N,
                           reinterpret_cast<hipComplex const*>(&alpha),reinterpret_cast<hipComplex const*>(A),lda,
                           reinterpret_cast<hipComplex const*>(x),incx,reinterpret_cast<hipComplex const*>(&beta),
                           reinterpret_cast<hipComplex *>(y),incy);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_gemv(hipblasHandle_t handle,
                          char Atrans, int M, int N,
                          const std::complex<double> alpha,
                          const std::complex<double> * A, int lda,
                          const std::complex<double> * x, int incx,
                          const std::complex<double> beta,
                          std::complex<double> * y, int incy)
  {
    hipblasStatus_t sucess =
                hipblasZgemv(handle,hipblasOperation(Atrans),M,N,
                           reinterpret_cast<hipDoubleComplex const*>(&alpha),
                           reinterpret_cast<hipDoubleComplex const*>(A),lda,
                           reinterpret_cast<hipDoubleComplex const*>(x),incx,
                           reinterpret_cast<hipDoubleComplex const*>(&beta),
                           reinterpret_cast<hipDoubleComplex *>(y),incy);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_gemv(hipblasHandle_t handle,
                          char Atrans, int M, int N,
                          const float alpha,
                          const float * A, int lda,
                          const std::complex<float> * x, int incx,
                          const float beta,
                          std::complex<float> *y, int incy)
  {
    hipblasStatus_t sucess;
    char Nt('N');
    char Tt('T');
    if(Atrans == 'n' || Atrans == 'N')
      sucess = hipblasSgemm(handle,hipblasOperation(Nt),hipblasOperation(Tt),2,M,N,&alpha,
                           reinterpret_cast<float const*>(x),2*incx,
                           A,lda,&beta,
                           reinterpret_cast<float *>(y),2*incy);
    else if(Atrans == 't' || Atrans == 'T')
      sucess = hipblasSgemm(handle,hipblasOperation(Nt),hipblasOperation(Nt),2,N,M,&alpha,
                           reinterpret_cast<float const*>(x),2*incx,
                           A,lda,&beta,
                           reinterpret_cast<float *>(y),2*incy);
    else
      assert(0);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_gemv(hipblasHandle_t handle,
                          char Atrans, int M, int N,
                          const double alpha,
                          const double * A, int lda,
                          const std::complex<double> * x, int incx,
                          const double beta,
                          std::complex<double> * y, int incy)
  {
    hipblasStatus_t sucess;
    char Nt('N');
    char Tt('T');
    if(Atrans == 'n' || Atrans == 'N')
      sucess = hipblasDgemm(handle,hipblasOperation(Nt),hipblasOperation(Tt),2,M,N,&alpha,
                           reinterpret_cast<double const*>(x),2*incx,
                           A,lda,&beta,
                           reinterpret_cast<double *>(y),2*incy);
    else if(Atrans == 't' || Atrans == 'T')
      sucess = hipblasDgemm(handle,hipblasOperation(Nt),hipblasOperation(Nt),2,N,M,&alpha,
                           reinterpret_cast<double const*>(x),2*incx,
                           A,lda,&beta,
                           reinterpret_cast<double *>(y),2*incy);
    else
      assert(0);
    hipDeviceSynchronize ();
    return sucess;
  }


  // Level-3
  inline hipblasStatus_t hipblas_gemm(hipblasHandle_t handle,
                          char Atrans, char Btrans, int M, int N, int K,
                          const float alpha,
                          const float * A, int lda,
                          const float * B, int ldb,
                          const float beta,
                          float * C, int ldc)
  {
    hipblasStatus_t sucess =
                hipblasSgemm(handle,
                           hipblasOperation(Atrans),hipblasOperation(Btrans),
                           M,N,K,&alpha,A,lda,B,ldb,&beta,C,ldc);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_gemm(hipblasHandle_t handle,
                          char Atrans, char Btrans, int M, int N, int K,
                          const double alpha,
                          const double * A, int lda,
                          const double * B, int ldb,
                          const double beta,
                          double * C, int ldc)
  {
    hipblasStatus_t sucess =
                hipblasDgemm(handle,
                           hipblasOperation(Atrans),hipblasOperation(Btrans),
                           M,N,K,&alpha,A,lda,B,ldb,&beta,C,ldc);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_gemm(hipblasHandle_t handle,
                          char Atrans, char Btrans, int M, int N, int K,
                          const std::complex<float> alpha,
                          const std::complex<float> * A, int lda,
                          const std::complex<float> * B, int ldb,
                          const std::complex<float> beta,
                          std::complex<float> * C, int ldc)
  {
    hipblasStatus_t sucess =
                hipblasCgemm(handle,
                           hipblasOperation(Atrans),hipblasOperation(Btrans),M,N,K,
                           reinterpret_cast<hipComplex const*>(&alpha),
                           reinterpret_cast<hipComplex const*>(A),lda,
                           reinterpret_cast<hipComplex const*>(B),ldb,
                           reinterpret_cast<hipComplex const*>(&beta),
                           reinterpret_cast<hipComplex *>(C),ldc);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_gemm(hipblasHandle_t handle,
                          char Atrans, char Btrans, int M, int N, int K,
                          const std::complex<double> alpha,
                          const std::complex<double> * A, int lda,
                          const std::complex<double> * B, int ldb,
                          const std::complex<double> beta,
                          std::complex<double> * C, int ldc)
  {
    hipblasStatus_t sucess =
                hipblasZgemm(handle,
                           hipblasOperation(Atrans),hipblasOperation(Btrans),M,N,K,
                           reinterpret_cast<hipDoubleComplex const*>(&alpha),
                           reinterpret_cast<hipDoubleComplex const*>(A),lda,
                           reinterpret_cast<hipDoubleComplex const*>(B),ldb,
                           reinterpret_cast<hipDoubleComplex const*>(&beta),
                           reinterpret_cast<hipDoubleComplex *>(C),ldc);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_gemm(hipblasHandle_t handle,
                          char Atrans, char Btrans, int M, int N, int K,
                          const float alpha,
                          const std::complex<float> * A, int lda,
                          const float * B, int ldb,
                          const float beta,
                          std::complex<float> * C, int ldc)
  {
    assert(Atrans=='n' || Atrans=='N');
    hipblasStatus_t sucess =
                hipblasSgemm(handle,
                           hipblasOperation(Atrans),hipblasOperation(Btrans),2*M,N,K,&alpha,
                           reinterpret_cast<float const*>(A),2*lda,
                           B,ldb,&beta,
                           reinterpret_cast<float *>(C),2*ldc);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_gemm(hipblasHandle_t handle,
                          char Atrans, char Btrans, int M, int N, int K,
                          double alpha,
                          const std::complex<double> * A, int lda,
                          const double * B, int ldb,
                          double beta,
                          std::complex<double> * C, int ldc)
  {
    assert(Atrans=='n' || Atrans=='N');
    hipblasStatus_t sucess =
                hipblasDgemm(handle,
                           hipblasOperation(Atrans),hipblasOperation(Btrans),2*M,N,K,&alpha,
                           reinterpret_cast<double const*>(A),2*lda,
                           B,ldb,&beta,
                           reinterpret_cast<double *>(C),2*ldc);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_gemm(hipblasHandle_t handle,
                          char Atrans, char Btrans, int M, int N, int K,
                          const hipDoubleComplex alpha,
                          const hipDoubleComplex * A, int lda,
                          const hipDoubleComplex * B, int ldb,
                          const hipDoubleComplex beta,
                          hipDoubleComplex * C, int ldc)
  {
    hipblasStatus_t sucess =
                hipblasZgemm(handle,
                           hipblasOperation(Atrans),hipblasOperation(Btrans),
                           M,N,K,&alpha,A,lda,B,ldb,&beta,C,ldc);
    hipDeviceSynchronize ();
    return sucess;
  }

  // Extensions
  inline hipblasStatus_t hipblas_getrfBatched(hipblasHandle_t handle,
                                   int n,
                                   float **Aarray,
                                   int lda,
                                   int *PivotArray,
                                   int *infoArray,
                                   int batchSize)
  {
    hipblasStatus_t sucess =
              hipblasSgetrfBatched(handle,n,Aarray,lda,PivotArray,infoArray,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_getrfBatched(hipblasHandle_t handle,
                                   int n,
                                   double **Aarray,
                                   int lda,
                                   int *PivotArray,
                                   int ldc,
                                   int *infoArray,
                                   int batchSize)
  {
    hipblasStatus_t sucess =
              hipblasDgetrfBatched(handle,n,Aarray,lda,PivotArray,infoArray,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_getrfBatched(hipblasHandle_t handle,
                                   int n,
                                   std::complex<double> **Aarray,
                                   int lda,
                                   int *PivotArray,
                                   int *infoArray,
                                   int batchSize)
  {
    hipblasStatus_t sucess =
              hipblasZgetrfBatched(handle,n,
                            reinterpret_cast<hipDoubleComplex *const *>(Aarray),lda,PivotArray,
                            infoArray,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_getrfBatched(hipblasHandle_t handle,
                                   int n,
                                   std::complex<float> **Aarray,
                                   int lda,
                                   int *PivotArray,
                                   int *infoArray,
                                   int batchSize)
  {
    hipblasStatus_t sucess =
              hipblasCgetrfBatched(handle,n,
                            reinterpret_cast<hipComplex *const *>(Aarray),lda,PivotArray,
                            infoArray,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_getriBatched(hipblasHandle_t handle,
                                   int n,
                                   float **Aarray,
                                   int lda,
                                   int *PivotArray,
                                   float **Carray,
                                   int ldc,
                                   int *infoArray,
                                   int batchSize)
  {
    hipblasStatus_t sucess =
              hipblasSgetriBatched(handle,n,Aarray,lda,PivotArray,Carray,ldc,infoArray,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_getriBatched(hipblasHandle_t handle,
                                   int n,
                                   double **Aarray,
                                   int lda,
                                   int *PivotArray,
                                   double **Carray,
                                   int ldc,
                                   int *infoArray,
                                   int batchSize)
  {
    hipblasStatus_t sucess =
              hipblasDgetriBatched(handle,n,Aarray,lda,PivotArray,Carray,ldc,infoArray,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_getriBatched(hipblasHandle_t handle,
                                   int n,
                                   std::complex<double> **Aarray,
                                   int lda,
                                   int *PivotArray,
                                   std::complex<double> **Carray,
                                   int ldc,
                                   int *infoArray,
                                   int batchSize)
  {
    hipblasStatus_t sucess =
              hipblasZgetriBatched(handle,n,
                            reinterpret_cast<const hipDoubleComplex *const *>(Aarray),lda,PivotArray,
                            reinterpret_cast<hipDoubleComplex *const *>(Carray),ldc,infoArray,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_getriBatched(hipblasHandle_t handle,
                                   int n,
                                   std::complex<float> **Aarray,
                                   int lda,
                                   int *PivotArray,
                                   std::complex<float> **Carray,
                                   int ldc,
                                   int *infoArray,
                                   int batchSize)
  {
    hipblasStatus_t sucess =
              hipblasCgetriBatched(handle,n,
                            reinterpret_cast<const hipComplex *const *>(Aarray),lda,PivotArray,
                            reinterpret_cast<hipComplex *const *>(Carray),ldc,infoArray,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_matinvBatched(hipblasHandle_t handle,
                                   int n,
                                   float **Aarray,
                                   int lda,
                                   float **Carray,
                                   int ldc,
                                   int *infoArray,
                                   int batchSize)
  {
    hipblasStatus_t sucess =
              hipblasSmatinvBatched(handle,n,Aarray,lda,Carray,ldc,infoArray,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_matinvBatched(hipblasHandle_t handle,
                                   int n,
                                   double **Aarray,
                                   int lda,
                                   double **Carray,
                                   int ldc,
                                   int *infoArray,
                                   int batchSize)
  {
    hipblasStatus_t sucess =
              hipblasDmatinvBatched(handle,n,Aarray,lda,Carray,ldc,infoArray,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_matinvBatched(hipblasHandle_t handle,
                                   int n,
                                   std::complex<float> **Aarray,
                                   int lda,
                                   std::complex<float> **Carray,
                                   int ldc,
                                   int *infoArray,
                                   int batchSize)
  {
    hipblasStatus_t sucess =
              hipblasCmatinvBatched(handle,n,
                        reinterpret_cast<const hipComplex * const*>(Aarray),lda,
                        reinterpret_cast<hipComplex **>(Carray),ldc,
                        infoArray,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_matinvBatched(hipblasHandle_t handle,
                                   int n,
                                   std::complex<double> **Aarray,
                                   int lda,
                                   std::complex<double> **Carray,
                                   int ldc,
                                   int *infoArray,
                                   int batchSize)
  {
    hipblasStatus_t sucess =
              hipblasZmatinvBatched(handle,n,
                        reinterpret_cast<const hipDoubleComplex *const*>(Aarray),lda,
                        reinterpret_cast<hipDoubleComplex **>(Carray),ldc,
                        infoArray,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_geam(hipblasHandle_t handle, char Atrans, char Btrans, int M, int N,
                         float const alpha,
                         float const* A, int lda,
                         float const beta,
                         float const* B, int ldb,
                         float *C, int ldc)
  {
    hipblasStatus_t sucess =
                hipblasSgeam(handle,hipblasOperation(Atrans),hipblasOperation(Btrans),
                           M,N,&alpha,A,lda,&beta,B,ldb,C,ldc);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_geam(hipblasHandle_t handle, char Atrans, char Btrans, int M, int N,
                         double const alpha,
                         double const* A, int lda,
                         double const beta,
                         double const* B, int ldb,
                         double *C, int ldc)
  {
    hipblasStatus_t sucess =
                hipblasDgeam(handle,hipblasOperation(Atrans),hipblasOperation(Btrans),
                           M,N,&alpha,A,lda,&beta,B,ldb,C,ldc);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_geam(hipblasHandle_t handle, char Atrans, char Btrans, int M, int N,
                         std::complex<float> const alpha,
                         std::complex<float> const* A, int lda,
                         std::complex<float> const beta,
                         std::complex<float> const* B, int ldb,
                         std::complex<float> *C, int ldc)
  {
    hipblasStatus_t sucess =
                hipblasCgeam(handle,hipblasOperation(Atrans),hipblasOperation(Btrans),M,N,
                           reinterpret_cast<hipComplex const*>(&alpha),
                           reinterpret_cast<hipComplex const*>(A),lda,
                           reinterpret_cast<hipComplex const*>(&beta),
                           reinterpret_cast<hipComplex const*>(B),ldb,
                           reinterpret_cast<hipComplex *>(C),ldc);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_geam(hipblasHandle_t handle, char Atrans, char Btrans, int M, int N,
                         std::complex<double> const alpha,
                         std::complex<double> const* A, int lda,
                         std::complex<double> const beta,
                         std::complex<double> const* B, int ldb,
                         std::complex<double> *C, int ldc)
  {
    hipblasStatus_t sucess =
                hipblasZgeam(handle,hipblasOperation(Atrans),hipblasOperation(Btrans),M,N,
                           reinterpret_cast<hipDoubleComplex const*>(&alpha),
                           reinterpret_cast<hipDoubleComplex const*>(A),lda,
                           reinterpret_cast<hipDoubleComplex const*>(&beta),
                           reinterpret_cast<hipDoubleComplex const*>(B),ldb,
                           reinterpret_cast<hipDoubleComplex *>(C),ldc);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_gemmStridedBatched(hipblasHandle_t handle,
                          char Atrans, char Btrans, int M, int N, int K,
                          const float alpha,
                          const float * A, int lda, int strideA,
                          const float * B, int ldb, int strideB,
                          const float beta,
                          float * C, int ldc, int strideC, int batchSize)
  {
    hipblasStatus_t sucess =
                hipblasSgemmStridedBatched(handle,
                           hipblasOperation(Atrans),hipblasOperation(Btrans),
                           M,N,K,&alpha,A,lda,strideA,B,ldb,strideB,&beta,C,ldc,strideC,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_gemmStridedBatched(hipblasHandle_t handle,
                          char Atrans, char Btrans, int M, int N, int K,
                          const double alpha,
                          const double * A, int lda, int strideA,
                          const double * B, int ldb, int strideB,
                          const double beta,
                          double * C, int ldc, int strideC, int batchSize)
  {
    hipblasStatus_t sucess =
                hipblasDgemmStridedBatched(handle,
                           hipblasOperation(Atrans),hipblasOperation(Btrans),
                           M,N,K,&alpha,A,lda,strideA,B,ldb,strideB,&beta,C,ldc,strideC,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_gemmStridedBatched(hipblasHandle_t handle,
                          char Atrans, char Btrans, int M, int N, int K,
                          const std::complex<float> alpha,
                          const std::complex<float> * A, int lda, int strideA,
                          const std::complex<float> * B, int ldb, int strideB,
                          const std::complex<float> beta,
                          std::complex<float> * C, int ldc, int strideC, int batchSize)
  {
    hipblasStatus_t sucess =
                hipblasCgemmStridedBatched(handle,
                           hipblasOperation(Atrans),hipblasOperation(Btrans),
                           M,N,K,
                           reinterpret_cast<hipComplex const*>(&alpha),
                           reinterpret_cast<hipComplex const*>(A),lda,strideA,
                           reinterpret_cast<hipComplex const*>(B),ldb,strideB,
                           reinterpret_cast<hipComplex const*>(&beta),
                           reinterpret_cast<hipComplex *>(C),ldc,strideC,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_gemmStridedBatched(hipblasHandle_t handle,
                          char Atrans, char Btrans, int M, int N, int K,
                          const std::complex<double> alpha,
                          const std::complex<double> * A, int lda, int strideA,
                          const std::complex<double> * B, int ldb, int strideB,
                          const std::complex<double> beta,
                          std::complex<double> * C, int ldc, int strideC, int batchSize)
  {
    hipblasStatus_t sucess =
                hipblasZgemmStridedBatched(handle,
                           hipblasOperation(Atrans),hipblasOperation(Btrans),
                           M,N,K,
                           reinterpret_cast<hipDoubleComplex const*>(&alpha),
                           reinterpret_cast<hipDoubleComplex const*>(A),lda,strideA,
                           reinterpret_cast<hipDoubleComplex const*>(B),ldb,strideB,
                           reinterpret_cast<hipDoubleComplex const*>(&beta),
                           reinterpret_cast<hipDoubleComplex *>(C),ldc,strideC,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_gemmBatched(hipblasHandle_t handle,
                          char Atrans, char Btrans, int M, int N, int K,
                          float alpha,
                          float ** A, int lda,
                          float ** B, int ldb,
                          float beta,
                          float ** C, int ldc, int batchSize)
  {
    hipblasStatus_t sucess =
                hipblasSgemmBatched(handle,
                           hipblasOperation(Atrans),hipblasOperation(Btrans),
                           M,N,K,&alpha,A,lda,B,ldb,&beta,C,ldc,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_gemmBatched(hipblasHandle_t handle,
                          char Atrans, char Btrans, int M, int N, int K,
                          double alpha,
                          double ** A, int lda,
                          double ** B, int ldb,
                          double beta,
                          double ** C, int ldc, int batchSize)
  {
    hipblasStatus_t sucess =
                hipblasDgemmBatched(handle,
                           hipblasOperation(Atrans),hipblasOperation(Btrans),
                           M,N,K,&alpha,A,lda,B,ldb,&beta,C,ldc,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_gemmBatched(hipblasHandle_t handle,
                          char Atrans, char Btrans, int M, int N, int K,
                          std::complex<float> alpha,
                          std::complex<float> ** A, int lda,
                          std::complex<float> ** B, int ldb,
                          std::complex<float> beta,
                          std::complex<float> ** C, int ldc, int batchSize)
  {
    hipblasStatus_t sucess =
                hipblasCgemmBatched(handle,
                        hipblasOperation(Atrans),hipblasOperation(Btrans),M,N,K,
                        reinterpret_cast<hipComplex *>(&alpha),
                        reinterpret_cast<hipComplex **>(A),lda,
                        reinterpret_cast<hipComplex **>(B),ldb,
                        reinterpret_cast<hipComplex *>(&beta),
                        reinterpret_cast<hipComplex **>(C),ldc,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_gemmBatched(hipblasHandle_t handle,
                          char Atrans, char Btrans, int M, int N, int K,
                          std::complex<double> alpha,
                          std::complex<double> ** A, int lda,
                          std::complex<double> ** B, int ldb,
                          std::complex<double> beta,
                          std::complex<double> ** C, int ldc, int batchSize)
  {
    hipblasStatus_t sucess =
                hipblasZgemmBatched(handle,
                        hipblasOperation(Atrans),hipblasOperation(Btrans),M,N,K,
                        reinterpret_cast<hipDoubleComplex *>(&alpha),
                        reinterpret_cast<hipDoubleComplex **>(A),lda,
                        reinterpret_cast<hipDoubleComplex **>(B),ldb,
                        reinterpret_cast<hipDoubleComplex *>(&beta),
                        reinterpret_cast<hipDoubleComplex **>(C),ldc,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_gemmBatched(hipblasHandle_t handle,
                          char Atrans, char Btrans, int M, int N, int K,
                          float alpha,
                          std::complex<float> ** A, int lda,
                          float ** B, int ldb,
                          float beta,
                          std::complex<float> ** C, int ldc, int batchSize)
  {
    hipblasStatus_t sucess =
                hipblasSgemmBatched(handle,
                        hipblasOperation(Atrans),hipblasOperation(Btrans),2*M,N,K,
                        &alpha,
                        reinterpret_cast<float **>(A),2*lda,
                        B, ldb, &beta,
                        reinterpret_cast<float **>(C),2*ldc,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_gemmBatched(hipblasHandle_t handle,
                          char Atrans, char Btrans, int M, int N, int K,
                          double alpha,
                          std::complex<double> ** A, int lda,
                          double ** B, int ldb,
                          double beta,
                          std::complex<double> ** C, int ldc, int batchSize)
  {
    hipblasStatus_t sucess =
                hipblasDgemmBatched(handle,
                        hipblasOperation(Atrans),hipblasOperation(Btrans),2*M,N,K,
                        &alpha,
                        reinterpret_cast<double **>(A),2*lda,
                        B, ldb, &beta,
                        reinterpret_cast<double **>(C),2*ldc,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_geqrfBatched( hipblasHandle_t handle,
                                    int m,
                                    int n,
                                    double **Aarray,
                                    int lda,
                                    double **TauArray,
                                    int *info,
                                    int batchSize)
  {
    hipblasStatus_t sucess = hipblasDgeqrfBatched(handle,m,n,Aarray,lda,TauArray,info,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_geqrfBatched( hipblasHandle_t handle,
                                    int m,
                                    int n,
                                    float **Aarray,
                                    int lda,
                                    float **TauArray,
                                    int *info,
                                    int batchSize)
  {
    hipblasStatus_t sucess = hipblasSgeqrfBatched(handle,m,n,Aarray,lda,TauArray,info,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }


  inline hipblasStatus_t hipblas_geqrfBatched( hipblasHandle_t handle,
                                    int m,
                                    int n,
                                    std::complex<double> **Aarray,
                                    int lda,
                                    std::complex<double> **TauArray,
                                    int *info,
                                    int batchSize)
  {
     hipblasStatus_t sucess = hipblasZgeqrfBatched(handle,m,n,
                        reinterpret_cast<hipDoubleComplex **>(Aarray),lda,
                        reinterpret_cast<hipDoubleComplex **>(TauArray),info,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline hipblasStatus_t hipblas_geqrfBatched( hipblasHandle_t handle,
                                    int m,
                                    int n,
                                    std::complex<float> **Aarray,
                                    int lda,
                                    std::complex<float> **TauArray,
                                    int *info,
                                    int batchSize)
  {
     hipblasStatus_t sucess = hipblasCgeqrfBatched(handle,m,n,
                        reinterpret_cast<hipComplex **>(Aarray),lda,
                        reinterpret_cast<hipComplex **>(TauArray),info,batchSize);
    hipDeviceSynchronize ();
    return sucess;
  }

}

#endif
