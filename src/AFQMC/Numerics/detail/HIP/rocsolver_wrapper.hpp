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

#ifndef HIPSOLVER_FUNCTIONDEFS_H
#define HIPSOLVER_FUNCTIONDEFS_H

#include<cassert>
#include <hip/hip_runtime.h>
#if defined(ENABLE_ROCM)
#include "rocsolver.h"
#endif
#include "AFQMC/Memory/HIP/hip_utilities.h"

// Abusing hip/rocm names
namespace rocsolver {

  using qmc_hip::hipblasOperation;
  using qmc_hip::rocsolverStatus_t;
  using qmc_hip::rocsolverHandle_t;

  // TODO: FDM These don't exist in rocsolver.
  inline rocsolverStatus_t
  rocsolver_getrf_bufferSize(rocsolverHandle_t handle,
                      int m,
                      int n,
                      float *A,
                      int lda,
                      int *Lwork )
  {
    rocsolverStatus_t sucess =
              rocsolver_sgetrf_bufferSize(handle,m,n,A,lda,Lwork);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_getrf_bufferSize(rocsolverHandle_t handle,
                      int m,
                      int n,
                      std::complex<float> *A,
                      int lda,
                      int *Lwork )
  {
    rocsolverStatus_t sucess =
              rocsolver_cgetrf_bufferSize(handle,m,n,
                                          reinterpret_cast<hipComplex *>(A),lda,Lwork);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_getrf_bufferSize(rocsolverHandle_t handle,
                      int m,
                      int n,
                      double *A,
                      int lda,
                      int *Lwork )
  {
    rocsolverStatus_t sucess =
              rocsolver_dgetrf_bufferSize(handle,m,n,A,lda,Lwork);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_getrf_bufferSize(rocsolverHandle_t handle,
                      int m,
                      int n,
                      std::complex<double> *A,
                      int lda,
                      int *Lwork )
  {
    rocsolverStatus_t sucess =
              rocsolver_zgetrf_bufferSize(handle,m,n,
                                          reinterpret_cast<hipDoubleComplex *>(A),lda,Lwork);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t rocsolver_getrf(rocsolverHandle_t handle,
                                   int m,
                                   int n,
                                   float *A,
                                   int lda,
                                   float *Work,
                                   int *devIpiv,
                                   int *devInfo )
  {
    rocsolverStatus_t sucess =
              rocsolver_sgetrf(handle,m,n,A,lda,Work,devIpiv,devInfo);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t rocsolver_getrf(rocsolverHandle_t handle,
                                   int m,
                                   int n,
                                   double *A,
                                   int lda,
                                   double *Work,
                                   int *devIpiv,
                                   int *devInfo )
  {
    rocsolverStatus_t sucess =
              rocsolver_dgetrf(handle,m,n,A,lda,Work,devIpiv,devInfo);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t rocsolver_getrf(rocsolverHandle_t handle,
                                   int m,
                                   int n,
                                   std::complex<float> *A,
                                   int lda,
                                   std::complex<float> *Work,
                                   int *devIpiv,
                                   int *devInfo )
  {
    rocsolverStatus_t sucess =
              rocsolver_cgetrf(handle,m,n,
                              reinterpret_cast<hipComplex *>(A),lda,
                              reinterpret_cast<hipComplex *>(Work),devIpiv,devInfo);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t rocsolver_getrf(rocsolverHandle_t handle,
                                   int m,
                                   int n,
                                   std::complex<double> *A,
                                   int lda,
                                   std::complex<double> *Work,
                                   int *devIpiv,
                                   int *devInfo )
  {
    rocsolverStatus_t sucess =
              rocsolver_zgetrf(handle,m,n,
                              reinterpret_cast<hipDoubleComplex *>(A),lda,
                              reinterpret_cast<hipDoubleComplex *>(Work),devIpiv,devInfo);
    hipDeviceSynchronize ();
    return sucess;
  }


  // getrs
  inline rocsolverStatus_t
  rocsolver_getrs(rocsolverHandle_t handle,
           hipblasOperation_t trans,
           int n,
           int nrhs,
           const float *A,
           int lda,
           const int *devIpiv,
           float *B,
           int ldb,
           int *devInfo )
  {
    rocsolverStatus_t sucess =
              rocsolver_sgetrs(handle,trans,n,nrhs,A,lda,devIpiv,B,ldb,devInfo);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_getrs(rocsolverHandle_t handle,
           hipblasOperation_t trans,
           int n,
           int nrhs,
           const double *A,
           int lda,
           const int *devIpiv,
           double *B,
           int ldb,
           int *devInfo )
  {
    rocsolverStatus_t sucess =
              rocsolver_dgetrs(handle,trans,n,nrhs,A,lda,devIpiv,B,ldb,devInfo);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_getrs(rocsolverHandle_t handle,
           hipblasOperation_t trans,
           int n,
           int nrhs,
           const std::complex<float> *A,
           int lda,
           const int *devIpiv,
           std::complex<float> *B,
           int ldb,
           int *devInfo )
  {
    rocsolverStatus_t sucess =
              rocsolver_cgetrs(handle,trans,n,nrhs,
                               reinterpret_cast<hipComplex const*>(A),lda,devIpiv,
                               reinterpret_cast<hipComplex *>(B),ldb,devInfo);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_getrs(rocsolverHandle_t handle,
           hipblasOperation_t trans,
           int n,
           int nrhs,
           const std::complex<double> *A,
           int lda,
           const int *devIpiv,
           std::complex<double> *B,
           int ldb,
           int *devInfo )
  {
    rocsolverStatus_t sucess =
              rocsolver_zgetrs(handle,trans,n,nrhs,
                               reinterpret_cast<hipDoubleComplex const*>(A),lda,devIpiv,
                               reinterpret_cast<hipDoubleComplex *>(B),ldb,devInfo);
    hipDeviceSynchronize ();
    return sucess;
  }

  //geqrf_bufferSize
  inline rocsolverStatus_t
  rocsolver_geqrf_bufferSize(rocsolverHandle_t handle,
                      int m,
                      int n,
                      float *A,
                      int lda,
                      int *Lwork )
  {
    rocsolverStatus_t sucess =
              rocsolver_sgeqrf_bufferSize(handle,m,n,A,lda,Lwork);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_geqrf_bufferSize(rocsolverHandle_t handle,
                      int m,
                      int n,
                      double *A,
                      int lda,
                      int *Lwork )
  {
    rocsolverStatus_t sucess =
              rocsolver_dgeqrf_bufferSize(handle,m,n,A,lda,Lwork);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_geqrf_bufferSize(rocsolverHandle_t handle,
                      int m,
                      int n,
                      std::complex<float> *A,
                      int lda,
                      int *Lwork )
  {
    rocsolverStatus_t sucess =
              rocsolver_cgeqrf_bufferSize(handle,m,n,reinterpret_cast<hipComplex *>(A),lda,Lwork);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_geqrf_bufferSize(rocsolverHandle_t handle,
                      int m,
                      int n,
                      std::complex<double> *A,
                      int lda,
                      int *Lwork )
  {
    rocsolverStatus_t sucess =
              rocsolver_zgeqrf_bufferSize(handle,m,n,reinterpret_cast<hipDoubleComplex *>(A),lda,Lwork);
    hipDeviceSynchronize ();
    return sucess;
  }

  //geqrf
  inline rocsolverStatus_t
  rocsolver_geqrf(rocsolverHandle_t handle,
           int m,
           int n,
           float *A,
           int lda,
           float *TAU,
           float *Workspace,
           int Lwork,
           int *devInfo )
  {
    rocsolverStatus_t sucess =
              rocsolver_sgeqrf(handle,m,n,A,lda,TAU,Workspace,Lwork,devInfo);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_geqrf(rocsolverHandle_t handle,
           int m,
           int n,
           double *A,
           int lda,
           double *TAU,
           double *Workspace,
           int Lwork,
           int *devInfo )
  {
    rocsolverStatus_t sucess =
              rocsolver_dgeqrf(handle,m,n,A,lda,TAU,Workspace,Lwork,devInfo);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_geqrf(rocsolverHandle_t handle,
           int m,
           int n,
           std::complex<float> *A,
           int lda,
           std::complex<float> *TAU,
           std::complex<float> *Workspace,
           int Lwork,
           int *devInfo )
  {
    rocsolverStatus_t sucess =
              rocsolver_cgeqrf(handle,m,n,reinterpret_cast<hipComplex *>(A),lda,
                               reinterpret_cast<hipComplex *>(TAU),reinterpret_cast<hipComplex *>(Workspace),Lwork,devInfo);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_geqrf(rocsolverHandle_t handle,
           int m,
           int n,
           std::complex<double> *A,
           int lda,
           std::complex<double> *TAU,
           std::complex<double> *Workspace,
           int Lwork,
           int *devInfo )
  {
    rocsolverStatus_t sucess =
              rocsolver_zgeqrf(handle,m,n,reinterpret_cast<hipDoubleComplex *>(A),lda,
                               reinterpret_cast<hipDoubleComplex *>(TAU),reinterpret_cast<hipDoubleComplex *>(Workspace),Lwork,devInfo);
    hipDeviceSynchronize ();
    return sucess;
  }



  //gqr_bufferSize
// MAM: There is an inconsistency between the dohipmentation online and the header files
// in hip.9.2.148 in the declaration of: rocsolvernXungqr_bufferSize and rocsolverXorgqr_bufferSize
//  I'm going to hack it right here, under the assumption that *tau is not used at all.
  inline rocsolverStatus_t
  rocsolver_gqr_bufferSize(
    rocsolverHandle_t handle,
    int m,
    int n,
    int k,
    const float *A,
    int lda,
    int *lwork)
  {
    rocsolverStatus_t sucess =
              rocsolver_sorgqr_bufferSize(handle,m,n,k,A,lda,A,lwork);
// HACK
//              rocsolver_sorgqr_bufferSize(handle,m,n,k,A,lda,lwork);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_gqr_bufferSize(
    rocsolverHandle_t handle,
    int m,
    int n,
    int k,
    const double *A,
    int lda,
    int *lwork)
  {
    rocsolverStatus_t sucess =
              rocsolver_dorgqr_bufferSize(handle,m,n,k,A,lda,A,lwork);
// HACK
//              rocsolver_dorgqr_bufferSize(handle,m,n,k,A,lda,lwork);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_gqr_bufferSize(
    rocsolverHandle_t handle,
    int m,
    int n,
    int k,
    const std::complex<float> *A,
    int lda,
    int *lwork)
  {
    rocsolverStatus_t sucess =
              rocsolver_cungqr_bufferSize(handle,m,n,k,reinterpret_cast<hipComplex const *>(A),lda,
                                          reinterpret_cast<hipComplex const *>(A),lwork);
// HACK
//                                          lwork);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_gqr_bufferSize(
    rocsolverHandle_t handle,
    int m,
    int n,
    int k,
    const std::complex<double> *A,
    int lda,
    int *lwork)
  {
    rocsolverStatus_t sucess =
              rocsolver_zungqr_bufferSize(handle,m,n,k,reinterpret_cast<hipDoubleComplex const*>(A),lda,
                                          reinterpret_cast<hipDoubleComplex const*>(A),lwork);
// HACK
//                                          lwork);
    hipDeviceSynchronize ();
    return sucess;
  }

  //gqr
  inline rocsolverStatus_t
  rocsolver_gqr(
    rocsolverHandle_t handle,
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
    rocsolverStatus_t sucess =
              rocsolver_sorgqr(handle,m,n,k,A,lda,tau,work,lwork,devInfo);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_gqr(
    rocsolverHandle_t handle,
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
    rocsolverStatus_t sucess =
              rocsolver_dorgqr(handle,m,n,k,A,lda,tau,work,lwork,devInfo);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_gqr(
    rocsolverHandle_t handle,
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
    rocsolverStatus_t sucess =
              rocsolver_cungqr(handle,m,n,k,reinterpret_cast<hipComplex *>(A),lda,
                               reinterpret_cast<hipComplex const *>(tau),reinterpret_cast<hipComplex *>(work),lwork,devInfo);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_gqr(
    rocsolverHandle_t handle,
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
    rocsolverStatus_t sucess =
              rocsolver_zungqr(handle,m,n,k,reinterpret_cast<hipDoubleComplex *>(A),lda,
                               reinterpret_cast<hipDoubleComplex const *>(tau),reinterpret_cast<hipDoubleComplex *>(work),lwork,devInfo);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t rocsolver_gqr_strided(
    rocsolverHandle_t handle,
    const int m,
    const int n,
    const int k,
    std::complex<double> *A,
    const int lda,
    const int Astride,
    const std::complex<double> *tau,
    const int tstride,
    std::complex<double> *work,
    int lwork,
    int *devInfo,
    int batchsize) {

    using qmc_hip::afqmc_hip_streams;
    if(afqmc_hip_streams.size() < batchsize) {
      int n0=afqmc_hip_streams.size();
      for(int n=n0; n<batchsize; n++) {
        afqmc_hip_streams.emplace_back(hipStream_t{});
        hipStreamCreate(&(afqmc_hip_streams.back()));
      }
    }

    hipStream_t s0;
    qmc_hip::rocsolver_check(rocsolverGetStream(handle,&s0), "rocsolverGetStream");

    for(int i=0; i<batchsize; i++) {
      qmc_hip::rocsolver_check(rocsolver_setStream(handle,afqmc_hip_streams[i]), "rocsolver_setStream");
      rocsolverStatus_t sucess =
              rocsolver_zungqr(handle,m,n,k,reinterpret_cast<hipDoubleComplex *>(A)+i*Astride,lda,
                        reinterpret_cast<hipDoubleComplex const *>(tau)+i*tstride,
                        reinterpret_cast<hipDoubleComplex *>(work)+i*lwork,lwork,devInfo+i);
      qmc_hip::rocsolver_check(sucess,"rocsolver_gqr_strided");
    }
    qmc_hip::hip_check(hipDeviceSynchronize(),"rocsolver_gqr_strided");
    qmc_hip::hip_check(hipGetLastError(),"rocsolver_gqr_strided");
    qmc_hip::rocsolver_check(rocsolver_setStream(handle,s0), "rocsolver_setStream");

    return CUSOLVER_STATUS_SUCCESS;

  }

  //gesvd_bufferSize
  inline rocsolverStatus_t
  rocsolver_gesvd_bufferSize(rocsolverHandle_t handle,
                      int m,
                      int n,
                      float* A,
                      int *Lwork )
  {
    rocsolverStatus_t sucess =
              rocsolver_sgesvd_bufferSize(handle,m,n,Lwork);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_gesvd_bufferSize(rocsolverHandle_t handle,
                      int m,
                      int n,
                      double* A,
                      int *Lwork )
  {
    rocsolverStatus_t sucess =
              rocsolver_dgesvd_bufferSize(handle,m,n,Lwork);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_gesvd_bufferSize(rocsolverHandle_t handle,
                      int m,
                      int n,
                      std::complex<float>* A,
                      int *Lwork )
  {
    rocsolverStatus_t sucess =
              rocsolver_cgesvd_bufferSize(handle,m,n,Lwork);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_gesvd_bufferSize(rocsolverHandle_t handle,
                      int m,
                      int n,
                      std::complex<double>* A,
                      int *Lwork )
  {
    rocsolverStatus_t sucess =
              rocsolver_zgesvd_bufferSize(handle,m,n,Lwork);
    hipDeviceSynchronize ();
    return sucess;
  }

  //gesvd
  inline rocsolverStatus_t
  rocsolver_gesvd(rocsolverHandle_t handle,
    signed char jobu,
    signed char jobvt,
    int m,
    int n,
    float *A,
    int lda,
    float *S,
    float *U,
    int ldu,
    float *VT,
    int ldvt,
    float *work,
    int lwork,
    int *devInfo)
  {
    rocsolverStatus_t sucess =
              rocsolver_sgesvd(handle,jobu,jobvt,m,n,A,lda,S,U,ldu,VT,ldvt,work,lwork,nullptr,devInfo);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_gesvd(rocsolverHandle_t handle,
    signed char jobu,
    signed char jobvt,
    int m,
    int n,
    double *A,
    int lda,
    double *S,
    double *U,
    int ldu,
    double *VT,
    int ldvt,
    double *work,
    int lwork,
    int *devInfo)
  {
    rocsolverStatus_t sucess =
              rocsolver_dgesvd(handle,jobu,jobvt,m,n,A,lda,S,U,ldu,VT,ldvt,work,lwork,nullptr,devInfo);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_gesvd(rocsolverHandle_t handle,
    signed char jobu,
    signed char jobvt,
    int m,
    int n,
    std::complex<float> *A,
    int lda,
    float *S,
    std::complex<float> *U,
    int ldu,
    std::complex<float> *VT,
    int ldvt,
    std::complex<float> *work,
    int lwork,
    int *devInfo)
  {
    rocsolverStatus_t sucess =
             rocsolver_cgesvd(handle,jobu,jobvt,m,n,
                               reinterpret_cast<hipComplex *>(A),lda,S,
                               reinterpret_cast<hipComplex *>(U),ldu,
                               reinterpret_cast<hipComplex *>(VT),ldvt,
                               reinterpret_cast<hipComplex *>(work),lwork,
                               nullptr,devInfo);
    hipDeviceSynchronize ();
    return sucess;
  }

  inline rocsolverStatus_t
  rocsolver_gesvd(rocsolverHandle_t handle,
    signed char jobu,
    signed char jobvt,
    int m,
    int n,
    std::complex<double> *A,
    int lda,
    double *S,
    std::complex<double> *U,
    int ldu,
    std::complex<double> *VT,
    int ldvt,
    std::complex<double> *work,
    int lwork,
    int *devInfo)
  {
    rocsolverStatus_t sucess =
              rocsolver_zgesvd(handle,jobu,jobvt,m,n,
                               reinterpret_cast<hipDoubleComplex *>(A),lda,S,
                               reinterpret_cast<hipDoubleComplex *>(U),ldu,
                               reinterpret_cast<hipDoubleComplex *>(VT),ldvt,
                               reinterpret_cast<hipDoubleComplex *>(work),lwork,
                               nullptr,devInfo);
    hipDeviceSynchronize ();
    return sucess;
  }

}

#endif
