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

#ifndef AFQMC_BLAS_CUDA_CATCH_ALL_HPP
#define AFQMC_BLAS_CUDA_CATCH_ALL_HPP

#include<cassert>
#include "AFQMC/config.0.h"

// Currently available:
// Lvl-1: dot, axpy, scal
// Lvl-2: gemv
// Lvl-3: gemm

namespace qmc_cuda 
{

/*
  template<class ptrA, class ptrB>
  inline static void copy(int n, ptrA x, int incx, ptrB y, int incy)
  {
    print_stacktrace
    throw std::runtime_error("Error: Calling qmc_cuda::copy catch all."); 
  }
*/
  // dot Specializations
  template<class ptrA, class ptrB>
  inline static auto dot(int const n, ptrA x, int const incx, ptrB y, int const incy)
  {
    print_stacktrace
    throw std::runtime_error("Error: Calling qmc_cuda::dot catch all."); 
    using element = typename std::decay<ptrA>::type::value_type;
    return element(0);
  }

  // axpy Specializations
  template<typename T, class ptrA, class ptrB>
  inline static void axpy(int n, T const a,
                          ptrA x, int incx, 
                          ptrB y, int incy)
  {
    print_stacktrace
    throw std::runtime_error("Error: Calling qmc_cuda::axpy catch all."); 
  }

  // GEMV Specializations
  template<typename T, class ptrA, class ptrB, class ptrC>
  inline static void gemv(char Atrans, int M, int N,
                          T alpha,
                          ptrA A, int lda,
                          ptrB x, int incx,
                          T beta,
                          ptrC y, int incy)
  {
    std::cout<<" types: (gemv) "
             <<"  T: " <<typeid(alpha).name() <<"\n"
             <<"  ptrA: " <<typeid(A).name() <<"\n"
             <<"  ptrB: " <<typeid(x).name() <<"\n"
             <<"  ptrC: " <<typeid(y).name() <<std::endl;
    print_stacktrace
    throw std::runtime_error("Error: Calling qmc_cuda::gemv catch all."); 
  }

  template<typename T, class ptrA, class ptrB, class ptrC>
  inline static void gemm(char Atrans, char Btrans, int M, int N, int K,
                          T alpha,
                          ptrA A, int lda,
                          ptrB B, int ldb,
                          T beta,
                          ptrC C, int ldc)
  {
    std::cout<<" types: (gemm)"
             <<"  T: " <<typeid(alpha).name() <<"\n"
             <<"  ptrA: " <<typeid(A).name() <<"\n"
             <<"  ptrB: " <<typeid(B).name() <<"\n"
             <<"  ptrC: " <<typeid(C).name() <<std::endl;
    print_stacktrace
    throw std::runtime_error("Error: Unimplemented qmc_cuda::gemm with mixed pointers."); 
  }


  // Blas Extensions
  // geam  
  template<typename T, class ptrA, class ptrB, class ptrC>
  inline static void geam(char Atrans, char Btrans, int M, int N,
                         T const alpha,
                         ptrA A, int lda,
                         T const beta,
                         ptrB B, int ldb,
                         ptrC C, int ldc)
  {
    std::cout<<" types: (geam)"
             <<"  T: " <<typeid(alpha).name() <<"\n"
             <<"  ptrA: " <<typeid(A).name() <<"\n"
             <<"  ptrB: " <<typeid(B).name() <<"\n"
             <<"  ptrC: " <<typeid(C).name() <<std::endl;
    print_stacktrace
    throw std::runtime_error("Error: Calling qmc_cuda::geam catch all."); 
  }

  // dot extension 
  template<typename T, typename Q, class ptrA, class ptrB, class ptrC>
  inline static void adotpby(int const n, T const alpha, ptrA x, int const incx, ptrB y, int const incy, Q const beta, ptrC result)
  {
    std::cout<<" types: (adotpby)"
             <<"  T: " <<typeid(alpha).name() <<"\n"
             <<"  Q: " <<typeid(beta).name() <<"\n"
             <<"  ptrA: " <<typeid(x).name() <<"\n"
             <<"  ptrB: " <<typeid(y).name() <<"\n"
             <<"  ptrC: " <<typeid(result).name() <<std::endl;
    print_stacktrace
    throw std::runtime_error("Error: Calling qmc_cuda::adotpby catch all."); 
  }

  // axty
  template<class T, class ptrA, class ptrB>
  inline static void axty(int n,
                         T const alpha,
                         ptrA const x, int incx,
                         ptrB y, int incy)
  {
    std::cout<<" types: (axty)"
             <<"  T: " <<typeid(alpha).name() <<"\n"
             <<"  ptrA: " <<typeid(x).name() <<"\n"
             <<"  ptrB: " <<typeid(y).name() <<std::endl;
    print_stacktrace
    throw std::runtime_error("Error: Calling qmc_cuda::axty catch all."); 
  }

  // acAxpbB
  template<class T, class ptrA, class ptrx, class ptrB>
  inline static void acAxpbB(int m, int n,
                             T const alpha,
                             ptrA const A, int lda,
                             ptrx const x, int incx,
                             T const beta,
                             ptrB B, int ldb)
  {
    std::cout<<" types: (acAxpbB)"
             <<"  T: " <<typeid(alpha).name() <<"\n"
             <<"  ptrA: " <<typeid(A).name() <<"\n"
             <<"  ptrx: " <<typeid(x).name() <<"\n"
             <<"  ptrB: " <<typeid(B).name() <<std::endl;
    print_stacktrace
    throw std::runtime_error("Error: Calling qmc_cuda::acAxpbB catch all."); 
  }

  // adiagApy
  template<class T,
           class ptrA,
           class ptrB>
  inline static void adiagApy(int n,
                         T const alpha,
                         ptrA const A, int lda,
                         ptrB y, int incy)
  {
    std::cout<<" types: "
             <<"  T: " <<typeid(alpha).name() <<"\n"
             <<"  ptrA: " <<typeid(A).name() <<"\n"
             <<"  ptrB: " <<typeid(y).name() <<std::endl;
    print_stacktrace
    throw std::runtime_error("Error: Calling qmc_cuda::adiagApy catch all."); 
  }

  template<class T,
           class ptrA,
           class ptrB,
           class ptrC>
  inline static void gemmStridedBatched(char Atrans, char Btrans, int M, int N, int K,
                          T const alpha, ptrA const A, int lda, int strideA,
                          ptrB const B, int ldb, int strideB, T beta,
                          ptrC C, int ldc, int strideC, int batchSize)
  {
    std::cout<<" types: (gemmStridedBatched)"
             <<"  T: " <<typeid(alpha).name() <<"\n"
             <<"  ptrA: " <<typeid(A).name() <<"\n"
             <<"  ptrB: " <<typeid(B).name() <<"\n"
             <<"  ptrC: " <<typeid(C).name() <<std::endl;
    print_stacktrace
    throw std::runtime_error("Error: Calling qmc_cuda::gemmStridedBatched catch all."); 
  }

  template<class T,
           class ptrA,
           class ptrB,
           class ptrC>
  inline static void gemmBatched(char Atrans, char Btrans, int M, int N, int K,
                          T const alpha, ptrA const* A, int lda, 
                          ptrB const* B, int ldb, T beta,
                          ptrC * C, int ldc, int batchSize)
  {
    std::cout<<" types: (gemmBatched) "
             <<"  T: " <<typeid(alpha).name() <<"\n"
             <<"  ptrA: " <<typeid(A).name() <<"\n"
             <<"  ptrB: " <<typeid(B).name() <<"\n"
             <<"  ptrC: " <<typeid(C).name() <<std::endl;
    print_stacktrace
    throw std::runtime_error("Error: Calling qmc_cuda::gemmBatched catch all."); 
  }

}

#endif
