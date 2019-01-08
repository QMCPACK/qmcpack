////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
// Alfredo Correa, correaa@llnl.gov 
//    Lawrence Livermore National Laboratory 
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////


#ifndef  AFQMC_DENSITY_MATRIX_HPP 
#define  AFQMC_DENSITY_MATRIX_HPP 

#include "AFQMC/config.h"
#include "AFQMC/Numerics/ma_operations.hpp"
#include <Utilities/FairDivide.h>

namespace qmcplusplus
{

namespace afqmc 
{

namespace SlaterDeterminantOperations
{

namespace base 
{

/*
 * Calculates the 1-body mixed density matrix:
 *   < A | c+i cj | B > / <A|B> = conj(A) * ( T(B) * conj(A) )^-1 * T(B) 
 *   If compact == True, returns [NEL x M] matrix:
 *   < A | c+i cj | B > / <A|B> = ( T(B) * conj(A) )^-1 * T(B) 
 * Parameters:
 *  - hermA = conjugateTranspose(A)
 *  - B
 *  - C = < A | c+i cj | B > / <A|B>
 *  - T1: [ NEL x NEL ] work matrix   
 *  - T2: (only used if compact = False) [ NEL x M ] work matrix   
 *  - IWORK: [ N ] integer biffer for invert. Dimensions must be at least NEL. 
 *  - WORK: [ >=NEL ] Work space for invert. Dimensions must be at least NEL.   
 *  - compact (default = True)
 *  returns:
 *  - <A|B> = det[ T(conj(A)) * B ]  
 */
// Serial Implementation
template< class Tp,
          class MatA,
          class MatB,
          class MatC,
          class Mat1,
          class Mat2,
          class IBuffer,
          class TBuffer 
        >
inline Tp MixedDensityMatrix(const MatA& hermA, const MatB& B, MatC&& C, Mat1&& T1, Mat2&& T2, IBuffer& IWORK, TBuffer& WORK, bool compact=true)
{
  // check dimensions are consistent
  assert( hermA.shape()[1] == B.shape()[0] );
  assert( hermA.shape()[0] == B.shape()[1] );
  assert( hermA.shape()[0] == T1.shape()[0] );
  assert( B.shape()[1] == T1.shape()[1] );
  if(compact) {
    assert( C.shape()[0] == T1.shape()[1] );
    assert( C.shape()[1] == B.shape()[0] );
  } else {
    assert( T2.shape()[1] == B.shape()[0] );
    assert( T2.shape()[0] == T1.shape()[1] );
    assert( C.shape()[0] == hermA.shape()[1] );
    assert( C.shape()[1] == T2.shape()[1] );
  }

  using ma::T;

  // T(B)*conj(A) 
  ma::product(hermA,B,std::forward<Mat1>(T1));  

  // NOTE: Using C as temporary 
  // T1 = T1^(-1)
  Tp ovlp = static_cast<Tp>(ma::invert(std::forward<Mat1>(T1),IWORK,WORK));
  if(ovlp == Tp(0.0)) return Tp(0.0); // don't bother calculating further

  if(compact) {

    // C = T(T1) * T(B)
    ma::product(T(T1),T(B),std::forward<MatC>(C)); 

  } else {

    // T2 = T(T1) * T(B)
    ma::product(T(T1),T(B),std::forward<Mat2>(T2)); 

    // C = conj(A) * T2
    ma::product(T(hermA),T2,std::forward<MatC>(C));

  }

  return ovlp;
}

template< class Tp,
          class integer,
          class MatA,
          class MatB,
          class MatC,
          class MatD,
          class Mat1,
          class Mat2,
          class Mat3,
          class IBuffer,
          class TBuffer 
        >
inline Tp MixedDensityMatrixForWoodbury(const MatA& hermA, const MatB& B, MatC&& C, MatD&& QQ0, integer *ref, Mat1&& TNN, Mat2&& TAB, Mat3&& TNM, IBuffer& IWORK, TBuffer& WORK, bool compact=true)
{
  // check dimensions are consistent
  int NEL = B.shape()[1];
  assert( hermA.shape()[1] == B.shape()[0] );
  assert( hermA.shape()[0] == TAB.shape()[0] );
  assert( B.shape()[1] == TAB.shape()[1] );
  assert( B.shape()[1] == TNN.shape()[0] );
  assert( B.shape()[1] == TNN.shape()[1] );
  assert( hermA.shape()[0] == QQ0.shape()[0] );
  assert( B.shape()[1] == QQ0.shape()[1] );
  if(compact) {
    assert( C.shape()[0] == TNN.shape()[1] );
    assert( C.shape()[1] == B.shape()[0] );
  } else {
    assert( TNM.shape()[1] == B.shape()[0] );
    assert( TNM.shape()[0] == TNN.shape()[1] );
    assert( C.shape()[0] == hermA.shape()[1] );
    assert( C.shape()[1] == TNM.shape()[1] );
  }

  using ma::T;

  // TAB = herm(A)*B
  ma::product(hermA,B,std::forward<Mat2>(TAB));  

  // TNN = TAB[ref,:] 
  for(int i=0; i<NEL; i++)
    std::copy_n(std::addressof(*TAB[*(ref+i)].origin()),NEL,std::addressof(*TNN[i].origin()));

  // TNN = TNN^(-1)
  Tp ovlp = static_cast<Tp>(ma::invert(std::forward<Mat1>(TNN),IWORK,WORK));
  if(ovlp == Tp(0.0)) return Tp(0.0); // don't bother calculating further

  // QQ0 = TAB * inv(TNN) 
  ma::product(TAB,TNN,std::forward<MatD>(QQ0));

  if(compact) {

    // C = T(TNN) * T(B)
    ma::product(T(TNN),T(B),std::forward<MatC>(C)); 

  } else {

    // TNM = T(TNN) * T(B)
    ma::product(T(TNN),T(B),std::forward<Mat3>(TNM)); 

    // C = conj(A) * TNM
    ma::product(T(hermA),TNM,std::forward<MatC>(C));

  }

  return ovlp;
}

template< class Tp,
          class integer,
          class MatA,
          class MatB,
          class MatC,
          class Mat1,
          class Mat2,
          class Mat3,
          class IBuffer,
          class TBuffer 
        >
inline Tp MixedDensityMatrixFromConfiguration(const MatA& hermA, const MatB& B, MatC&& C, integer *ref, Mat1&& TNN, Mat2&& TAB, Mat3&& TNM, IBuffer& IWORK, TBuffer& WORK, bool compact=true)
{
  // check dimensions are consistent
  int NEL = B.shape()[1];
  assert( hermA.shape()[1] == B.shape()[0] );
  assert( hermA.shape()[0] == TAB.shape()[0] );
  assert( B.shape()[1] == TAB.shape()[1] );
  assert( B.shape()[1] == TNN.shape()[0] );
  assert( B.shape()[1] == TNN.shape()[1] );
  if(compact) {
    assert( C.shape()[0] == TNN.shape()[1] );
    assert( C.shape()[1] == B.shape()[0] );
  } else {
    assert( TNM.shape()[1] == B.shape()[0] );
    assert( TNM.shape()[0] == TNN.shape()[1] );
    assert( C.shape()[0] == hermA.shape()[1] );
    assert( C.shape()[1] == TNM.shape()[1] );
  }

  using ma::T;

  // TAB = herm(A)*B
  ma::product(hermA,B,std::forward<Mat2>(TAB));  

  // TNN = TAB[ref,:] 
  for(int i=0; i<NEL; i++)
    std::copy_n(std::addressof(*TAB[*(ref+i)].origin()),NEL,std::addressof(*TNN[i].origin()));

  // TNN = TNN^(-1)
  Tp ovlp = static_cast<Tp>(ma::invert(std::forward<Mat1>(TNN),IWORK,WORK));
  if(ovlp == Tp(0.0)) return Tp(0.0); // don't bother calculating further

  if(compact) {

    // C = T(TNN) * T(B)
    ma::product(T(TNN),T(B),std::forward<MatC>(C)); 

  } else {

    // TNM = T(TNN) * T(B)
    ma::product(T(TNN),T(B),std::forward<Mat3>(TNM)); 

    // C = conj(A) * TNM
    ma::product(T(hermA),TNM,std::forward<MatC>(C));

  }

  return ovlp;
}

/*
 * Calculates the 1-body mixed density matrix:
 *   < A | c+i cj | B > / <A|B> = conj(A) * ( T(B) * conj(A) )^-1 * T(B) 
 *   If compact == True, returns [NEL x M] matrix:
 *   < A | c+i cj | B > / <A|B> = ( T(B) * conj(A) )^-1 * T(B) 
 * Parameters:
 *  - A = A
 *  - B
 *  - C = < A | c+i cj | B > / <A|B>
 *  - T1: [ NEL x NEL ] work matrix   
 *  - T2: (only used if compact = False) [ NEL x M ] work matrix   
 *  - IWORK: [ N ] integer biffer for invert. Dimensions must be at least NEL. 
 *  - WORK: [ >=NEL ] Work space for invert. Dimensions must be at least NEL.   
 *  - compact (default = True)
 *  returns:
 *  - <A|B> = det[ T(conj(A)) * B ]  
 */
// Serial Implementation
template< class Tp,
          class MatA,
          class MatB,
          class MatC,
          class Mat1,
          class Mat2,
          class IBuffer,
          class TBuffer 
        >
inline Tp MixedDensityMatrix_noHerm(const MatA& A, const MatB& B, MatC&& C, Mat1&& T1, Mat2&& T2, IBuffer& IWORK, TBuffer& WORK, bool compact=true)
{
  // check dimensions are consistent
  assert( A.shape()[0] == B.shape()[0] );
  assert( A.shape()[1] == B.shape()[1] );
  assert( A.shape()[1] == T1.shape()[0] );
  assert( B.shape()[1] == T1.shape()[1] );
  if(compact) {
    assert( C.shape()[0] == T1.shape()[1] );
    assert( C.shape()[1] == B.shape()[0] );
  } else {
    assert( T2.shape()[1] == B.shape()[0] );
    assert( T2.shape()[0] == T1.shape()[1] );
    assert( C.shape()[0] == A.shape()[0] );
    assert( C.shape()[1] == T2.shape()[1] );
  }

  using ma::T;
  using ma::H;

  // T1 = H(A)*B
  ma::product(H(A),B,std::forward<Mat1>(T1));  

  // NOTE: Using C as temporary 
  // T1 = T1^(-1)
  Tp ovlp = static_cast<Tp>(ma::invert(std::forward<Mat1>(T1),IWORK,WORK));
  if(ovlp == Tp(0.0)) return Tp(0.0); // don't bother calculating further

  if(compact) {

    // C = T(T1) * T(B)
    ma::product(T(T1),T(B),std::forward<MatC>(C)); 

  } else {

    // T2 = T1 * H(A) 
    ma::product(T1,H(A),std::forward<Mat2>(T2)); 

    // C = T( B * T2) = T(T2) * T(B)
    ma::product(T(T2),T(B),std::forward<MatC>(C));

  }

  return ovlp;
}


/*
 * Returns the overlap of 2 Slater determinants:  <A|B> = det[ T(B) * conj(A) ]  
 * Parameters:
 *  - hermA = conjugateTranspose(A)
 *  - B
 *  - IWORK: [ M ] integer work matrix   
 *  - T1: [ NEL x NEL ] work matrix   
 *  returns:
 *  - <A|B> = det[ T(conj(A)) * B ]  
 */
// Serial Implementation
template< class Tp,
          class MatA,
          class MatB,
          class Mat,
          class Buffer,
          class IBuffer
        >
inline Tp Overlap(const MatA& hermA, const MatB& B, Mat&& T1, IBuffer& IWORK, Buffer& WORK)
{
  // check dimensions are consistent
  assert( hermA.shape()[1] == B.shape()[0] );
  assert( hermA.shape()[0] == B.shape()[1] );
  assert( hermA.shape()[0] == T1.shape()[0] );
  assert( B.shape()[1] == T1.shape()[1] );

  using ma::T;

  // T(B)*conj(A) 
  ma::product(hermA,B,std::forward<Mat>(T1));   

  return static_cast<Tp>(ma::determinant(std::forward<Mat>(T1),IWORK,WORK));
}

template< class Tp,
          class integer,
          class MatA,
          class MatB,
          class MatC,
          class MatD,
          class MatE,
          class IBuffer,
          class TBuffer
        >
inline Tp OverlapForWoodbury(const MatA& hermA, const MatB& B, MatC&& QQ0, integer *ref, 
                             MatD && TNN, MatE && TMN, IBuffer& IWORK, TBuffer& WORK)
{
  // check dimensions are consistent
  int NEL = B.shape()[1];
  assert( hermA.shape()[1] == B.shape()[0] );
  assert( hermA.shape()[0] == TMN.shape()[0] );
  assert( B.shape()[1] == TMN.shape()[1] );
  assert( B.shape()[1] == TNN.shape()[0] );
  assert( B.shape()[1] == TNN.shape()[1] );
  assert( hermA.shape()[0] == QQ0.shape()[0] );
  assert( B.shape()[1] == QQ0.shape()[1] );

  using ma::T;

  // TMN = herm(A)*B 
  ma::product(hermA,B,std::forward<MatD>(TMN));

  // TNN = TMN[ref,:]
  for(int i=0; i<NEL; i++)
    std::copy_n(std::addressof(*TMN[*(ref+i)].origin()),NEL,std::addressof(*TNN[i].origin())); 
 
  // TNN -> inv(TNN)
  Tp res =  static_cast<Tp>(ma::invert(std::forward<MatD>(TNN),IWORK,WORK));

  // QQ0 = TMN * inv(TNN) 
  ma::product(TMN,TNN,std::forward<MatC>(QQ0));  

  return res;
}

/*
 * Returns the overlap of 2 Slater determinants:  <A|B> = det[ T(B) * conj(A) ]  
 * Parameters:
 *  - A = conj(A)
 *  - B
 *  - IWORK: [ M ] integer work matrix   
 *  - T1: [ NEL x NEL ] work matrix   
 *  returns:
 *  - <A|B> = det[ T(conj(A)) * B ]  
 */
// Serial Implementation
template< class Tp,
          class MatA,
          class MatB,
          class Mat,
          class Buffer,
          class IBuffer
        >
inline Tp Overlap_noHerm(const MatA& A, const MatB& B, Mat&& T1, IBuffer& IWORK, Buffer& WORK)
{
  // check dimensions are consistent
  assert( A.shape()[0] == B.shape()[0] );
  assert( A.shape()[1] == B.shape()[1] );
  assert( A.shape()[1] == T1.shape()[0] );
  assert( B.shape()[1] == T1.shape()[1] );

  using ma::T;
  using ma::H;

  // T(B)*conj(A) 
  ma::product(H(A),B,std::forward<Mat>(T1));    

  return static_cast<Tp>(ma::determinant(std::forward<Mat>(T1),IWORK,WORK));
}

} // namespace base

namespace shm 
{

/*
 * Calculates the 1-body mixed density matrix:
 *   < A | c+i cj | B > / <A|B> = conj(A) * ( T(B) * conj(A) )^-1 * T(B) 
 *   If compact == True, returns [NEL x M] matrix:
 *   < A | c+i cj | B > / <A|B> = ( T(B) * conj(A) )^-1 * T(B) 
 * Parameters:
 *  - hermA = conjugateTranspose(A)
 *  - B
 *  - C = < A | c+i cj | B > / <A|B>
 *  - T1: [ NEL x NEL ] work matrix   
 *  - T2: (only used if compact = False) [ NEL x M ] work matrix   
 *  - IWORK: [ N ] integer biffer for invert. Dimensions must be at least NEL. 
 *  - WORK: [ >=NEL ] Work space for invert. Dimensions must be at least NEL.   
 *  - compact (default = True)
 *  returns:
 *  - <A|B> = det[ T(B) * conj(A) ]  
 */
// Serial Implementation
template< class Tp,
          class MatA,
          class MatB,
          class MatC,
          class Mat,
          class IBuffer,
          class TBuffer,
          class communicator 
        >
inline Tp MixedDensityMatrix(const MatA& hermA, const MatB& B, MatC&& C, Mat&& T1, Mat&& T2, IBuffer& IWORK, TBuffer& WORK, communicator& comm, bool compact=true)
{
  // check dimensions are consistent
  assert( hermA.shape()[1] == B.shape()[0] );
  assert( hermA.shape()[0] == B.shape()[1] );
  assert( hermA.shape()[0] == T1.shape()[0] );
  assert( B.shape()[1] == T1.shape()[1] );
  if(compact) {
    assert( C.shape()[0] == T1.shape()[1] );
    assert( C.shape()[1] == B.shape()[0] );
  } else {
    assert( T2.shape()[1] == B.shape()[0] );
    assert( T2.shape()[0] == T1.shape()[1] );
    assert( C.shape()[0] == hermA.shape()[1] );
    assert( C.shape()[1] == T2.shape()[1] );
  }

  using ma::T;

  int N0,Nn,sz=B.shape()[1];
  std::tie(N0,Nn) = FairDivideBoundary(comm.rank(),sz,comm.size());

  // T(B)*conj(A) 
  if(N0!=Nn)
    ma::product(hermA,
              B(B.extension(0),{N0,Nn}),
              T1(T1.extension(0),{N0,Nn}));  

  comm.barrier();

  // NOTE: Using C as temporary 
  // T1 = T1^(-1)
  Tp ovlp=Tp(0.);
  if(comm.rank()==0)
   ovlp = static_cast<Tp>(ma::invert(std::forward<Mat>(T1),IWORK,WORK));
  comm.broadcast_n(&ovlp,1,0);  
  if(ovlp == Tp(0.0)) return Tp(0.0); // don't bother calculating further

  if(compact) {

    // C = T(T1) * T(B)
    //ma::product(T1.sliced(N0,Nn),
    //            T(B),
    //            C.sliced(N0,Nn)); 
    if(N0!=Nn)
      ma::product(T(T1(T1.extension(0),{N0,Nn})),
                T(B),
                C.sliced(N0,Nn)); 

  } else {

    // T2 = T(T1) * T(B)
    //ma::product(T1.sliced(N0,Nn),
    //            T(B),
    //            T2.sliced(N0,Nn)); 
    if(N0!=Nn)
      ma::product(T(T1(T1.extension(0),{N0,Nn})),
                T(B),
                T2.sliced(N0,Nn)); 

    comm.barrier();
    
    sz=T2.shape()[1];
    std::tie(N0,Nn) = FairDivideBoundary(comm.rank(),sz,comm.size());

    // C = conj(A) * T2
    if(N0!=Nn)
      ma::product(T(hermA),
                T2(T2.extension(0),{N0,Nn}),
                C(C.extension(0),{N0,Nn}));

  }

  comm.barrier();

  return ovlp;
}


/*
 * Returns the overlap of 2 Slater determinants:  <A|B> = det[ T(conj(A)) * B ]  
 * Parameters:
 *  - hermA = conjugateTranspose(A)
 *  - B
 *  - IWORK: [ M ] integer work matrix   
 *  - T1: [ NEL x NEL ] work matrix   
 *  returns:
 *  - <A|B> = det[ hermA * B ]  
 */
// Serial Implementation
template< class Tp,
          class MatA,
          class MatB,
          class Mat,
          class IBuffer,
          class Buffer,
          class communicator
        >
inline Tp Overlap(const MatA& hermA, const MatB& B, Mat&& T1, IBuffer& IWORK, Buffer& WORK, communicator& comm)
{
  // check dimensions are consistent
  assert( hermA.shape()[1] == B.shape()[0] );
  assert( hermA.shape()[0] == B.shape()[1] );
  assert( hermA.shape()[0] == T1.shape()[0] );
  assert( B.shape()[1] == T1.shape()[1] );

  using ma::T;

  int N0,Nn,sz = B.shape()[1];
  std::tie(N0,Nn) = FairDivideBoundary(comm.rank(),sz,comm.size());

  // T(B)*conj(A) 
  if(N0!=Nn)
    ma::product(hermA,
              B(B.extension(0),{N0,Nn}),
              T1(T1.extension(0),{N0,Nn}));

  comm.barrier();

  Tp ovlp=Tp(0.);
  if(comm.rank()==0)
   ovlp = static_cast<Tp>(ma::determinant(std::forward<Mat>(T1),IWORK,WORK));
  comm.broadcast_n(&ovlp,1,0);  

  return ovlp; 
}

// Serial Implementation
template< class Tp,
          class integer,
          class MatA,
          class MatB,
          class MatC,
          class MatD,
          class MatE,
          class IBuffer,
          class TBuffer,
          class communicator
        >
inline Tp OverlapForWoodbury(const MatA& hermA, const MatB& B, MatC&& QQ0, integer *ref, MatD&& TNN, MatE& TMN, IBuffer& IWORK, TBuffer& WORK, communicator& comm)
{
  // check dimensions are consistent
  int NEL = B.shape()[1];
  assert( hermA.shape()[1] == B.shape()[0] );
  assert( hermA.shape()[0] == TMN.shape()[0] );
  assert( B.shape()[1] == TMN.shape()[1] );
  assert( B.shape()[1] == TNN.shape()[0] );
  assert( B.shape()[1] == TNN.shape()[1] );
  assert( hermA.shape()[0] == QQ0.shape()[0] );
  assert( B.shape()[1] == QQ0.shape()[1] );

  using ma::T;

  int N0,Nn,sz = B.shape()[1];
  std::tie(N0,Nn) = FairDivideBoundary(comm.rank(),sz,comm.size());

  // T(B)*conj(A) 
  if(N0!=Nn)
    ma::product(hermA,
              B(B.extension(0),{N0,Nn}),
              TMN(TMN.extension(0),{N0,Nn}));
  comm.barrier();
  Tp ovlp=Tp(0.);
  if(comm.rank()==0) {
    for(int i=0; i<NEL; i++)
      std::copy_n(std::addressof(*TMN[*(ref+i)].origin()),NEL,std::addressof(*TNN[i].origin()));
   ovlp = static_cast<Tp>(ma::invert(std::forward<MatD>(TNN),IWORK,WORK));
  }
  comm.broadcast_n(&ovlp,1,0);

  int M0,Mn;
  sz = TMN.shape()[0];
  std::tie(M0,Mn) = FairDivideBoundary(comm.rank(),sz,comm.size());

  // QQ0 = TMN * inv(TNN) 
  ma::product(TMN.sliced(M0,Mn),TNN,
              QQ0({M0,Mn},QQ0.extension(1))); //.sliced(M0,Mn)); 
              //QQ0.sliced(M0,Mn)); 
  comm.barrier();
  return ovlp;
}

template< class Tp,
          class integer,
          class MatA,
          class MatB,
          class MatC,
          class MatD,
          class Mat1,
          class Mat2,
          class Mat3,
          class IBuffer,
          class TBuffer,
          class communicator 
        >
inline Tp MixedDensityMatrixForWoodbury(const MatA& hermA, const MatB& B, MatC&& C, MatD&& QQ0, integer *ref, Mat1&& TNN, Mat2&& TAB, Mat3&& TNM, IBuffer& IWORK, TBuffer& WORK, communicator& comm, bool compact=true)
{
  // check dimensions are consistent
  int NEL = B.shape()[1];
  assert( hermA.shape()[1] == B.shape()[0] );
  assert( hermA.shape()[0] == TAB.shape()[0] );
  assert( B.shape()[1] == TAB.shape()[1] );
  assert( B.shape()[1] == TNN.shape()[0] );
  assert( B.shape()[1] == TNN.shape()[1] );
  assert( hermA.shape()[0] == QQ0.shape()[0] );
  assert( B.shape()[1] == QQ0.shape()[1] );
  if(compact) {
    assert( C.shape()[0] == TNN.shape()[1] );
    assert( C.shape()[1] == B.shape()[0] );
  } else {
    assert( TNM.shape()[1] == B.shape()[0] );
    assert( TNM.shape()[0] == TNN.shape()[1] );
    assert( C.shape()[0] == hermA.shape()[1] );
    assert( C.shape()[1] == TNM.shape()[1] );
  }

  using ma::T;

  int N0,Nn;
  std::tie(N0,Nn) = FairDivideBoundary(comm.rank(),NEL,comm.size());

  // TAB = herm(A)*B
  if(N0!=Nn) {
    ma::product(hermA,
              B(B.extension(0),{N0,Nn}),
              TAB(TAB.extension(0),{N0,Nn}));  

    // TNN = TAB[ref,:] 
    for(int i=0; i<NEL; i++)
      std::copy_n(std::addressof(*TAB[*(ref+i)].origin())+N0,Nn-N0,std::addressof(*TNN[i].origin())+N0);
  }

  comm.barrier();

  // TNN = TNN^(-1)
  Tp ovlp=Tp(0.);
  if(comm.rank()==0)
    ovlp = static_cast<Tp>(ma::invert(std::forward<Mat1>(TNN),IWORK,WORK));
  comm.broadcast_n(&ovlp,1,0);
  if(ovlp == Tp(0.0)) return Tp(0.0); // don't bother calculating further

  int P0,Pn;
  std::tie(P0,Pn) = FairDivideBoundary(comm.rank(),int(TAB.shape()[0]),comm.size());

  // QQ0 = TAB * inv(TNN) 
  if(P0!=Pn)  
    ma::product(TAB.sliced(P0,Pn),
              TNN,
              QQ0({P0,Pn},QQ0.extension(1)));  
              //QQ0.sliced(P0,Pn));
  if(compact) {

    // C = T(TNN) * T(B)
    if(N0!=Nn)
      ma::product(T(TNN(TNN.extension(0),{N0,Nn})),
                  T(B),
                  C({N0,Nn},C.extension(1))); 

  } else {

    // TNM = T(TNN) * T(B)
    if(N0!=Nn)    
      ma::product(T(TNN(TNN.extension(0),{N0,Nn})),
                  T(B),
                  TNM.sliced(N0,Nn)); 

    int sz=TNM.shape()[1];
    std::tie(N0,Nn) = FairDivideBoundary(comm.rank(),sz,comm.size());   
    comm.barrier();

    // C = conj(A) * TNM
    ma::product(T(hermA),
                TNM(TNM.extension(0),{N0,Nn}),
                C(C.extension(0), {N0,Nn})); 

  }

  comm.barrier();

  return ovlp;
}

} // namespace shm 

} // namespace SlaterDeterminantOperations

} // namespace afqmc

} // namespace qmcplusplus 

#endif

