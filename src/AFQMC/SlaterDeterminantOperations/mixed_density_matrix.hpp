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
inline void MixedDensityMatrix(const MatA& hermA, const MatB& B, MatC&& C, Tp* ovlp, Mat1&& T1, Mat2&& T2, IBuffer& IWORK, TBuffer& WORK, bool compact=true)
{
  // check dimensions are consistent
  assert( hermA.size(1) == B.size(0) );
  assert( hermA.size(0) == B.size(1) );
  assert( hermA.size(0) == T1.size(0) );
  assert( B.size(1) == T1.size(1) );
  if(compact) {
    assert( C.size(0) == T1.size(1) );
    assert( C.size(1) == B.size(0) );
  } else {
    assert( T2.size(1) == B.size(0) );
    assert( T2.size(0) == T1.size(1) );
    assert( C.size(0) == hermA.size(1) );
    assert( C.size(1) == T2.size(1) );
  }

  using ma::T;

  // T(B)*conj(A) 
  ma::product(hermA,B,std::forward<Mat1>(T1));  

  // NOTE: Using C as temporary 
  // T1 = T1^(-1)
  ma::invert(std::forward<Mat1>(T1),IWORK,WORK,ovlp);

  if(compact) {

    // C = T(T1) * T(B)
    ma::product(T(T1),T(B),std::forward<MatC>(C)); 

  } else {

    // T2 = T(T1) * T(B)
    ma::product(T(T1),T(B),std::forward<Mat2>(T2)); 

    // C = conj(A) * T2
    ma::product(T(hermA),T2,std::forward<MatC>(C));

  }
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
inline void MixedDensityMatrixForWoodbury(const MatA& hermA, const MatB& B, MatC&& C, Tp* ovlp, MatD&& QQ0, integer *ref, Mat1&& TNN, Mat2&& TAB, Mat3&& TNM, IBuffer& IWORK, TBuffer& WORK, bool compact=true)
{
  // check dimensions are consistent
  int NEL = B.size(1);
  assert( hermA.size(1) == B.size(0) );
  assert( hermA.size(0) == TAB.size(0) );
  assert( B.size(1) == TAB.size(1) );
  assert( B.size(1) == TNN.size(0) );
  assert( B.size(1) == TNN.size(1) );
  assert( hermA.size(0) == QQ0.size(0) );
  assert( B.size(1) == QQ0.size(1) );
  if(compact) {
    assert( C.size(0) == TNN.size(1) );
    assert( C.size(1) == B.size(0) );
  } else {
    assert( TNM.size(1) == B.size(0) );
    assert( TNM.size(0) == TNN.size(1) );
    assert( C.size(0) == hermA.size(1) );
    assert( C.size(1) == TNM.size(1) );
  }

  using ma::T;

  // TAB = herm(A)*B
  ma::product(hermA,B,std::forward<Mat2>(TAB));  

  // TNN = TAB[ref,:] 
  for(int i=0; i<NEL; i++)
    std::copy_n(to_address(TAB[*(ref+i)].origin()),NEL,to_address(TNN[i].origin()));

  // TNN = TNN^(-1)
  ma::invert(std::forward<Mat1>(TNN),IWORK,WORK,ovlp);

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
inline void MixedDensityMatrixFromConfiguration(const MatA& hermA, const MatB& B, MatC&& C, Tp* ovlp, integer *ref, Mat1&& TNN, Mat2&& TAB, Mat3&& TNM, IBuffer& IWORK, TBuffer& WORK, bool compact=true)
{
  // check dimensions are consistent
  int NEL = B.size(1);
  assert( hermA.size(1) == B.size(0) );
  assert( hermA.size(0) == TAB.size(0) );
  assert( B.size(1) == TAB.size(1) );
  assert( B.size(1) == TNN.size(0) );
  assert( B.size(1) == TNN.size(1) );
  if(compact) {
    assert( C.size(0) == TNN.size(1) );
    assert( C.size(1) == B.size(0) );
  } else {
    assert( TNM.size(1) == B.size(0) );
    assert( TNM.size(0) == TNN.size(1) );
    assert( C.size(0) == hermA.size(1) );
    assert( C.size(1) == TNM.size(1) );
  }

  using ma::T;

  // TAB = herm(A)*B
  ma::product(hermA,B,std::forward<Mat2>(TAB));  

  // TNN = TAB[ref,:] 
  for(int i=0; i<NEL; i++)
    std::copy_n(to_address(TAB[*(ref+i)].origin()),NEL,to_address(TNN[i].origin()));

  // TNN = TNN^(-1)
  ma::invert(std::forward<Mat1>(TNN),IWORK,WORK,ovlp);

  if(compact) {

    // C = T(TNN) * T(B)
    ma::product(T(TNN),T(B),std::forward<MatC>(C)); 

  } else {

    // TNM = T(TNN) * T(B)
    ma::product(T(TNN),T(B),std::forward<Mat3>(TNM)); 

    // C = conj(A) * TNM
    ma::product(T(hermA),TNM,std::forward<MatC>(C));

  }
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
inline void MixedDensityMatrix_noHerm(const MatA& A, const MatB& B, MatC&& C, Tp* ovlp, Mat1&& T1, Mat2&& T2, IBuffer& IWORK, TBuffer& WORK, bool compact=true)
{
  // check dimensions are consistent
  assert( A.size(0) == B.size(0) );
  assert( A.size(1) == B.size(1) );
  assert( A.size(1) == T1.size(0) );
  assert( B.size(1) == T1.size(1) );
  if(compact) {
    assert( C.size(0) == T1.size(1) );
    assert( C.size(1) == B.size(0) );
  } else {
    assert( T2.size(1) == B.size(0) );
    assert( T2.size(0) == T1.size(1) );
    assert( C.size(0) == A.size(0) );
    assert( C.size(1) == T2.size(1) );
  }

  using ma::T;
  using ma::H;

  // T1 = H(A)*B
  ma::product(H(A),B,std::forward<Mat1>(T1));  

  // NOTE: Using C as temporary 
  // T1 = T1^(-1)
  ma::invert(std::forward<Mat1>(T1),IWORK,WORK,ovlp);

  if(compact) {

    // C = T(T1) * T(B)
    ma::product(T(T1),T(B),std::forward<MatC>(C)); 

  } else {

    // T2 = T1 * H(A) 
    ma::product(T1,H(A),std::forward<Mat2>(T2)); 

    // C = T( B * T2) = T(T2) * T(B)
    ma::product(T(T2),T(B),std::forward<MatC>(C));

  }
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
inline void Overlap(const MatA& hermA, const MatB& B, Tp* ovlp, Mat&& T1, IBuffer& IWORK, Buffer& WORK)
{
  // check dimensions are consistent
  assert( hermA.size(1) == B.size(0) );
  assert( hermA.size(0) == B.size(1) );
  assert( hermA.size(0) == T1.size(0) );
  assert( B.size(1) == T1.size(1) );

  using ma::T;

  // T(B)*conj(A) 
  ma::product(hermA,B,std::forward<Mat>(T1));   

  ma::determinant(std::forward<Mat>(T1),IWORK,WORK,ovlp);
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
inline void OverlapForWoodbury(const MatA& hermA, const MatB& B, Tp* ovlp, MatC&& QQ0, integer *ref, 
                             MatD && TNN, MatE && TMN, IBuffer& IWORK, TBuffer& WORK)
{
  // check dimensions are consistent
  int NEL = B.size(1);
  assert( hermA.size(1) == B.size(0) );
  assert( hermA.size(0) == TMN.size(0) );
  assert( B.size(1) == TMN.size(1) );
  assert( B.size(1) == TNN.size(0) );
  assert( B.size(1) == TNN.size(1) );
  assert( hermA.size(0) == QQ0.size(0) );
  assert( B.size(1) == QQ0.size(1) );

  using ma::T;

  // TMN = herm(A)*B 
  ma::product(hermA,B,std::forward<MatD>(TMN));

  // TNN = TMN[ref,:]
  for(int i=0; i<NEL; i++)
    std::copy_n(to_address(TMN[*(ref+i)].origin()),NEL,to_address(TNN[i].origin())); 
 
  // TNN -> inv(TNN)
  ma::invert(std::forward<MatD>(TNN),IWORK,WORK,ovlp);

  // QQ0 = TMN * inv(TNN) 
  ma::product(TMN,TNN,std::forward<MatC>(QQ0));  
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
inline void Overlap_noHerm(const MatA& A, const MatB& B, Tp* ovlp, Mat&& T1, IBuffer& IWORK, Buffer& WORK)
{
  // check dimensions are consistent
  assert( A.size(0) == B.size(0) );
  assert( A.size(1) == B.size(1) );
  assert( A.size(1) == T1.size(0) );
  assert( B.size(1) == T1.size(1) );

  using ma::T;
  using ma::H;

  // T(B)*conj(A) 
  ma::product(H(A),B,std::forward<Mat>(T1));    

  ma::determinant(std::forward<Mat>(T1),IWORK,WORK,ovlp);
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
inline void MixedDensityMatrix(const MatA& hermA, const MatB& B, MatC&& C, Tp* ovlp, Mat&& T1, Mat&& T2, IBuffer& IWORK, TBuffer& WORK, communicator& comm, bool compact=true)
{
  // check dimensions are consistent
  assert( hermA.size(1) == B.size(0) );
  assert( hermA.size(0) == B.size(1) );
  assert( hermA.size(0) == T1.size(0) );
  assert( B.size(1) == T1.size(1) );
  if(compact) {
    assert( C.size(0) == T1.size(1) );
    assert( C.size(1) == B.size(0) );
  } else {
    assert( T2.size(1) == B.size(0) );
    assert( T2.size(0) == T1.size(1) );
    assert( C.size(0) == hermA.size(1) );
    assert( C.size(1) == T2.size(1) );
  }

  using ma::T;

  int N0,Nn,sz=B.size(1);
  std::tie(N0,Nn) = FairDivideBoundary(comm.rank(),sz,comm.size());

  // T(B)*conj(A) 
  if(N0!=Nn)
    ma::product(hermA,
              B(B.extension(0),{N0,Nn}),
              T1(T1.extension(0),{N0,Nn}));  

  comm.barrier();

  // NOTE: Using C as temporary 
  // T1 = T1^(-1)
  if(comm.rank()==0)
   ma::invert(std::forward<Mat>(T1),IWORK,WORK,ovlp);
  comm.broadcast_n(ovlp,1,0);  

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
    
    sz=T2.size(1);
    std::tie(N0,Nn) = FairDivideBoundary(comm.rank(),sz,comm.size());

    // C = conj(A) * T2
    if(N0!=Nn)
      ma::product(T(hermA),
                T2(T2.extension(0),{N0,Nn}),
                C(C.extension(0),{N0,Nn}));

  }
  comm.barrier();
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
inline void Overlap(const MatA& hermA, const MatB& B, Tp* ovlp, Mat&& T1, IBuffer& IWORK, Buffer& WORK, communicator& comm)
{
  // check dimensions are consistent
  assert( hermA.size(1) == B.size(0) );
  assert( hermA.size(0) == B.size(1) );
  assert( hermA.size(0) == T1.size(0) );
  assert( B.size(1) == T1.size(1) );

  using ma::T;

  int N0,Nn,sz = B.size(1);
  std::tie(N0,Nn) = FairDivideBoundary(comm.rank(),sz,comm.size());

  // T(B)*conj(A) 
  if(N0!=Nn)
    ma::product(hermA,
              B(B.extension(0),{N0,Nn}),
              T1(T1.extension(0),{N0,Nn}));

  comm.barrier();

  if(comm.rank()==0)
   ma::determinant(std::forward<Mat>(T1),IWORK,WORK,ovlp);
  comm.broadcast_n(ovlp,1,0);  
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
inline void OverlapForWoodbury(const MatA& hermA, const MatB& B, Tp* ovlp, MatC&& QQ0, integer *ref, MatD&& TNN, MatE& TMN, IBuffer& IWORK, TBuffer& WORK, communicator& comm)
{
  // check dimensions are consistent
  int NEL = B.size(1);
  assert( hermA.size(1) == B.size(0) );
  assert( hermA.size(0) == TMN.size(0) );
  assert( B.size(1) == TMN.size(1) );
  assert( B.size(1) == TNN.size(0) );
  assert( B.size(1) == TNN.size(1) );
  assert( hermA.size(0) == QQ0.size(0) );
  assert( B.size(1) == QQ0.size(1) );

  using ma::T;

  int N0,Nn,sz = B.size(1);
  std::tie(N0,Nn) = FairDivideBoundary(comm.rank(),sz,comm.size());

  // T(B)*conj(A) 
  if(N0!=Nn)
    ma::product(hermA,
              B(B.extension(0),{N0,Nn}),
              TMN(TMN.extension(0),{N0,Nn}));
  comm.barrier();
  if(comm.rank()==0) {
    for(int i=0; i<NEL; i++)
      std::copy_n(to_address(TMN[*(ref+i)].origin()),NEL,to_address(TNN[i].origin()));
    ma::invert(std::forward<MatD>(TNN),IWORK,WORK,ovlp);
  }
  comm.broadcast_n(ovlp,1,0);

  int M0,Mn;
  sz = TMN.size(0);
  std::tie(M0,Mn) = FairDivideBoundary(comm.rank(),sz,comm.size());

  // QQ0 = TMN * inv(TNN) 
  ma::product(TMN.sliced(M0,Mn),TNN,
              QQ0({M0,Mn},QQ0.extension(1))); //.sliced(M0,Mn)); 
              //QQ0.sliced(M0,Mn)); 
  comm.barrier();
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
inline void MixedDensityMatrixForWoodbury(const MatA& hermA, const MatB& B, MatC&& C, Tp* ovlp, MatD&& QQ0, integer *ref, Mat1&& TNN, Mat2&& TAB, Mat3&& TNM, IBuffer& IWORK, TBuffer& WORK, communicator& comm, bool compact=true)
{
  // check dimensions are consistent
  int NEL = B.size(1);
  assert( hermA.size(1) == B.size(0) );
  assert( hermA.size(0) == TAB.size(0) );
  assert( B.size(1) == TAB.size(1) );
  assert( B.size(1) == TNN.size(0) );
  assert( B.size(1) == TNN.size(1) );
  assert( hermA.size(0) == QQ0.size(0) );
  assert( B.size(1) == QQ0.size(1) );
  if(compact) {
    assert( C.size(0) == TNN.size(1) );
    assert( C.size(1) == B.size(0) );
  } else {
    assert( TNM.size(1) == B.size(0) );
    assert( TNM.size(0) == TNN.size(1) );
    assert( C.size(0) == hermA.size(1) );
    assert( C.size(1) == TNM.size(1) );
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
      std::copy_n(to_address(TAB[*(ref+i)].origin())+N0,Nn-N0,to_address(TNN[i].origin())+N0);
  }

  comm.barrier();

  // TNN = TNN^(-1)
  if(comm.rank()==0)
    ma::invert(std::forward<Mat1>(TNN),IWORK,WORK,ovlp);
  comm.broadcast_n(ovlp,1,0);

  int P0,Pn;
  std::tie(P0,Pn) = FairDivideBoundary(comm.rank(),int(TAB.size(0)),comm.size());

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

    int sz=TNM.size(1);
    std::tie(N0,Nn) = FairDivideBoundary(comm.rank(),sz,comm.size());   
    comm.barrier();

    // C = conj(A) * TNM
    ma::product(T(hermA),
                TNM(TNM.extension(0),{N0,Nn}),
                C(C.extension(0), {N0,Nn})); 

  }

  comm.barrier();
}

} // namespace shm 

namespace batched
{

template< class WlkIt,
          class MatA,
          class MatC,
          class Mat,
          class TVec,
          class IBuffer,
          class TBuffer
        >
inline void MixedDensityMatrix(int nbatch, WlkIt wit, SpinTypes spin, const MatA& hermA, MatC&& C, TVec&& ovlp, Mat&& TNN3D, Mat&& TNM3D, IBuffer& IWORK, TBuffer& WORK, bool compact=true)
{
  static_assert( std::decay<TVec>::type::dimensionality == 1, " TVec::dimensionality == 1" );
  static_assert( std::decay<MatC>::type::dimensionality == 3, " MatC::dimensionality == 3" );
  static_assert( std::decay<Mat>::type::dimensionality == 3, "std::decay<Mat>::type::dimensionality == 3" );

  using ma::T;
  using ma::gemmBatched;
  using ma::getrfBatched;
  using ma::getriBatched;

  int NEL = hermA.size(0);
  int NMO = hermA.size(1);
  int wsz = ma::invert_optimal_workspace_size(TNN3D[0]);

  assert( wit->SlaterMatrix(spin).size(0) == NMO );  
  assert( wit->SlaterMatrix(spin).size(1) == NEL );  
  assert( C.size(0) == nbatch );
  assert( C.size(2) == NMO );
  if(compact)
    assert( C.size(1) == NEL );
  else
    assert( C.size(1) == NMO );
  assert( ovlp.size() == nbatch ); 
  assert( TNN3D.size(1) == NEL );
  assert( TNN3D.size(2) == NEL );
  if( not compact) {
    assert( TNM3D.size(0) == nbatch );
    assert( TNM3D.size(1) == NEL );
    assert( TNM3D.size(2) == NMO );
  }
  assert( WORK.num_elements() >= nbatch*wsz );
  assert( IWORK.num_elements() >= nbatch*(NEL+1) );

  using pointer = typename std::decay<MatC>::type::element_ptr;

  int ldw = wit->SlaterMatrix(spin).stride(0);
  int ldN = TNN3D.stride(1);
  int ldC = C.stride(1);
  std::vector<pointer> Carray;
  std::vector<pointer> Warray;
  std::vector<pointer> workArray;
  std::vector<pointer> NNarray;
  Carray.reserve(nbatch);
  Warray.reserve(nbatch);
  workArray.reserve(nbatch);
  NNarray.reserve(nbatch);
  for(int i=0; i<nbatch; i++) {
    NNarray.emplace_back(TNN3D[i].origin());
    workArray.emplace_back(WORK.origin()+i*wsz);
    Carray.emplace_back(C[i].origin());
    Warray.emplace_back((wit+i)->SlaterMatrix(spin).origin());
  }

    // T(conj(A))*B 
    for(int b=0; b<nbatch; ++b)
      ma::product(hermA,(wit+b)->SlaterMatrix(spin),TNN3D[b]);
  

    // T1 = T1^(-1)
//    for(int b=0; b<nbatch; ++b) 
//      ma::invert(TNN3D[b],IWORK,WORK,to_address(ovlp.origin())+b);
    // Invert
    getrfBatched(NEL,NNarray.data(),ldN, IWORK.origin(), IWORK.origin()+nbatch*NEL, nbatch);

    for(int i=0; i<nbatch; i++) {

      using ma::determinant_from_getrf;
      determinant_from_getrf(NEL, NNarray[i], ldN, IWORK.origin()+i*NEL,
                             to_address(ovlp.origin())+i); 

    }

    getriBatched(NEL,NNarray.data(),ldN, IWORK.origin(), workArray.data(), wsz, 
                 IWORK.origin()+nbatch*NEL, nbatch);

    if(compact) {

      // C = T(T1) * T(B)
//      for(int b=0; b<nbatch; ++b)
//        ma::product(T(TNN3D[b]),T((wit+b)->SlaterMatrix(spin)),C[b]);
      // careful with fortan ordering
      gemmBatched('T','T',NMO,NEL,NEL,ComplexType(1.0),Warray.data(),ldw,NNarray.data(),ldN,
                  ComplexType(0.0),Carray.data(),ldC,nbatch);

    } else {

      int ldM = TNM3D.stride(1);
      std::vector<pointer> NMarray;
      NMarray.reserve(nbatch);
      for(int i=0; i<nbatch; i++) 
        NMarray.emplace_back(TNM3D[i].origin());

      // T2 = T(T1) * T(B)
      //for(int b=0; b<nbatch; ++b)
      //  ma::product(T(TNN3D[b]),T((wit+b)->SlaterMatrix(spin)),TNM3D[b]);
      gemmBatched('T','T',NMO,NEL,NEL,ComplexType(1.0),Warray.data(),ldw,NNarray.data(),ldN,
                  ComplexType(0.0),NMarray.data(),ldM,nbatch);

      // C = conj(A) * T2
      for(int b=0; b<nbatch; ++b)
        ma::product(T(hermA),TNM3D[b],C[b]);

    }


}

template< class WlkIt,
          class MatA,
          class Mat,
          class TVec,
          class IBuffer
        >
inline void Overlap(int nbatch, WlkIt wit, SpinTypes spin, const MatA& hermA, TVec&& ovlp, Mat&& TNN3D, IBuffer& IWORK)
{
  static_assert( std::decay<TVec>::type::dimensionality == 1, " TVec::dimensionality == 1" );
  static_assert( std::decay<Mat>::type::dimensionality == 3, "std::decay<Mat>::type::dimensionality == 3" );

  using ma::T;
  using ma::gemmBatched;
  using ma::getrfBatched;

  int NEL = hermA.size(0);
  int NMO = hermA.size(1);

  assert( wit->SlaterMatrix(spin).size(0) == NMO );  
  assert( wit->SlaterMatrix(spin).size(1) == NEL );  
  assert( ovlp.size() == nbatch ); 
  assert( TNN3D.size(1) == NEL );
  assert( TNN3D.size(2) == NEL );
  assert( IWORK.num_elements() >= nbatch*(NEL+1) );

  using pointer = typename std::decay<Mat>::type::element_ptr;

  int ldw = wit->SlaterMatrix(spin).stride(0);
  int ldN = TNN3D.stride(1);
  std::vector<pointer> Warray;
  std::vector<pointer> NNarray;
  Warray.reserve(nbatch);
  NNarray.reserve(nbatch);
  for(int i=0; i<nbatch; i++) {
    NNarray.emplace_back(TNN3D[i].origin());
    Warray.emplace_back((wit+i)->SlaterMatrix(spin).origin());
  }

    // T(conj(A))*B 
    for(int b=0; b<nbatch; ++b)
      ma::product(hermA,(wit+b)->SlaterMatrix(spin),TNN3D[b]);
  

    // T1 = T1^(-1)
//    for(int b=0; b<nbatch; ++b) 
//      ma::invert(TNN3D[b],IWORK,WORK,to_address(ovlp.origin())+b);
    // Invert
    getrfBatched(NEL,NNarray.data(),ldN, IWORK.origin(), IWORK.origin()+nbatch*NEL, nbatch);

    for(int i=0; i<nbatch; i++) {

      using ma::determinant_from_getrf;
      determinant_from_getrf(NEL, NNarray[i], ldN, IWORK.origin()+i*NEL,
                             to_address(ovlp.origin())+i); 

    }

}


} // namespace batched

} // namespace SlaterDeterminantOperations

} // namespace afqmc

} // namespace qmcplusplus 

#endif

