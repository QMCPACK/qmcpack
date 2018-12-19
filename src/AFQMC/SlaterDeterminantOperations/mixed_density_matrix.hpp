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
          class IBuffer
        >
inline Tp Overlap(const MatA& hermA, const MatB& B, Mat&& T1, IBuffer& IWORK)
{
  // check dimensions are consistent
  assert( hermA.shape()[1] == B.shape()[0] );
  assert( hermA.shape()[0] == B.shape()[1] );
  assert( hermA.shape()[0] == T1.shape()[0] );
  assert( B.shape()[1] == T1.shape()[1] );

  using ma::T;

  // T(B)*conj(A) 
  ma::product(hermA,B,std::forward<Mat>(T1));   

  return static_cast<Tp>(ma::determinant(std::forward<Mat>(T1),IWORK));
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
          class IBuffer
        >
inline Tp Overlap_noHerm(const MatA& A, const MatB& B, Mat&& T1, IBuffer& IWORK)
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

  return static_cast<Tp>(ma::determinant(std::forward<Mat>(T1),IWORK));
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
//  using boost::indices;
//  using range_t = boost::multi_array_types::index_range;

  int N0,Nn,sz=B.shape()[1];
  std::tie(N0,Nn) = FairDivideBoundary(comm.rank(),sz,comm.size());

  // T(B)*conj(A) 
  if(N0!=Nn)
    ma::product(hermA,
              B[indices[range_t()][range_t(N0,Nn)]],
              T1[indices[range_t()][range_t(N0,Nn)]]);  

  comm.barrier();

  // NOTE: Using C as temporary 
  // T1 = T1^(-1)
  Tp ovlp=Tp(0.);
  if(comm.rank()==0)
   ovlp = static_cast<Tp>(ma::invert(std::forward<Mat>(T1),IWORK,WORK));
  comm.broadcast_value(ovlp);  

  if(compact) {

    // C = T(T1) * T(B)
    //ma::product(T1[indices[range_t(N0,Nn)][range_t()]],
    //            T(B),
    //            C[indices[range_t(N0,Nn)][range_t()]]); 
    if(N0!=Nn)
      ma::product(T(T1[indices[range_t()][range_t(N0,Nn)]]),
                T(B),
                C[indices[range_t(N0,Nn)][range_t()]]); 

  } else {

    // T2 = T(T1) * T(B)
    //ma::product(T1[indices[range_t(N0,Nn)][range_t()]],
    //            T(B),
    //            T2[indices[range_t(N0,Nn)][range_t()]]); 
    if(N0!=Nn)
      ma::product(T(T1[indices[range_t()][range_t(N0,Nn)]]),
                T(B),
                T2[indices[range_t(N0,Nn)][range_t()]]); 

    comm.barrier();
    
    sz=T2.shape()[1];
    std::tie(N0,Nn) = FairDivideBoundary(comm.rank(),sz,comm.size());

    // C = conj(A) * T2
    if(N0!=Nn)
      ma::product(T(hermA),
                T2[indices[range_t()][range_t(N0,Nn)]],
                C[indices[range_t()][range_t(N0,Nn)]]);

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
          class communicator
        >
inline Tp Overlap(const MatA& hermA, const MatB& B, Mat&& T1, IBuffer& IWORK, communicator& comm)
{
  // check dimensions are consistent
  assert( hermA.shape()[1] == B.shape()[0] );
  assert( hermA.shape()[0] == B.shape()[1] );
  assert( hermA.shape()[0] == T1.shape()[0] );
  assert( B.shape()[1] == T1.shape()[1] );

  using ma::T;
//  using boost::indices;
//  using range_t = boost::multi_array_types::index_range;

  int N0,Nn,sz = B.shape()[1];
  std::tie(N0,Nn) = FairDivideBoundary(comm.rank(),sz,comm.size());

  // T(B)*conj(A) 
  if(N0!=Nn)
    ma::product(hermA,
              B[indices[range_t()][range_t(N0,Nn)]],
              T1[indices[range_t()][range_t(N0,Nn)]]);

  comm.barrier();

  Tp ovlp=Tp(0.);
  if(comm.rank()==0)
   ovlp = static_cast<Tp>(ma::determinant(std::forward<Mat>(T1),IWORK));
  comm.broadcast_value(ovlp);  

  return ovlp; 
}

} // namespace shm 

} // namespace SlaterDeterminantOperations

} // namespace afqmc

} // namespace qmcplusplus 

#endif

