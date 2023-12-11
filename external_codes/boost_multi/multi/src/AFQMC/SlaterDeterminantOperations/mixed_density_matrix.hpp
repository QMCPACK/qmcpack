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


#ifndef AFQMC_DENSITY_MATRIX_HPP
#define AFQMC_DENSITY_MATRIX_HPP

#include "AFQMC/config.h"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Numerics/tensor_operations.hpp"
#include "Utilities/FairDivide.h"

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
template<class Tp, class MatA, class MatB, class MatC, class Mat1, class Mat2, class IBuffer, class TBuffer>
Tp MixedDensityMatrix(const MatA& hermA,
                      const MatB& B,
                      MatC&& C,
                      Tp LogOverlapFactor,
                      Mat1&& T1,
                      Mat2&& T2,
                      IBuffer& IWORK,
                      TBuffer& WORK,
                      bool compact = true,
                      bool herm    = true)
{
  // check dimensions are consistent
  int NMO = (herm ? std::get<1>(hermA.sizes()) : std::get<0>(hermA.sizes()));
  int NEL = (herm ? std::get<0>(hermA.sizes()) : std::get<1>(hermA.sizes()));
  assert(NMO == std::get<0>(B.sizes()));
  assert(NEL == std::get<1>(B.sizes()));
  assert(NEL == T1.size());
  assert(std::get<1>(B.sizes()) == std::get<1>(T1.sizes()));
  if (compact)
  {
    assert(std::get<0>(C.sizes()) == std::get<1>(T1.sizes()));
    assert(std::get<1>(C.sizes()) == std::get<0>(B.sizes()));
  }
  else
  {
    assert(std::get<1>(T2.sizes()) == B.size());
    assert(T2.size() == std::get<1>(T1.sizes()));
    assert(C.size() == NMO);
    assert(std::get<1>(C.sizes()) == std::get<1>(T2.sizes()));
  }

  using ma::H;
  using ma::T;

  // H(A) * B
  if (herm)
    ma::product(hermA, B, std::forward<Mat1>(T1));
  else
    ma::product(H(hermA), B, std::forward<Mat1>(T1));

  // NOTE: Using C as temporary
  // T1 = T1^(-1)
  Tp ovlp = ma::invert(std::forward<Mat1>(T1), IWORK, WORK, LogOverlapFactor);
  if (std::abs(ovlp) == Tp(0.0))
  {
    using element = typename std::decay<MatC>::type::element;
    ma::fill(C, element(0.0));
    return ovlp;
  }

  if (compact)
  {
    // C = T(T1) * T(B)
    ma::product(T(T1), T(B), std::forward<MatC>(C));
  }
  else
  {
    if (herm)
    {
      // T2 = T(T1) * T(B)
      ma::product(T(T1), T(B), std::forward<Mat2>(T2));

      // C = conj(A) * T2
      ma::product(T(hermA), T2, std::forward<MatC>(C));
    }
    else
    {
      // T2 = T1 * H(A)
      ma::product(T1, H(hermA), std::forward<Mat2>(T2));

      // C = T( B * T2) = T(T2) * T(B)
      ma::product(T(T2), T(B), std::forward<MatC>(C));
    }
  }
  return ovlp;
}

template<class Tp,
         class integer,
         class MatA,
         class MatB,
         class MatC,
         class MatD,
         class Mat1,
         class Mat2,
         class Mat3,
         class IBuffer,
         class TBuffer>
Tp MixedDensityMatrixForWoodbury(const MatA& hermA,
                                 const MatB& B,
                                 MatC&& C,
                                 Tp LogOverlapFactor,
                                 MatD&& QQ0,
                                 integer* ref,
                                 Mat1&& TNN,
                                 Mat2&& TAB,
                                 Mat3&& TNM,
                                 IBuffer& IWORK,
                                 TBuffer& WORK,
                                 bool compact = true)
{
  // check dimensions are consistent
  int NEL = std::get<1>(B.sizes());
  assert(std::get<1>(hermA.sizes()) == std::get<0>(B.sizes()));
  assert(std::get<0>(hermA.sizes()) == std::get<0>(TAB.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<1>(TAB.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<0>(TNN.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<1>(TNN.sizes()));
  assert(std::get<0>(hermA.sizes()) == std::get<0>(QQ0.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<1>(QQ0.sizes()));
  if (compact)
  {
    assert(std::get<0>(C.sizes()) == std::get<1>(TNN.sizes()));
    assert(std::get<1>(C.sizes()) == std::get<0>(B.sizes()));
  }
  else
  {
    assert(std::get<1>(TNM.sizes()) == std::get<0>(B.sizes()));
    assert(std::get<0>(TNM.sizes()) == std::get<1>(TNN.sizes()));
    assert(std::get<0>(C.sizes()) == std::get<1>(hermA.sizes()));
    assert(std::get<1>(C.sizes()) == std::get<1>(TNM.sizes()));
  }

  using ma::T;

  // TAB = herm(A)*B
  ma::product(hermA, B, std::forward<Mat2>(TAB));

  // TNN = TAB[ref,:]
  for (int i = 0; i < NEL; i++)
    std::copy_n(to_address(TAB[*(ref + i)].origin()), NEL, to_address(TNN[i].origin()));

  // TNN = TNN^(-1)
  Tp ovlp = ma::invert(std::forward<Mat1>(TNN), IWORK, WORK, LogOverlapFactor);
  if (std::abs(ovlp) == Tp(0.0))
  {
    using element = typename std::decay<MatC>::type::element;
    ma::fill(C, element(0.0));
    return ovlp;
  }

  // QQ0 = TAB * inv(TNN)
  ma::product(TAB, TNN, std::forward<MatD>(QQ0));

  if (compact)
  {
    // C = T(TNN) * T(B)
    ma::product(T(TNN), T(B), std::forward<MatC>(C));
  }
  else
  {
    // TNM = T(TNN) * T(B)
    ma::product(T(TNN), T(B), std::forward<Mat3>(TNM));

    // C = conj(A) * TNM
    ma::product(T(hermA), TNM, std::forward<MatC>(C));
  }
  return ovlp;
}

template<class Tp,
         class integer,
         class MatA,
         class MatB,
         class MatC,
         class Mat1,
         class Mat2,
         class Mat3,
         class IBuffer,
         class TBuffer>
Tp MixedDensityMatrixFromConfiguration(const MatA& hermA,
                                       const MatB& B,
                                       MatC&& C,
                                       Tp LogOverlapFactor,
                                       integer* ref,
                                       Mat1&& TNN,
                                       Mat2&& TAB,
                                       Mat3&& TNM,
                                       IBuffer& IWORK,
                                       TBuffer& WORK,
                                       bool compact = true)
{
  // check dimensions are consistent
  int NEL = std::get<1>(B.sizes());
  assert(std::get<1>(hermA.sizes()) == std::get<0>(B.sizes()));
  assert(std::get<0>(hermA.sizes()) == std::get<0>(TAB.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<1>(TAB.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<0>(TNN.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<1>(TNN.sizes()));
  if (compact)
  {
    assert(std::get<0>(C.sizes()) == std::get<1>(TNN.sizes()));
    assert(std::get<1>(C.sizes()) == std::get<0>(B.sizes()));
  }
  else
  {
    assert(std::get<1>(TNM.sizes()) == std::get<0>(B.sizes()));
    assert(std::get<0>(TNM.sizes()) == std::get<1>(TNN.sizes()));
    assert(std::get<0>(C.sizes()) == std::get<1>(hermA.sizes()));
    assert(std::get<1>(C.sizes()) == std::get<1>(TNM.sizes()));
  }

  using ma::T;

  // TAB = herm(A)*B
  ma::product(hermA, B, std::forward<Mat2>(TAB));

  // TNN = TAB[ref,:]
  for (int i = 0; i < NEL; i++)
    std::copy_n(to_address(TAB[*(ref + i)].origin()), NEL, to_address(TNN[i].origin()));

  // TNN = TNN^(-1)
  Tp ovlp = ma::invert(std::forward<Mat1>(TNN), IWORK, WORK, LogOverlapFactor);
  if (std::abs(ovlp) == Tp(0.0))
  {
    using element = typename std::decay<MatC>::type::element;
    ma::fill(C, element(0.0));
    return ovlp;
  }

  if (compact)
  {
    // C = T(TNN) * T(B)
    ma::product(T(TNN), T(B), std::forward<MatC>(C));
  }
  else
  {
    // TNM = T(TNN) * T(B)
    ma::product(T(TNN), T(B), std::forward<Mat3>(TNM));

    // C = conj(A) * TNM
    ma::product(T(hermA), TNM, std::forward<MatC>(C));
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
template<class Tp, class MatA, class MatB, class MatC, class Mat1, class Mat2, class IBuffer, class TBuffer>
Tp MixedDensityMatrix_noHerm(const MatA& A,
                             const MatB& B,
                             MatC&& C,
                             Tp LogOverlapFactor,
                             Mat1&& T1,
                             Mat2&& T2,
                             IBuffer& IWORK,
                             TBuffer& WORK,
                             bool compact = true)
{
  // check dimensions are consistent
  assert(std::get<0>(A.sizes()) == std::get<0>(B.sizes()));
  assert(std::get<1>(A.sizes()) == std::get<1>(B.sizes()));
  assert(std::get<1>(A.sizes()) == std::get<0>(T1.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<1>(T1.sizes()));
  if (compact)
  {
    assert(std::get<0>(C.sizes()) == std::get<1>(T1.sizes()));
    assert(std::get<1>(C.sizes()) == std::get<0>(B.sizes()));
  }
  else
  {
    assert(std::get<1>(T2.sizes()) == std::get<0>(B.sizes()));
    assert(std::get<0>(T2.sizes()) == std::get<1>(T1.sizes()));
    assert(std::get<0>(C.sizes()) == std::get<0>(A.sizes()));
    assert(std::get<1>(C.sizes()) == std::get<1>(T2.sizes()));
  }

  using ma::H;
  using ma::T;

  // T1 = H(A)*B
  ma::product(H(A), B, std::forward<Mat1>(T1));

  // NOTE: Using C as temporary
  // T1 = T1^(-1)
  Tp ovlp = ma::invert(std::forward<Mat1>(T1), IWORK, WORK, LogOverlapFactor);
  if (std::abs(ovlp) == Tp(0.0))
  {
    using element = typename std::decay<MatC>::type::element;
    ma::fill(C, element(0.0));
    return ovlp;
  }

  if (compact)
  {
    // C = T(T1) * T(B)
    ma::product(T(T1), T(B), std::forward<MatC>(C));
  }
  else
  {
    // T2 = T1 * H(A)
    ma::product(T1, H(A), std::forward<Mat2>(T2));

    // C = T( B * T2) = T(T2) * T(B)
    ma::product(T(T2), T(B), std::forward<MatC>(C));
  }
  return ovlp;
}

template<class Tp, class MatA, class MatB, class MatC, class Mat1, class RVec, class IBuffer, class TBuffer>
Tp MixedDensityMatrix_noHerm_wSVD(const MatA& A,
                                  const MatB& B,
                                  MatC&& C,
                                  Tp LogOverlapFactor,
                                  RVec&& S,
                                  Mat1&& U,
                                  Mat1&& VT,
                                  Mat1&& BV,
                                  Mat1&& UA,
                                  IBuffer& IWORK,
                                  TBuffer& WORK,
                                  bool compact = true)
{
  // check dimensions are consistent
  assert(std::get<0>(A.sizes()) == std::get<0>(B.sizes()));
  assert(std::get<1>(A.sizes()) == std::get<1>(B.sizes()));
  assert(std::get<1>(A.sizes()) == std::get<0>(U.sizes())); // [U] = [NxN]
  assert(std::get<1>(A.sizes()) == std::get<1>(U.sizes()));
  assert(std::get<1>(A.sizes()) == std::get<0>(VT.sizes())); // [V] = [NxN]
  assert(std::get<1>(A.sizes()) == std::get<1>(VT.sizes()));
  assert(std::get<1>(A.sizes()) <= (6 * S.size() + 1)); // [S] = [N+1]
  assert(std::get<1>(A.sizes()) == std::get<0>(UA.sizes()));          // [UA] = [NxM]
  assert(std::get<0>(A.sizes()) == std::get<1>(UA.sizes()));
  if (compact)
  {
    assert(std::get<0>(C.sizes()) == std::get<1>(B.sizes()));
    assert(std::get<1>(C.sizes()) == std::get<0>(B.sizes()));
  }
  else
  {
    assert( std::get<0>(A.sizes()) == std::get<0>(BV.sizes()) ); // [BV] = [MxN]
    assert( std::get<1>(A.sizes()) == std::get<1>(BV.sizes()) );
    assert( std::get<0>(C.sizes()) == std::get<0>(A.sizes()) );
    assert( std::get<1>(C.sizes()) == std::get<0>(A.sizes()) );
  }

  using ma::determinant_from_geqrf;
  using ma::H;
  using ma::real;
  using ma::T;
  using ma::term_by_term_matrix_vector;

  int N(U.size());

  // T1 = H(A)*B
  ma::product(H(A), B, U);

  // keep a copy of U in UA temporarily
  using std::copy_n;
  copy_n(U.origin(), U.num_elements(), UA.origin());

  // get determinant, since you can't get the phase trivially from SVD
  Tp ovlp = ma::determinant(U, IWORK, WORK, LogOverlapFactor);
  //  ma::geqrf(U,VT[0],WORK);
  //  determinant_from_geqrf(N,U.origin(),U.stride(0),VT[1].origin(),LogOverlapFactor,ovlp);
  //  if you want the correct phase of the determinant
  //  ma::gqr(U,S.sliced(0,N),WORK);
  //  ComplexType ovQ = ma::determinant(U,IWORK,WORK,0.0);
  //  *ovlp *= ovQ;

  // restore U
  copy_n(UA.origin(), U.num_elements(), U.origin());

  // H(A)*B = AtB = U * S * VT
  // inv(H(A)*B) = inv(AtB) = H(VT) * inv(S) * H(U)
  ma::gesvd('O', 'A', U, S.sliced(0, N), U, VT, WORK, S.sliced(N, 6 * N));

  // testing
  boost::multi::array<double, 1> Sh(S.sliced(0, N));
  double ov_(0.0);
  for (int i = 0; i < N; i++)
    ov_ += std::log(Sh[i]);
  ov_ = std::exp(ov_ - real(LogOverlapFactor));
  //  std::cout<<" SVD: " <<ov0 <<" " <<ov_ <<" " <<ov0/ov_ <<" " <<Sh[0] <<" " <<Sh[N-1] <<" " <<Sh[0]/Sh[N-1] <<std::endl;

  // mod Sh
  // S = Sh;

  if (compact)
  {
    // G = T( B * inv(AtB) )
    //   = T( B * H(VT) * inv(S) * H(U) )
    //   = T(H(VT) * inv(S) * H(U)) * T(B)


    // VT = VT * inv(S), which works since S is diagonal and real
    term_by_term_matrix_vector(ma::TOp_DIV, 0, std::get<0>(VT.sizes()), std::get<1>(VT.sizes()), ma::pointer_dispatch(VT.origin()), VT.stride(0),
                               ma::pointer_dispatch(S.origin()), 1);

    // BV = H(VT) * H(U)
    ma::product(H(VT), H(U), BV.sliced(0, N));

    // G = T(BV) * T(B)
    product(T(BV.sliced(0, N)), T(B), C);
  }
  else
  {
    // G = T( B * inv(AtB) * H(A) )
    //   = T( B * H(VT) * inv(S) * H(U) * H(A) )
    //   = T( BV * inv(S) * UA )
    //   = T(UA) * T(BV * inv(S))

    // BV = B * H(VT)
    ma::product(B, H(VT), BV);

    // BV = BV * inv(S), which works since S is diagonal and real
    term_by_term_matrix_vector(ma::TOp_DIV, 1, std::get<0>(BV.sizes()), std::get<1>(BV.sizes()), ma::pointer_dispatch(BV.origin()), BV.stride(0),
                               ma::pointer_dispatch(S.origin()), 1);

    // UA = H(U) * H(A)
    ma::product(H(U), H(A), UA);

    // G = T(UA) * T(BV * inv(S))
    product(T(UA), T(BV), C);
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
template<class Tp, class MatA, class MatB, class Mat, class Buffer, class IBuffer>
Tp Overlap(const MatA& hermA,
           const MatB& B,
           Tp LogOverlapFactor,
           Mat&& T1,
           IBuffer&& IWORK,
           Buffer&& WORK,
           bool herm = true)
{
  int NMO = (herm ? std::get<1>(hermA.sizes()) : std::get<0>(hermA.sizes()));
  int NEL = (herm ? std::get<0>(hermA.sizes()) : std::get<1>(hermA.sizes()));
  // check dimensions are consistent
  assert(NMO == std::get<0>(B.sizes()));
  assert(NEL == std::get<1>(B.sizes()));
  assert(NEL == std::get<0>(T1.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<1>(T1.sizes()));

  using ma::H;
  using ma::T;

  // T(B)*conj(A)
  if (herm)
    ma::product(hermA, B, std::forward<Mat>(T1));
  else
    ma::product(H(hermA), B, std::forward<Mat>(T1));

  return ma::determinant(std::forward<Mat>(T1), IWORK, WORK, LogOverlapFactor);
}

template<class Tp,
         class integer,
         class MatA,
         class MatB,
         class MatC,
         class MatD,
         class MatE,
         class IBuffer,
         class TBuffer>
Tp OverlapForWoodbury(const MatA& hermA,
                      const MatB& B,
                      Tp LogOverlapFactor,
                      MatC&& QQ0,
                      integer* ref,
                      MatD&& TNN,
                      MatE&& TMN,
                      IBuffer& IWORK,
                      TBuffer& WORK)
{
  // check dimensions are consistent
  int NEL = std::get<1>(B.sizes());
  assert(std::get<1>(hermA.sizes()) == std::get<0>(B.sizes()));
  assert(std::get<0>(hermA.sizes()) == std::get<0>(TMN.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<1>(TMN.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<0>(TNN.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<1>(TNN.sizes()));
  assert(std::get<0>(hermA.sizes()) == std::get<0>(QQ0.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<1>(QQ0.sizes()));

  using ma::T;

  // TMN = herm(A)*B
  ma::product(hermA, B, std::forward<MatD>(TMN));

  // TNN = TMN[ref,:]
  for (int i = 0; i < NEL; i++)
    std::copy_n(to_address(TMN[*(ref + i)].origin()), NEL, to_address(TNN[i].origin()));

  // TNN -> inv(TNN)
  Tp ovlp = ma::invert(std::forward<MatD>(TNN), IWORK, WORK, LogOverlapFactor);
  if (std::abs(ovlp) == Tp(0.0))
  {
    using element = typename std::decay<MatC>::type::element;
    ma::fill(QQ0, element(0.0));
    return ovlp;
  }

  // QQ0 = TMN * inv(TNN)
  ma::product(TMN, TNN, std::forward<MatC>(QQ0));

  return ovlp;
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
template<class Tp, class MatA, class MatB, class Mat, class Buffer, class IBuffer>
Tp Overlap_noHerm(const MatA& A, const MatB& B, Tp LogOverlapFactor, Mat&& T1, IBuffer& IWORK, Buffer& WORK)
{
  // check dimensions are consistent
  assert(std::get<0>(A.sizes()) == std::get<0>(B.sizes()));
  assert(std::get<1>(A.sizes()) == std::get<1>(B.sizes()));
  assert(std::get<1>(A.sizes()) == std::get<0>(T1.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<1>(T1.sizes()));

  using ma::H;
  using ma::T;

  // T(B)*conj(A)
  ma::product(H(A), B, std::forward<Mat>(T1));

  return ma::determinant(std::forward<Mat>(T1), IWORK, WORK, LogOverlapFactor);
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
template<class Tp, class MatA, class MatB, class MatC, class Mat, class IBuffer, class TBuffer, class communicator>
Tp MixedDensityMatrix(const MatA& hermA,
                      const MatB& B,
                      MatC&& C,
                      Tp LogOverlapFactor,
                      Mat&& T1,
                      Mat&& T2,
                      IBuffer& IWORK,
                      TBuffer& WORK,
                      communicator& comm,
                      bool compact = true,
                      bool herm    = true)
{
  int NMO = (herm ? std::get<1>(hermA.sizes()) : std::get<0>(hermA.sizes()));
  int NEL = (herm ? std::get<0>(hermA.sizes()) : std::get<1>(hermA.sizes()));
  // check dimensions are consistent
  assert(NMO == std::get<0>(B.sizes()));
  assert(NEL == std::get<1>(B.sizes()));
  assert(NEL == std::get<0>(T1.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<1>(T1.sizes()));
  if (compact)
  {
    assert(std::get<0>(C.sizes()) == std::get<1>(T1.sizes()));
    assert(std::get<1>(C.sizes()) == std::get<0>(B.sizes()));
  }
  else
  {
    assert(std::get<1>(T2.sizes()) == std::get<0>(B.sizes()));
    assert(std::get<0>(T2.sizes()) == std::get<1>(T1.sizes()));
    assert(std::get<0>(C.sizes()) == NMO);
    assert(std::get<1>(C.sizes()) == std::get<1>(T2.sizes()));
  }

  using ma::H;
  using ma::T;

  int N0, Nn, sz = std::get<1>(B.sizes());
  std::tie(N0, Nn) = FairDivideBoundary(comm.rank(), sz, comm.size());

  // T(B)*conj(A)
  if (N0 != Nn)
  {
    if (herm)
      ma::product(hermA, B(B.extension(0), {N0, Nn}), T1(T1.extension(0), {N0, Nn}));
    else
      ma::product(H(hermA), B(B.extension(0), {N0, Nn}), T1(T1.extension(0), {N0, Nn}));
  }

  comm.barrier();

  // NOTE: Using C as temporary
  // T1 = T1^(-1)
  Tp ovlp;
  if (comm.rank() == 0)
  {
    ovlp = ma::invert(std::forward<Mat>(T1), IWORK, WORK, LogOverlapFactor);
    if (std::abs(ovlp) == Tp(0.0))
    {
      using element = typename std::decay<MatC>::type::element;
      ma::fill(C, element(0.0));
    }
  }
  comm.broadcast_n(&ovlp, 1, 0);
  if (std::abs(ovlp) == Tp(0.0))
    return ovlp;

  if (compact)
  {
    // C = T(T1) * T(B)
    //ma::product(T1.sliced(N0,Nn),
    //            T(B),
    //            C.sliced(N0,Nn));
    if (N0 != Nn)
      ma::product(T(T1(T1.extension(0), {N0, Nn})), T(B), C.sliced(N0, Nn));
  }
  else
  {
    if (herm)
    {
      // T2 = T(T1) * T(B)
      //ma::product(T1.sliced(N0,Nn),
      //            T(B),
      //            T2.sliced(N0,Nn));

      if (N0 != Nn)
        ma::product(ComplexType(1.0), T(T1(T1.extension(0), {N0, Nn})), T(B), ComplexType(0.0), T2.sliced(N0, Nn));

      comm.barrier();

      sz               = std::get<1>(T2.sizes());
      std::tie(N0, Nn) = FairDivideBoundary(comm.rank(), sz, comm.size());

      // C = conj(A) * T2
      if (N0 != Nn)
        ma::product(T(hermA), T2(T2.extension(0), {N0, Nn}), C(C.extension(0), {N0, Nn}));
    }
    else
    {
      // T2 = T(T1) * T(B)
      // T2 = T1 * H(A)
      if (N0 != Nn)
        ma::product(T1.sliced(N0, Nn), H(hermA), T2.sliced(N0, Nn));

      comm.barrier();

      sz               = std::get<1>(T2.sizes());
      std::tie(N0, Nn) = FairDivideBoundary(comm.rank(), sz, comm.size());

      // C = T( B * T2) = T(T2) * T(B)
      if (N0 != Nn)
        ma::product(T(T2(T2.extension(0), {N0, Nn})), T(B), C.sliced(N0, Nn));
    }
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
template<class Tp, class MatA, class MatB, class Mat, class IBuffer, class Buffer, class communicator>
Tp Overlap(const MatA& hermA,
           const MatB& B,
           Tp LogOverlapFactor,
           Mat&& T1,
           IBuffer&& IWORK,
           Buffer&& WORK,
           communicator& comm,
           bool herm = true)
{
  int NMO = (herm ? std::get<1>(hermA.sizes()) : std::get<0>(hermA.sizes()));
  int NEL = (herm ? std::get<0>(hermA.sizes()) : std::get<1>(hermA.sizes()));
  // check dimensions are consistent
  assert(NMO == std::get<0>(B.sizes()));
  assert(NEL == std::get<1>(B.sizes()));
  assert(NEL == std::get<0>(T1.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<1>(T1.sizes()));

  using ma::H;
  using ma::T;

  int N0, Nn, sz = std::get<1>(B.sizes());
  std::tie(N0, Nn) = FairDivideBoundary(comm.rank(), sz, comm.size());

  // T(B)*conj(A)
  if (N0 != Nn)
  {
    if (herm)
      ma::product(hermA, B(B.extension(0), {N0, Nn}), T1(T1.extension(0), {N0, Nn}));
    else
      ma::product(H(hermA), B(B.extension(0), {N0, Nn}), T1(T1.extension(0), {N0, Nn}));
  }

  comm.barrier();

  Tp ovlp;
  if (comm.rank() == 0)
    ovlp = ma::determinant(std::forward<Mat>(T1), IWORK, WORK, LogOverlapFactor);
  comm.broadcast_n(&ovlp, 1, 0);
  return ovlp;
}

// Serial Implementation
template<class Tp,
         class integer,
         class MatA,
         class MatB,
         class MatC,
         class MatD,
         class MatE,
         class IBuffer,
         class TBuffer,
         class communicator>
Tp OverlapForWoodbury(const MatA& hermA,
                      const MatB& B,
                      Tp LogOverlapFactor,
                      MatC&& QQ0,
                      integer* ref,
                      MatD&& TNN,
                      MatE& TMN,
                      IBuffer& IWORK,
                      TBuffer& WORK,
                      communicator& comm)
{
  // check dimensions are consistent
  int NEL = std::get<1>(B.sizes());
  assert(std::get<1>(hermA.sizes()) == std::get<0>(B.sizes()));
  assert(std::get<0>(hermA.sizes()) == std::get<0>(TMN.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<1>(TMN.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<0>(TNN.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<1>(TNN.sizes()));
  assert(std::get<0>(hermA.sizes()) == std::get<0>(QQ0.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<1>(QQ0.sizes()));

  using ma::T;

  int N0, Nn, sz = std::get<1>(B.sizes());
  std::tie(N0, Nn) = FairDivideBoundary(comm.rank(), sz, comm.size());

  Tp ovlp;
  // T(B)*conj(A)
  if (N0 != Nn)
    ma::product(hermA, B(B.extension(0), {N0, Nn}), TMN(TMN.extension(0), {N0, Nn}));
  comm.barrier();
  if (comm.rank() == 0)
  {
    for (int i = 0; i < NEL; i++)
      std::copy_n(to_address(TMN[*(ref + i)].origin()), NEL, to_address(TNN[i].origin()));
    ovlp = ma::invert(std::forward<MatD>(TNN), IWORK, WORK, LogOverlapFactor);
  }
  comm.broadcast_n(&ovlp, 1, 0);

  int M0, Mn;
  sz               = TMN.size();
  std::tie(M0, Mn) = FairDivideBoundary(comm.rank(), sz, comm.size());

  // QQ0 = TMN * inv(TNN)
  ma::product(TMN.sliced(M0, Mn), TNN, QQ0({M0, Mn}, QQ0.extension(1))); //.sliced(M0,Mn));
                                                                         //QQ0.sliced(M0,Mn));
  comm.barrier();
  return ovlp;
}

template<class Tp,
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
         class communicator>
Tp MixedDensityMatrixForWoodbury(const MatA& hermA,
                                 const MatB& B,
                                 MatC&& C,
                                 Tp LogOverlapFactor,
                                 MatD&& QQ0,
                                 integer* ref,
                                 Mat1&& TNN,
                                 Mat2&& TAB,
                                 Mat3&& TNM,
                                 IBuffer& IWORK,
                                 TBuffer& WORK,
                                 communicator& comm,
                                 bool compact = true)
{
  // check dimensions are consistent
  int NEL = std::get<1>(B.sizes());
  assert(std::get<1>(hermA.sizes()) == std::get<0>(B.sizes()));
  assert(std::get<0>(hermA.sizes()) == std::get<0>(TAB.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<1>(TAB.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<0>(TNN.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<1>(TNN.sizes()));
  assert(std::get<0>(hermA.sizes()) == std::get<0>(QQ0.sizes()));
  assert(std::get<1>(B.sizes()) == std::get<1>(QQ0.sizes()));
  if (compact)
  {
    assert(std::get<0>(C.sizes()) == std::get<1>(TNN.sizes()));
    assert(std::get<1>(C.sizes()) == std::get<0>(B.sizes()));
  }
  else
  {
    assert(std::get<1>(TNM.sizes()) == std::get<0>(B.sizes()));
    assert(std::get<0>(TNM.sizes()) == std::get<1>(TNN.sizes()));
    assert(std::get<0>(C.sizes()) == std::get<1>(hermA.sizes()));
    assert(std::get<1>(C.sizes()) == std::get<1>(TNM.sizes()));
  }

  using ma::T;

  int N0, Nn;
  std::tie(N0, Nn) = FairDivideBoundary(comm.rank(), NEL, comm.size());

  // TAB = herm(A)*B
  if (N0 != Nn)
  {
    ma::product(hermA, B(B.extension(0), {N0, Nn}), TAB(TAB.extension(0), {N0, Nn}));

    // TNN = TAB[ref,:]
    for (int i = 0; i < NEL; i++)
      std::copy_n(to_address(TAB[*(ref + i)].origin()) + N0, Nn - N0, to_address(TNN[i].origin()) + N0);
  }

  comm.barrier();

  // TNN = TNN^(-1)
  Tp ovlp;
  if (comm.rank() == 0)
    ovlp = ma::invert(std::forward<Mat1>(TNN), IWORK, WORK, LogOverlapFactor);
  comm.broadcast_n(&ovlp, 1, 0);

  int P0, Pn;
  std::tie(P0, Pn) = FairDivideBoundary(comm.rank(), int(std::get<0>(TAB.sizes())), comm.size());

  // QQ0 = TAB * inv(TNN)
  if (P0 != Pn)
    ma::product(TAB.sliced(P0, Pn), TNN, QQ0({P0, Pn}, QQ0.extension(1)));
  //QQ0.sliced(P0,Pn));
  if (compact)
  {
    // C = T(TNN) * T(B)
    if (N0 != Nn)
      ma::product(T(TNN(TNN.extension(0), {N0, Nn})), T(B), C({N0, Nn}, C.extension(1)));
  }
  else
  {
    // TNM = T(TNN) * T(B)
    if (N0 != Nn)
      ma::product(T(TNN(TNN.extension(0), {N0, Nn})), T(B), TNM.sliced(N0, Nn));

    int sz           = std::get<1>(TNM.sizes());
    std::tie(N0, Nn) = FairDivideBoundary(comm.rank(), sz, comm.size());
    comm.barrier();

    // C = conj(A) * TNM
    ma::product(T(hermA), TNM(TNM.extension(0), {N0, Nn}), C(C.extension(0), {N0, Nn}));
  }

  comm.barrier();
  return ovlp;
}

} // namespace shm

namespace batched
{
template<class MatA, class MatB, class MatC, class Mat, class TVec, class IBuffer, class Tp>
void MixedDensityMatrix(std::vector<MatA>& hermA,
                        std::vector<MatB>& Bi,
                        MatC&& C,
                        Tp LogOverlapFactor,
                        TVec&& ovlp,
                        Mat&& TNN3D,
                        Mat&& TNM3D,
                        IBuffer& IWORK,
                        bool compact = true,
                        bool herm    = true)
{
  static_assert(std::decay<TVec>::type::dimensionality == 1, " TVec::dimensionality == 1");
  static_assert((pointedType<MatA>::dimensionality == 2 or pointedType<MatA>::dimensionality == -2),
                " MatB::dimensionality == 2");
  static_assert(pointedType<MatB>::dimensionality == 2, " MatB::dimensionality == 2");
  static_assert(std::decay<MatC>::type::dimensionality == 3, " MatC::dimensionality == 3");
  static_assert(std::decay<Mat>::type::dimensionality == 3, "std::decay<Mat>::type::dimensionality == 3");

  using ma::gemmBatched;
  using ma::getrfBatched;
  using ma::getriBatched;
  using ma::H;
  using ma::T;

  int nbatch = Bi.size();
  int NMO    = (herm ? std::get<1>((*hermA[0]).sizes()) : std::get<0>((*hermA[0]).sizes()));
  int NEL    = (herm ? std::get<0>((*hermA[0]).sizes()) : std::get<1>((*hermA[0]).sizes()));

  assert(std::get<0>((*Bi[0]).sizes()) == NMO);
  assert(std::get<1>((*Bi[0]).sizes()) == NEL);
  assert(C.size() == nbatch);
  assert(std::get<2>(C.sizes()) == NMO);
  if (compact)
    assert(std::get<1>(C.sizes()) == NEL);
  else
    assert(std::get<1>(C.sizes()) == NMO);
  assert(ovlp.size() == nbatch);
  assert(std::get<1>(TNN3D.sizes()) == NEL);
  assert(std::get<2>(TNN3D.sizes()) == NEL);
  if (not compact)
  {
    assert(std::get<0>(TNM3D.sizes()) == nbatch);
    assert(std::get<1>(TNM3D.sizes()) == NEL);
    assert(std::get<2>(TNM3D.sizes()) == NMO);
  }
  assert(IWORK.num_elements() >= nbatch * (NEL + 1));
  assert(TNN3D.stride(1) == NEL); // needed by getriBatched

  using element = typename std::decay<MatC>::type::element;
  using pointer = typename std::decay<MatC>::type::element_ptr;

  int ldw = (*Bi[0]).stride(0);
  int ldN = TNN3D.stride(1);
  int ldC = C.stride(1);
  std::vector<pointer> Carray;
  std::vector<pointer> Warray;
  std::vector<pointer> NNarray;
  std::vector<decltype(&C[0])> Ci;
  //  std::vector<decltype(&TNN3D[0])> TNNi;
  Carray.reserve(nbatch);
  Warray.reserve(nbatch);
  NNarray.reserve(nbatch);
  Ci.reserve(nbatch);
  //  TNNi.reserve(nbatch);
  for (int i = 0; i < nbatch; i++)
  {
    NNarray.emplace_back(TNN3D[i].origin());
    Carray.emplace_back(C[i].origin());
    Warray.emplace_back((*Bi[i]).origin());
    Ci.emplace_back(&C[i]);
    //    TNNi.emplace_back(TNN3D[i]);
  }

  // using C for temporary storage, since getriBatched is out-of-place
  std::vector<decltype(&C[0]({0, NEL}, {0, NEL}))> Ct;
  Ct.reserve(nbatch);
  for (int i = 0; i < nbatch; i++)
    Ct.emplace_back(&C[i]({0, NEL}, {0, NEL}));

  //T(conj(A))*B
  if (herm)
    ma::BatchedProduct('N', 'N', hermA, Bi, Ct);
  else
    ma::BatchedProduct('C', 'N', hermA, Bi, Ct);

  // Invert Ct into TNN3D
  getrfBatched(NEL, Carray.data(), ldC, ma::pointer_dispatch(IWORK.origin()),
               ma::pointer_dispatch(IWORK.origin()) + nbatch * NEL, nbatch);

  using ma::strided_determinant_from_getrf;
  strided_determinant_from_getrf(NEL, ma::pointer_dispatch(Carray[0]), ldC, C.stride(0),
                                 ma::pointer_dispatch(IWORK.origin()), NEL, LogOverlapFactor, to_address(ovlp.origin()),
                                 nbatch);

  getriBatched(NEL, Carray.data(), ldC, ma::pointer_dispatch(IWORK.origin()), NNarray.data(), ldN,
               ma::pointer_dispatch(IWORK.origin()) + nbatch * NEL, nbatch);


  if (compact)
  {
    gemmBatched('T', 'T', NMO, NEL, NEL, ComplexType(1.0), Warray.data(), ldw, NNarray.data(), ldN, ComplexType(0.0),
                Carray.data(), ldC, nbatch);
  }
  else
  {
    if (herm)
    {
      int ldM = TNM3D.stride(1);
      std::vector<pointer> NMarray;
      std::vector<decltype(&TNM3D[0])> TNMi;
      NMarray.reserve(nbatch);
      TNMi.reserve(nbatch);
      for (int i = 0; i < nbatch; i++)
      {
        NMarray.emplace_back(TNM3D[i].origin());
        TNMi.emplace_back(&TNM3D[i]);
      }

      // T2 = T(T1) * T(B)
      gemmBatched('T', 'T', NMO, NEL, NEL, ComplexType(1.0), Warray.data(), ldw, NNarray.data(), ldN, ComplexType(0.0),
                  NMarray.data(), ldM, nbatch);

      // C = conj(A) * T2
      ma::BatchedProduct('T', 'N', hermA, TNMi, Ci);
    }
    else
    {
      /*
      std::vector<decltype(&TNM3D[0])> TNMi;
      TNMi.reserve(nbatch);
      for(int i=0; i<nbatch; i++) 
        TNMi.emplace_back(&TNM3D[i]);
      ma::BatchedProduct('N','C',TNNi,hermA,TNMi);
*/
      // T2 = T1 * H(A)
      for (int b = 0; b < nbatch; ++b)
        ma::product(TNN3D[b], H(*hermA[b]), TNM3D[b]);

      int ldM = TNM3D.stride(1);
      std::vector<pointer> NMarray;
      NMarray.reserve(nbatch);
      for (int i = 0; i < nbatch; i++)
        NMarray.emplace_back(TNM3D[i].origin());

      // T2 = T(T1) * T(B)
      // C = T( B * T2) = T(T2) * T(B)
      gemmBatched('T', 'T', NMO, NMO, NEL, ComplexType(1.0), Warray.data(), ldw, NMarray.data(), ldM, ComplexType(0.0),
                  Carray.data(), ldC, nbatch);
    }
  }
}

// only takes dense arrays now
template<class MatA, class MatB, class MatC, class Mat, class TVec, class IBuffer, class Tp>
void DensityMatrices(std::vector<MatA> const& Left,
                     std::vector<MatB> const& Right,
                     std::vector<MatC>& G,
                     Tp LogOverlapFactor,
                     TVec&& ovlp,
                     Mat&& TNN3D,
                     Mat&& TNM3D,
                     IBuffer& IWORK,
                     bool compact = true,
                     bool herm    = true)
{
  static_assert(std::decay<TVec>::type::dimensionality == 1, " TVec::dimensionality == 1");
  //  static_assert( (pointedType<MatA>::dimensionality == 2 or
  //                  pointedType<MatA>::dimensionality == -2), " MatA::dimensionality == 2" );
  static_assert(pointedType<MatA>::dimensionality == 2, " MatA::dimensionality == 2");
  static_assert(pointedType<MatB>::dimensionality == 2, " MatB::dimensionality == 2");
  static_assert(pointedType<MatC>::dimensionality == 2, " MatC::dimensionality == 2");
  static_assert(std::decay<Mat>::type::dimensionality == 3, "std::decay<Mat>::type::dimensionality == 3");

  using ma::batched_determinant_from_getrf;
  using ma::gemmBatched;
  using ma::getrfBatched;
  using ma::getriBatched;
  using ma::H;
  using ma::T;

  int nbatch = Right.size();
  int NMO    = (herm ? std::get<1>((*Left[0]).sizes()) : std::get<0>((*Left[0]).sizes()));
  int NEL    = (herm ? std::get<0>((*Left[0]).sizes()) : std::get<1>((*Left[0]).sizes()));

  assert(std::get<0>((*Right[0]).sizes()) == NMO);
  assert(std::get<1>((*Right[0]).sizes()) == NEL);
  assert(G.size() == nbatch);
  assert(std::get<1>((*G[0]).sizes()) == NMO);
  if (compact)
    assert((*G[0]).size() == NEL);
  else
    assert((*G[0]).size() == NMO);
  assert(ovlp.size() == nbatch);
  assert(std::get<1>(TNN3D.sizes()) == NEL);
  assert(std::get<2>(TNN3D.sizes()) == NEL);
  if (not compact)
  {
    assert(std::get<0>(TNM3D.sizes()) == nbatch);
    assert(std::get<1>(TNM3D.sizes()) == NEL);
    assert(std::get<2>(TNM3D.sizes()) == NMO);
  }
  assert(IWORK.num_elements() >= nbatch * (NEL + 1));

  using pointer = typename pointedType<MatC>::element_ptr;

  int ldR = (*Right[0]).stride(0);
  int ldL = (*Left[0]).stride(0);
  int ldN = TNN3D.stride(1);
  int ldG = (*G[0]).stride(0);
  std::vector<pointer> Garray;
  std::vector<pointer> Rarray;
  std::vector<pointer> Larray;
  std::vector<pointer> NNarray;
  Garray.reserve(nbatch);
  Rarray.reserve(nbatch);
  Larray.reserve(nbatch);
  NNarray.reserve(nbatch);
  for (int i = 0; i < nbatch; i++)
  {
    assert((*Right[i]).stride(0) == ldR);
    assert((*Left[i]).stride(0) == ldL);
    NNarray.emplace_back(TNN3D[i].origin());
    Garray.emplace_back((*G[i]).origin());
    Rarray.emplace_back((*Right[i]).origin());
    Larray.emplace_back((*Left[i]).origin());
  }

  // T(conj(A))*B
  if (herm)
    gemmBatched('N', 'N', NEL, NEL, NMO, ComplexType(1.0), Rarray.data(), ldR, Larray.data(), ldL, ComplexType(0.0),
                Garray.data(), NEL, nbatch);
  else
    gemmBatched('N', 'C', NEL, NEL, NMO, ComplexType(1.0), Rarray.data(), ldR, Larray.data(), ldL, ComplexType(0.0),
                Garray.data(), NEL, nbatch);

  // T1 = T1^(-1)
  // Invert
  getrfBatched(NEL, Garray.data(), NEL, ma::pointer_dispatch(IWORK.origin()),
               ma::pointer_dispatch(IWORK.origin()) + nbatch * NEL, nbatch);

  batched_determinant_from_getrf(NEL, Garray.data(), NEL, IWORK.origin(), NEL, LogOverlapFactor,
                                 to_address(ovlp.origin()), nbatch);

  getriBatched(NEL, Garray.data(), NEL, ma::pointer_dispatch(IWORK.origin()), ma::pointer_dispatch(NNarray.data()), ldN,
               ma::pointer_dispatch(IWORK.origin()) + nbatch * NEL, nbatch);

  if (compact)
  {
    // C = T(T1) * T(B)
    // careful with fortan ordering
    gemmBatched('T', 'T', NMO, NEL, NEL, ComplexType(1.0), Rarray.data(), ldR, NNarray.data(), ldN, ComplexType(0.0),
                Garray.data(), ldG, nbatch);
  }
  else
  {
    int ldM = TNM3D.stride(1);
    std::vector<pointer> NMarray;
    NMarray.reserve(nbatch);
    for (int i = 0; i < nbatch; i++)
      NMarray.emplace_back(TNM3D[i].origin());

    if (herm)
    {
      // T2 = T(T1) * T(B)
      gemmBatched('T', 'T', NMO, NEL, NEL, ComplexType(1.0), Rarray.data(), ldR, NNarray.data(), ldN, ComplexType(0.0),
                  NMarray.data(), ldM, nbatch);

      // C = conj(A) * T2
      gemmBatched('N', 'T', NMO, NMO, NEL, ComplexType(1.0), NMarray.data(), ldM, Larray.data(), ldL, ComplexType(0.0),
                  Garray.data(), ldG, nbatch);
    }
    else
    {
      // T2 = T1 * H(A)
      gemmBatched('C', 'N', NMO, NEL, NEL, ComplexType(1.0), Larray.data(), ldL, NNarray.data(), ldN, ComplexType(0.0),
                  NMarray.data(), ldM, nbatch);

      // T2 = T(T1) * T(B)
      // C = T( B * T2) = T(T2) * T(B)
      gemmBatched('T', 'T', NMO, NMO, NEL, ComplexType(1.0), Rarray.data(), ldR, NMarray.data(), ldM, ComplexType(0.0),
                  Garray.data(), ldG, nbatch);
    }
  }
}

template<class MatA, class MatB, class Mat, class TVec, class IBuffer, class Tp>
void Overlap(std::vector<MatA>& hermA,
             std::vector<MatB>& Bi,
             Tp LogOverlapFactor,
             TVec&& ovlp,
             Mat&& TNN3D,
             IBuffer& IWORK,
             bool herm = true)
{
  static_assert((pointedType<MatA>::dimensionality == 2 or pointedType<MatA>::dimensionality == -2),
                " MatA::dimensionality == 2");
  static_assert(pointedType<MatB>::dimensionality == 2, " MatB::dimensionality == 2");
  static_assert(std::decay<TVec>::type::dimensionality == 1, " TVec::dimensionality == 1");
  static_assert(std::decay<Mat>::type::dimensionality == 3, "std::decay<Mat>::type::dimensionality == 3");

  using ma::gemmBatched;
  using ma::getrfBatched;
  using ma::H;
  using ma::T;

  int nbatch = Bi.size();
  assert(hermA.size() >= nbatch);
  int NMO = (herm ? std::get<1>((*hermA[0]).sizes()) : std::get<0>((*hermA[0]).sizes()));
  int NEL = (herm ? std::get<0>((*hermA[0]).sizes()) : std::get<1>((*hermA[0]).sizes()));

  assert(std::get<0>((*Bi[0]).sizes()) == NMO);
  assert(std::get<1>((*Bi[0]).sizes()) == NEL);
  assert(ovlp.size() == nbatch);
  assert(std::get<1>(TNN3D.sizes()) == NEL);
  assert(std::get<2>(TNN3D.sizes()) == NEL);
  assert(IWORK.num_elements() >= nbatch * (NEL + 1));

  using pointer = typename std::decay<Mat>::type::element_ptr;

  int ldw = (*Bi[0]).stride(0);
  int ldN = TNN3D.stride(1);
  std::vector<pointer> Warray;
  std::vector<pointer> NNarray;
  std::vector<decltype(&TNN3D[0])> Ci;
  Warray.reserve(nbatch);
  NNarray.reserve(nbatch);
  Ci.reserve(nbatch);
  for (int i = 0; i < nbatch; i++)
  {
    NNarray.emplace_back(TNN3D[i].origin());
    Warray.emplace_back((*Bi[i]).origin());
    Ci.emplace_back(&TNN3D[i]);
  }

  // T(conj(A))*B
  if (herm)
  {
    ma::BatchedProduct('N', 'N', hermA, Bi, Ci);
  }
  else
  {
    ma::BatchedProduct('C', 'N', hermA, Bi, Ci);
  }

  // T1 = T1^(-1)
  // Invert
  getrfBatched(NEL, NNarray.data(), ldN, IWORK.origin(), IWORK.origin() + nbatch * NEL, nbatch);

  using ma::strided_determinant_from_getrf;
  strided_determinant_from_getrf(NEL, NNarray[0], ldN, TNN3D.stride(0), IWORK.origin(), NEL, LogOverlapFactor,
                                 to_address(ovlp.origin()), nbatch);
}


} // namespace batched

} // namespace SlaterDeterminantOperations

} // namespace afqmc

} // namespace qmcplusplus

#endif
