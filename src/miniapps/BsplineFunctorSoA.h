//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp. 
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_BSPLINE2_FUNCTOR_H
#define QMCPLUSPLUS_BSPLINE2_FUNCTOR_H
#include "Numerics/OptimizableFunctorBase.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "Numerics/LinearFit.h"
#include "simd/allocator.hpp"
#include <cstdio>

namespace qmcplusplus {

template<class T>
struct BsplineFunctorSoA
#if QMC_BUILD_LEVEL < 5
: public OptimizableFunctorBase
#endif
{
  //real_type is not useful but keep it here
#if QMC_BUILD_LEVEL==5
  using real_type=T;
#endif
  const TinyVector<T,16> A, dA, d2A, d3A;
  aligned_vector<T> SplineCoefs;
  int NumParams;
  int sizeSplineCoefs;
  real_type cutoff_radius;
  real_type DeltaR, DeltaRInv;
  real_type CuspValue;
  real_type Y, dY, d2Y;
   
  ///Stores the derivatives w.r.t. SplineCoefs of the u, du/dr, and d2u/dr2
  std::vector<TinyVector<real_type,3> > SplineDerivs;
  std::vector<real_type> Parameters;
  std::vector<std::string> ParameterNames;
  std::string elementType, pairType;
  std::string fileName;

  int ResetCount;
  int ReportLevel;
  bool notOpt;
  bool periodic;

  ///constructor
  BsplineFunctorSoA(real_type cusp=real_type(0));

  /*@{ basic evaluation functions used by Jastrow functions before SoA distance tables */
  real_type evaluate(real_type r) const;
  real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2) const;
  real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2, real_type &d3udr3) const;
  bool evaluateDerivatives(real_type r, std::vector<TinyVector<real_type,3> >& derivs);
  bool evaluateDerivatives(real_type r, std::vector<real_type>& derivs) const;
  /*@}*/

  /** compute value, gradient and laplacian for [iStart, iEnd) pairs
   * @param iStart starting particle index
   * @param iEnd ending particle index
   * @param _distArray distance arrUay
   * @param _valArray  u(r_j) for j=[iStart,iEnd)
   * @param _gradArray  du(r_j)/dr /r_j for j=[iStart,iEnd)
   * @param _lapArray  d2u(r_j)/dr2 for j=[iStart,iEnd)
   * @param distArrayCompressed temp storage to filter r_j < cutoff_radius
   * @param distIndices temp storage for the compressed index
   */
  void evaluateVGL(const int iStart, const int iEnd, 
      const real_type* _distArray,  
      real_type* restrict _valArray,
      real_type* restrict _gradArray, 
      real_type* restrict _laplArray, 
      real_type* restrict distArrayCompressed, int* restrict distIndices ) const;

  /** evaluate sum of the pair potentials for [iStart,iEnd)
   * @param iStart starting particle index
   * @param iEnd ending particle index
   * @param _distArray distance arrUay
   * @param distArrayCompressed temp storage to filter r_j < cutoff_radius
   * @return \f$\sum u(r_j)\f$ for r_j < cutoff_radius
   */
  T evaluateU(const int iStart, const int iEnd, 
      const real_type* restrict _distArray, 
      real_type* restrict distArrayCompressed) const;

  void resize(int n);
  void reset();
  void initialize(int numPoints, std::vector<real_type>& x, std::vector<real_type>& y
                  , real_type cusp, real_type rcut, std::string& id, std::string& optimize );

#if QMC_BUILD_LEVEL<5
  bool put(xmlNodePtr cur);
  inline real_type evaluate(real_type r, real_type rinv) 
  { return Y=evaluate(r,dY,d2Y); }

  inline void evaluateAll(real_type r, real_type rinv)
  { Y=evaluate(r,dY,d2Y); }

  OptimizableFunctorBase* makeClone() const
  {
    return new BsplineFunctorSoA(*this);
  }

  void setReportLevel(int i, const std::string& fname);
  void reportStatus(ostream& os);
  void print();
  void print(std::ostream& os);
#endif
};

template<typename T>
  BsplineFunctorSoA<T>::BsplineFunctorSoA(real_type cusp): NumParams(0),
  A(-1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
      3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
      -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
      1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0),
  dA(0.0, -0.5,  1.0, -0.5,
      0.0,  1.5, -2.0,  0.0,
      0.0, -1.5,  1.0,  0.5,
      0.0,  0.5,  0.0,  0.0),
  d2A(0.0, 0.0, -1.0,  1.0,
      0.0, 0.0,  3.0, -2.0,
      0.0, 0.0, -3.0,  1.0,
      0.0, 0.0,  1.0,  0.0),
  d3A(0.0, 0.0,  0.0, -1.0,
      0.0, 0.0,  0.0,  3.0,
      0.0, 0.0,  0.0, -3.0,
      0.0, 0.0,  0.0,  1.0),
  CuspValue(cusp), ResetCount(0), ReportLevel(0), notOpt(false), periodic(true)
{
  cutoff_radius = real_type(0);
}


template<typename T> inline void BsplineFunctorSoA<T>::resize(int n)
{
  NumParams = n;
  int numCoefs = NumParams + 4;
  int numKnots = numCoefs - 2;
  DeltaR = cutoff_radius / (real_type)(numKnots - 1);
  DeltaRInv = 1.0/DeltaR;
  Parameters.resize(n);

  int nc=getAlignedSize<real_type>(numCoefs);
  SplineCoefs.resize(nc);
  SplineDerivs.resize(nc);

  std::fill(SplineCoefs.begin(),SplineCoefs.end(),real_type(0));
  sizeSplineCoefs = numCoefs;
}

template<typename T>
inline void BsplineFunctorSoA<T>::reset()
{
  int numCoefs = NumParams + 4;
  int numKnots = numCoefs - 2;
  DeltaR = cutoff_radius / (real_type)(numKnots - 1);
  DeltaRInv = 1.0/DeltaR;
  std::fill(SplineCoefs.begin(),SplineCoefs.end(),real_type(0));
  // Ensure that cusp conditions is satsified at the origin
  SplineCoefs[1] = Parameters[0];
  SplineCoefs[2] = Parameters[1];
  SplineCoefs[0] = Parameters[1] - 2.0*DeltaR * CuspValue;
  for (int i=2; i<Parameters.size(); i++)
    SplineCoefs[i+1] = Parameters[i];
}

template<typename T>
inline T
BsplineFunctorSoA<T>::evaluate(real_type r) const
{
  if (r >= cutoff_radius)
    return 0.0;

  r *= DeltaRInv;
  real_type ipart, t;
  t = std::modf(r, &ipart);
  int i = (int) ipart;
  real_type tp[4];
  tp[0] = t*t*t;
  tp[1] = t*t;
  tp[2] = t;
  tp[3] = 1.0;
  return
    (SplineCoefs[i+0]*(A[ 0]*tp[0] + A[ 1]*tp[1] + A[ 2]*tp[2] + A[ 3]*tp[3])+
     SplineCoefs[i+1]*(A[ 4]*tp[0] + A[ 5]*tp[1] + A[ 6]*tp[2] + A[ 7]*tp[3])+
     SplineCoefs[i+2]*(A[ 8]*tp[0] + A[ 9]*tp[1] + A[10]*tp[2] + A[11]*tp[3])+
     SplineCoefs[i+3]*(A[12]*tp[0] + A[13]*tp[1] + A[14]*tp[2] + A[15]*tp[3]));
}


template<typename T> inline T 
BsplineFunctorSoA<T>::evaluate(real_type r, real_type& dudr, real_type& d2udr2) const
{
  constexpr T czero(0);
  if (r >= cutoff_radius)
  {
    dudr = d2udr2 = czero;
    return czero;
  }
  r *= DeltaRInv;
  real_type ipart, t;
  t = std::modf(r, &ipart);
  int i = (int) ipart;
  real_type tp[4];
  tp[0] = t*t*t;
  tp[1] = t*t;
  tp[2] = t;
  tp[3] = 1.0;
  d2udr2 = DeltaRInv * DeltaRInv *
    (SplineCoefs[i+0]*(d2A[ 0]*tp[0] + d2A[ 1]*tp[1] + d2A[ 2]*tp[2] + d2A[ 3]*tp[3])+
     SplineCoefs[i+1]*(d2A[ 4]*tp[0] + d2A[ 5]*tp[1] + d2A[ 6]*tp[2] + d2A[ 7]*tp[3])+
     SplineCoefs[i+2]*(d2A[ 8]*tp[0] + d2A[ 9]*tp[1] + d2A[10]*tp[2] + d2A[11]*tp[3])+
     SplineCoefs[i+3]*(d2A[12]*tp[0] + d2A[13]*tp[1] + d2A[14]*tp[2] + d2A[15]*tp[3]));
  dudr = DeltaRInv *
    (SplineCoefs[i+0]*(dA[ 0]*tp[0] + dA[ 1]*tp[1] + dA[ 2]*tp[2] + dA[ 3]*tp[3])+
     SplineCoefs[i+1]*(dA[ 4]*tp[0] + dA[ 5]*tp[1] + dA[ 6]*tp[2] + dA[ 7]*tp[3])+
     SplineCoefs[i+2]*(dA[ 8]*tp[0] + dA[ 9]*tp[1] + dA[10]*tp[2] + dA[11]*tp[3])+
     SplineCoefs[i+3]*(dA[12]*tp[0] + dA[13]*tp[1] + dA[14]*tp[2] + dA[15]*tp[3]));
  return
    (SplineCoefs[i+0]*(A[ 0]*tp[0] + A[ 1]*tp[1] + A[ 2]*tp[2] + A[ 3]*tp[3])+
     SplineCoefs[i+1]*(A[ 4]*tp[0] + A[ 5]*tp[1] + A[ 6]*tp[2] + A[ 7]*tp[3])+
     SplineCoefs[i+2]*(A[ 8]*tp[0] + A[ 9]*tp[1] + A[10]*tp[2] + A[11]*tp[3])+
     SplineCoefs[i+3]*(A[12]*tp[0] + A[13]*tp[1] + A[14]*tp[2] + A[15]*tp[3]));
}

template<typename T> inline T 
BsplineFunctorSoA<T>::evaluate(real_type r, real_type& dudr, real_type& d2udr2, real_type &d3udr3) const
{
  if (r >= cutoff_radius)
  {
    dudr = d2udr2 = d3udr3 = 0.0;
    return 0.0;
  }
  r *= DeltaRInv;
  real_type ipart, t;
  t = std::modf(r, &ipart);
  int i = (int) ipart;
  real_type tp[4];
  tp[0] = t*t*t;
  tp[1] = t*t;
  tp[2] = t;
  tp[3] = 1.0;
  d3udr3 = DeltaRInv * DeltaRInv * DeltaRInv *
    (SplineCoefs[i+0]*(d3A[ 0]*tp[0] + d3A[ 1]*tp[1] + d3A[ 2]*tp[2] + d3A[ 3]*tp[3])+
     SplineCoefs[i+1]*(d3A[ 4]*tp[0] + d3A[ 5]*tp[1] + d3A[ 6]*tp[2] + d3A[ 7]*tp[3])+
     SplineCoefs[i+2]*(d3A[ 8]*tp[0] + d3A[ 9]*tp[1] + d3A[10]*tp[2] + d3A[11]*tp[3])+
     SplineCoefs[i+3]*(d3A[12]*tp[0] + d3A[13]*tp[1] + d3A[14]*tp[2] + d3A[15]*tp[3]));
  d2udr2 = DeltaRInv * DeltaRInv *
    (SplineCoefs[i+0]*(d2A[ 0]*tp[0] + d2A[ 1]*tp[1] + d2A[ 2]*tp[2] + d2A[ 3]*tp[3])+
     SplineCoefs[i+1]*(d2A[ 4]*tp[0] + d2A[ 5]*tp[1] + d2A[ 6]*tp[2] + d2A[ 7]*tp[3])+
     SplineCoefs[i+2]*(d2A[ 8]*tp[0] + d2A[ 9]*tp[1] + d2A[10]*tp[2] + d2A[11]*tp[3])+
     SplineCoefs[i+3]*(d2A[12]*tp[0] + d2A[13]*tp[1] + d2A[14]*tp[2] + d2A[15]*tp[3]));
  dudr = DeltaRInv *
    (SplineCoefs[i+0]*(dA[ 0]*tp[0] + dA[ 1]*tp[1] + dA[ 2]*tp[2] + dA[ 3]*tp[3])+
     SplineCoefs[i+1]*(dA[ 4]*tp[0] + dA[ 5]*tp[1] + dA[ 6]*tp[2] + dA[ 7]*tp[3])+
     SplineCoefs[i+2]*(dA[ 8]*tp[0] + dA[ 9]*tp[1] + dA[10]*tp[2] + dA[11]*tp[3])+
     SplineCoefs[i+3]*(dA[12]*tp[0] + dA[13]*tp[1] + dA[14]*tp[2] + dA[15]*tp[3]));
  return
    (SplineCoefs[i+0]*(A[ 0]*tp[0] + A[ 1]*tp[1] + A[ 2]*tp[2] + A[ 3]*tp[3])+
     SplineCoefs[i+1]*(A[ 4]*tp[0] + A[ 5]*tp[1] + A[ 6]*tp[2] + A[ 7]*tp[3])+
     SplineCoefs[i+2]*(A[ 8]*tp[0] + A[ 9]*tp[1] + A[10]*tp[2] + A[11]*tp[3])+
     SplineCoefs[i+3]*(A[12]*tp[0] + A[13]*tp[1] + A[14]*tp[2] + A[15]*tp[3]));
}

template<typename T>
inline bool
BsplineFunctorSoA<T>::evaluateDerivatives(real_type r, std::vector<TinyVector<real_type,3> >& derivs) 
{
  if (r >= cutoff_radius) return false;
  r *= DeltaRInv;
  real_type ipart, t;
  t = std::modf(r, &ipart);
  int i = (int) ipart;
  real_type tp[4];
  tp[0] = t*t*t;
  tp[1] = t*t;
  tp[2] = t;
  tp[3] = 1.0;

  SplineDerivs[0] = TinyVector<real_type,3>(0.0);
  // d/dp_i u(r)
  SplineDerivs[i+0][0] = A[ 0]*tp[0] + A[ 1]*tp[1] + A[ 2]*tp[2] + A[ 3]*tp[3];
  SplineDerivs[i+1][0] = A[ 4]*tp[0] + A[ 5]*tp[1] + A[ 6]*tp[2] + A[ 7]*tp[3];
  SplineDerivs[i+2][0] = A[ 8]*tp[0] + A[ 9]*tp[1] + A[10]*tp[2] + A[11]*tp[3];
  SplineDerivs[i+3][0] = A[12]*tp[0] + A[13]*tp[1] + A[14]*tp[2] + A[15]*tp[3];
  // d/dp_i du/dr
  SplineDerivs[i+0][1] = DeltaRInv * (dA[ 1]*tp[1] + dA[ 2]*tp[2] + dA[ 3]*tp[3]);
  SplineDerivs[i+1][1] = DeltaRInv * (dA[ 5]*tp[1] + dA[ 6]*tp[2] + dA[ 7]*tp[3]);
  SplineDerivs[i+2][1] = DeltaRInv * (dA[ 9]*tp[1] + dA[10]*tp[2] + dA[11]*tp[3]);
  SplineDerivs[i+3][1] = DeltaRInv * (dA[13]*tp[1] + dA[14]*tp[2] + dA[15]*tp[3]);
  // d/dp_i d2u/dr2
  SplineDerivs[i+0][2] = DeltaRInv * DeltaRInv * (d2A[ 2]*tp[2] + d2A[ 3]*tp[3]);
  SplineDerivs[i+1][2] = DeltaRInv * DeltaRInv * (d2A[ 6]*tp[2] + d2A[ 7]*tp[3]);
  SplineDerivs[i+2][2] = DeltaRInv * DeltaRInv * (d2A[10]*tp[2] + d2A[11]*tp[3]);
  SplineDerivs[i+3][2] = DeltaRInv * DeltaRInv * (d2A[14]*tp[2] + d2A[15]*tp[3]);

  int imin=std::max(i,1);
  int imax=std::min(i+4,NumParams+1);
  for (int n=imin; n<imax; ++n)
    derivs[n-1] = SplineDerivs[n];
  derivs[1]+=SplineDerivs[0];

  return true;
}

template<typename T>
  inline bool 
BsplineFunctorSoA<T>::evaluateDerivatives(real_type r, std::vector<real_type>& derivs) const
{
  if (r >= cutoff_radius) return false;
  real_type tp[4],v[4],ipart,t;
  t = std::modf(r*DeltaRInv, &ipart);
  tp[0] = t*t*t;
  tp[1] = t*t;
  tp[2] = t;
  tp[3] = 1.0;
  v[0] = A[ 0]*tp[0] + A[ 1]*tp[1] + A[ 2]*tp[2] + A[ 3]*tp[3];
  v[1] = A[ 4]*tp[0] + A[ 5]*tp[1] + A[ 6]*tp[2] + A[ 7]*tp[3];
  v[2] = A[ 8]*tp[0] + A[ 9]*tp[1] + A[10]*tp[2] + A[11]*tp[3];
  v[3] = A[12]*tp[0] + A[13]*tp[1] + A[14]*tp[2] + A[15]*tp[3];
  int i = (int) ipart;
  int imin=std::max(i,1);
  int imax=std::min(i+4,NumParams+1)-1;
  int n=imin-1, j=imin-i;
  while(n<imax && j<4)
  {
    derivs[n] = v[j];
    n++; j++;
  }
  if(i==0) derivs[1]+= v[0];
  return true;
}

template<typename T>
inline T 
BsplineFunctorSoA<T>::evaluateU(const int iStart, const int iEnd,
    const real_type* restrict _distArray, real_type* restrict distArrayCompressed ) const
{
  const real_type* restrict distArray = _distArray + iStart;

  ASSUME_ALIGNED(distArrayCompressed);
  int iCount = 0;
  const int iLimit = iEnd-iStart;

#pragma vector always 
  for ( int jat = 0; jat < iLimit; jat++ ) {
    real_type r = distArray[jat];
    if ( r < cutoff_radius )
      distArrayCompressed[iCount++] = distArray[jat];
  }

  real_type d = 0.0;
#pragma simd reduction (+:d )
  for ( int jat = 0; jat < iCount; jat++ ) {
    real_type r = distArrayCompressed[jat];
    r *= DeltaRInv;
    int i = (int)r;
    real_type t = r - real_type(i);
    real_type tp0 = t*t*t;
    real_type tp1 = t*t;
    real_type tp2 = t;

    real_type d1 = SplineCoefs[i+0]*(A[ 0]*tp0 + A[ 1]*tp1 + A[ 2]*tp2 + A[ 3]);
    real_type d2 = SplineCoefs[i+1]*(A[ 4]*tp0 + A[ 5]*tp1 + A[ 6]*tp2 + A[ 7]);
    real_type d3 = SplineCoefs[i+2]*(A[ 8]*tp0 + A[ 9]*tp1 + A[10]*tp2 + A[11]);
    real_type d4 = SplineCoefs[i+3]*(A[12]*tp0 + A[13]*tp1 + A[14]*tp2 + A[15]);
    d += ( d1 + d2 + d3 + d4 );
  }
  return d;
}

template<typename T> 
inline void BsplineFunctorSoA<T>::evaluateVGL(const int iStart, const int iEnd, 
    const real_type* _distArray,  real_type* restrict _valArray, 
    real_type* restrict _gradArray, real_type* restrict _laplArray, 
    real_type* restrict distArrayCompressed, int* restrict distIndices ) const
{

  real_type dSquareDeltaRinv = DeltaRInv * DeltaRInv;
  constexpr real_type cZero(0); 
  constexpr real_type cOne(1);
  constexpr real_type cMOne(-1); 

  //    START_MARK_FIRST();

  ASSUME_ALIGNED(distIndices);
  ASSUME_ALIGNED(distArrayCompressed);
  int iCount = 0;
  int iLimit = iEnd-iStart;
  const real_type* distArray = _distArray + iStart;
  real_type* valArray = _valArray + iStart;
  real_type* gradArray = _gradArray + iStart;
  real_type* laplArray = _laplArray + iStart;

#pragma vector always
  for ( int jat = 0; jat < iLimit; jat++ ) {
    real_type r = distArray[jat];
    if ( r < cutoff_radius ) {
      distIndices[iCount] = jat;
      distArrayCompressed[iCount] = r;
      iCount++;
    }
  }

#pragma omp simd 
  for ( int j = 0; j < iCount; j++ ) {

    real_type r = distArrayCompressed[j];
    int iScatter = distIndices[j]; 
    real_type rinv = cOne/r; 
    r *= DeltaRInv;
    int iGather = (int)r;
    real_type t = r - real_type(iGather);
    real_type tp0 = t*t*t;
    real_type tp1 = t*t;
    real_type tp2 = t;

    real_type sCoef0 = SplineCoefs[iGather+0];
    real_type sCoef1 = SplineCoefs[iGather+1];
    real_type sCoef2 = SplineCoefs[iGather+2];
    real_type sCoef3 = SplineCoefs[iGather+3];

    laplArray[iScatter] = dSquareDeltaRinv *
      (sCoef0*( d2A[ 2]*tp2 + d2A[ 3])+
       sCoef1*( d2A[ 6]*tp2 + d2A[ 7])+
       sCoef2*( d2A[10]*tp2 + d2A[11])+
       sCoef3*( d2A[14]*tp2 + d2A[15]));

    gradArray[iScatter] = DeltaRInv * rinv *
      (sCoef0*( dA[ 1]*tp1 + dA[ 2]*tp2 + dA[ 3])+
       sCoef1*( dA[ 5]*tp1 + dA[ 6]*tp2 + dA[ 7])+
       sCoef2*( dA[ 9]*tp1 + dA[10]*tp2 + dA[11])+
       sCoef3*( dA[13]*tp1 + dA[14]*tp2 + dA[15]));

    valArray[iScatter] = (sCoef0*(A[ 0]*tp0 + A[ 1]*tp1 + A[ 2]*tp2 + A[ 3])+
        sCoef1*(A[ 4]*tp0 + A[ 5]*tp1 + A[ 6]*tp2 + A[ 7])+
        sCoef2*(A[ 8]*tp0 + A[ 9]*tp1 + A[10]*tp2 + A[11])+
        sCoef3*(A[12]*tp0 + A[13]*tp1 + A[14]*tp2 + A[15]));
  }

#pragma simd 
  for ( int j = 0; j < iLimit; j++ ) {
    if ( distArray[j] > cutoff_radius ) {
      valArray[j] = cZero;
      gradArray[j] = cZero;
      laplArray[j] = cZero;
    }
  }
}


template<typename T> void 
  BsplineFunctorSoA<T>::initialize(int numPoints, std::vector<real_type>& x, std::vector<real_type>& y
                  , real_type cusp, real_type rcut, std::string& id, std::string& optimize )
  {
    ReportEngine PRE("BsplineFunctorSoA","initialize");
    NumParams = numPoints;
    cutoff_radius = rcut;
    CuspValue = cusp;
    if (NumParams == 0)
    {
      PRE.error("You must specify a positive number of parameters for the Bspline jastrow function.",true);
    }
    app_log() << "Initializing BsplineFunctorSoA from array. \n";
    app_log() << " size = " << NumParams << " parameters " << endl;
    app_log() << " cusp = " << CuspValue << endl;
    app_log() << " rcut = " << cutoff_radius << endl;
    resize(NumParams);
    int npts = x.size();
    Matrix<real_type> basis(npts,NumParams);
    std::vector<TinyVector<real_type,3> > derivs(NumParams);
    for (int i=0; i<npts; i++)
    {
      real_type r = x[i];
      if (r > cutoff_radius)
      {
        PRE.error("Error in BsplineFunctorSoA::initialize: r > cutoff_radius.",true);
      }
      evaluateDerivatives(r, derivs);
      for (int j=0; j<NumParams; j++)
        basis(i,j) = derivs[j][0];
    }
    resize(NumParams);
    LinearFit(y, basis, Parameters);
    app_log() << "New parameters are:\n";
    for (int i=0; i < Parameters.size(); i++)
      app_log() << "   " << Parameters[i] << endl;
#if QMC_BUILD_LEVEL < 5
    if(optimize == "yes")
    {
      // Setup parameter names
      for (int i=0; i< NumParams; i++)
      {
        std::stringstream sstr;
        sstr << id << "_" << i;
        myVars.insert(sstr.str(),Parameters[i],true,optimize::LOGLINEAR_P);
      }
      app_log() << "Parameter     Name      Value\n";
      myVars.print(app_log());
    }
    else
#endif
    {
      notOpt=true;
      app_log() << "Parameters of BsplineFunctorSoA id:"
                <<id <<" are not being optimized.\n";
    }
    reset();
  }
}
#endif
