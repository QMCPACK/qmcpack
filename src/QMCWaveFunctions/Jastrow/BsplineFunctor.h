//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: John R. Gergely,  University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_BSPLINE_FUNCTOR_H
#define QMCPLUSPLUS_BSPLINE_FUNCTOR_H
#include "Numerics/OptimizableFunctorBase.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "Numerics/LinearFit.h"
#include "simd/allocator.hpp"
#include <cstdio>

namespace qmcplusplus
{

template<class T>
struct BsplineFunctor: public OptimizableFunctorBase
{

  typedef real_type value_type;
  int NumParams;
  int Dummy;
  const real_type A[16], dA[16], d2A[16], d3A[16];
  aligned_vector<real_type> SplineCoefs;

  //static const real_type A[16], dA[16], d2A[16];
  real_type DeltaR, DeltaRInv;
  real_type CuspValue;
  real_type Y, dY, d2Y;
  // Stores the derivatives w.r.t. SplineCoefs
  // of the u, du/dr, and d2u/dr2
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
  BsplineFunctor(real_type cusp=0.0) :
    NumParams(0),
    A{-1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
       3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
      -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
       1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0},
    dA{0.0, -0.5,  1.0, -0.5,
       0.0,  1.5, -2.0,  0.0,
       0.0, -1.5,  1.0,  0.5,
       0.0,  0.5,  0.0,  0.0},
    d2A{0.0, 0.0, -1.0,  1.0,
        0.0, 0.0,  3.0, -2.0,
        0.0, 0.0, -3.0,  1.0,
        0.0, 0.0,  1.0,  0.0},
    d3A{0.0, 0.0,  0.0, -1.0,
        0.0, 0.0,  0.0,  3.0,
        0.0, 0.0,  0.0, -3.0,
        0.0, 0.0,  0.0,  1.0},
    CuspValue(cusp), ResetCount(0), ReportLevel(0), notOpt(false), periodic(true)
  {
    cutoff_radius = 0.0;
  }

  OptimizableFunctorBase* makeClone() const
  {
    return new BsplineFunctor(*this);
  }

  inline void setReportLevel(int i, const std::string& fname)
  {
    ReportLevel=i;
    fileName=(ReportLevel)? fname:"0";
  }

  void resize(int n)
  {
    NumParams = n;
    int numCoefs = NumParams + 4;
    int numKnots = numCoefs - 2;
    DeltaR = cutoff_radius / (real_type)(numKnots - 1);
    DeltaRInv = 1.0/DeltaR;
    Parameters.resize(n);
    SplineCoefs.resize(numCoefs);
    SplineDerivs.resize(numCoefs);
  }

  void reset()
  {
    int numCoefs = NumParams + 4;
    int numKnots = numCoefs - 2;
    DeltaR = cutoff_radius / (real_type)(numKnots - 1);
    DeltaRInv = 1.0/DeltaR;
    for (int i=0; i<SplineCoefs.size(); i++)
      SplineCoefs[i] = 0.0;
    // Ensure that cusp conditions is satsified at the origin
    SplineCoefs[1] = Parameters[0];
    SplineCoefs[2] = Parameters[1];
    SplineCoefs[0] = Parameters[1] - 2.0*DeltaR * CuspValue;
    for (int i=2; i<Parameters.size(); i++)
      SplineCoefs[i+1] = Parameters[i];
  }

  /** compute value, gradient and laplacian for [iStart, iEnd) pairs
   * @param iat dummy
   * @param iStart starting particle index
   * @param iEnd ending particle index
   * @param _distArray distance arrUay
   * @param _valArray  u(r_j) for j=[iStart,iEnd)
   * @param _gradArray  du(r_j)/dr /r_j for j=[iStart,iEnd)
   * @param _lapArray  d2u(r_j)/dr2 for j=[iStart,iEnd)
   * @param distArrayCompressed temp storage to filter r_j < cutoff_radius
   * @param distIndices temp storage for the compressed index
   */
  void evaluateVGL(const int iat, const int iStart, const int iEnd,
      const T* _distArray,  
      T* restrict _valArray,
      T* restrict _gradArray, 
      T* restrict _laplArray, 
      T* restrict distArrayCompressed, int* restrict distIndices ) const;

  /** evaluate sum of the pair potentials for [iStart,iEnd)
   * @param iat dummy
   * @param iStart starting particle index
   * @param iEnd ending particle index
   * @param _distArray distance arrUay
   * @param distArrayCompressed temp storage to filter r_j < cutoff_radius
   * @return \f$\sum u(r_j)\f$ for r_j < cutoff_radius
   */
  T evaluateV(const int iat, const int iStart, const int iEnd,
      const T* restrict _distArray, 
      T* restrict distArrayCompressed) const;

  inline real_type evaluate(real_type r)
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
  inline real_type evaluate(real_type r, real_type rinv)
  {
    return Y=evaluate(r,dY,d2Y);
  }

  inline void evaluateAll(real_type r, real_type rinv)
  {
    Y=evaluate(r,dY,d2Y);
  }

  inline real_type
  evaluate(real_type r, real_type& dudr, real_type& d2udr2)
  {
    if (r >= cutoff_radius)
    {
      dudr = d2udr2 = 0.0;
      return 0.0;
    }
//       real_type eps = 1.0e-5;
//       real_type dudr_FD = (evaluate(r+eps)-evaluate(r-eps))/(2.0*eps);
//       real_type d2udr2_FD = (evaluate(r+eps)+evaluate(r-eps)-2.0*evaluate(r))/(eps*eps);
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
//       if (std::abs(dudr_FD-dudr) > 1.0e-8)
//  std::cerr << "Error in BsplineFunction:  dudr = " << dudr
//       << "  dudr_FD = " << dudr_FD << std::endl;
//       if (std::abs(d2udr2_FD-d2udr2) > 1.0e-4)
//  std::cerr << "Error in BsplineFunction:  r = " << r << "  d2udr2 = " << dudr
//       << "  d2udr2_FD = " << d2udr2_FD << "  rcut = " << cutoff_radius << std::endl;
    return
      (SplineCoefs[i+0]*(A[ 0]*tp[0] + A[ 1]*tp[1] + A[ 2]*tp[2] + A[ 3]*tp[3])+
       SplineCoefs[i+1]*(A[ 4]*tp[0] + A[ 5]*tp[1] + A[ 6]*tp[2] + A[ 7]*tp[3])+
       SplineCoefs[i+2]*(A[ 8]*tp[0] + A[ 9]*tp[1] + A[10]*tp[2] + A[11]*tp[3])+
       SplineCoefs[i+3]*(A[12]*tp[0] + A[13]*tp[1] + A[14]*tp[2] + A[15]*tp[3]));
  }


  inline real_type
  evaluate(real_type r, real_type& dudr, real_type& d2udr2, real_type &d3udr3)
  {
    if (r >= cutoff_radius)
    {
      dudr = d2udr2 = d3udr3 = 0.0;
      return 0.0;
    }
    // real_type eps = 1.0e-5;
//       real_type dudr_FD = (evaluate(r+eps)-evaluate(r-eps))/(2.0*eps);
//       real_type d2udr2_FD = (evaluate(r+eps)+evaluate(r-eps)-2.0*evaluate(r))/(eps*eps);
    // real_type d3udr3_FD = (-1.0*evaluate(r+1.0*eps)
    //         +2.0*evaluate(r+0.5*eps)
    //         -2.0*evaluate(r-0.5*eps)
    //         +1.0*evaluate(r-1.0*eps))/(eps*eps*eps);
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
//       if (std::abs(dudr_FD-dudr) > 1.0e-8)
//  std::cerr << "Error in BsplineFunction:  dudr = " << dudr
//       << "  dudr_FD = " << dudr_FD << std::endl;
//       if (std::abs(d2udr2_FD-d2udr2) > 1.0e-4)
//  std::cerr << "Error in BsplineFunction:  r = " << r << "  d2udr2 = " << dudr
//       << "  d2udr2_FD = " << d2udr2_FD << "  rcut = " << cutoff_radius << std::endl;
    // if (std::abs(d3udr3_FD-d3udr3) > 1.0e-4)
    //  std::cerr << "Error in BsplineFunction:  r = " << r << "  d3udr3 = " << dudr
    //       << "  d3udr3_FD = " << d3udr3_FD << "  rcut = " << cutoff_radius << std::endl;
    return
      (SplineCoefs[i+0]*(A[ 0]*tp[0] + A[ 1]*tp[1] + A[ 2]*tp[2] + A[ 3]*tp[3])+
       SplineCoefs[i+1]*(A[ 4]*tp[0] + A[ 5]*tp[1] + A[ 6]*tp[2] + A[ 7]*tp[3])+
       SplineCoefs[i+2]*(A[ 8]*tp[0] + A[ 9]*tp[1] + A[10]*tp[2] + A[11]*tp[3])+
       SplineCoefs[i+3]*(A[12]*tp[0] + A[13]*tp[1] + A[14]*tp[2] + A[15]*tp[3]));
  }



  inline bool
  evaluateDerivatives(real_type r, std::vector<TinyVector<real_type,3> >& derivs)
  {
    if (r >= cutoff_radius)
      return false;
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

    //real_type v[4],dv[4],d2v[4];
    //v[0] = A[ 0]*tp[0] + A[ 1]*tp[1] + A[ 2]*tp[2] + A[ 3]*tp[3];
    //v[1] = A[ 4]*tp[0] + A[ 5]*tp[1] + A[ 6]*tp[2] + A[ 7]*tp[3];
    //v[2] = A[ 8]*tp[0] + A[ 9]*tp[1] + A[10]*tp[2] + A[11]*tp[3];
    //v[3] = A[12]*tp[0] + A[13]*tp[1] + A[14]*tp[2] + A[15]*tp[3];
    //// d/dp_i du/dr
    //dv[0] = DeltaRInv * (dA[ 1]*tp[1] + dA[ 2]*tp[2] + dA[ 3]*tp[3]);
    //dv[1] = DeltaRInv * (dA[ 5]*tp[1] + dA[ 6]*tp[2] + dA[ 7]*tp[3]);
    //dv[2] = DeltaRInv * (dA[ 9]*tp[1] + dA[10]*tp[2] + dA[11]*tp[3]);
    //dv[3] = DeltaRInv * (dA[13]*tp[1] + dA[14]*tp[2] + dA[15]*tp[3]);
    //// d/dp_i d2u/dr2
    //d2v[0] = DeltaRInv * DeltaRInv * (d2A[ 2]*tp[2] + d2A[ 3]*tp[3]);
    //d2v[1] = DeltaRInv * DeltaRInv * (d2A[ 6]*tp[2] + d2A[ 7]*tp[3]);
    //d2v[2] = DeltaRInv * DeltaRInv * (d2A[10]*tp[2] + d2A[11]*tp[3]);
    //d2v[3] = DeltaRInv * DeltaRInv * (d2A[14]*tp[2] + d2A[15]*tp[3]);

    //int imin=std::max(i,1);
    //int imax=std::min(i+4,NumParams+1)-1;
    //int n=imin-1, j=imin-i;
    //while(n<imax && j<4)
    //{
    //  derivs[n] = TinyVector<real_type,3>(v[j],dv[j],d2v[j]);
    //  n++; j++;
    //}
    //if(i==0) derivs[1]+= TinyVector<real_type,3>(v[0],dv[0],d2v[0]);

    return true;
  }

  inline bool evaluateDerivatives(real_type r, std::vector<real_type>& derivs)
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

  inline real_type f(real_type r)
  {
    if (r>=cutoff_radius)
      return 0.0;
    return evaluate(r);
  }
  inline real_type df(real_type r)
  {
    if (r>=cutoff_radius)
      return 0.0;
    real_type du, d2u;
    evaluate(r, du, d2u);
    return du;
  }


  bool put(xmlNodePtr cur)
  {
    ReportEngine PRE("BsplineFunctor","put(xmlNodePtr)");
    //CuspValue = -1.0e10;
    NumParams = 0;
    //cutoff_radius = 0.0;
    OhmmsAttributeSet rAttrib;
    real_type radius = -1.0;
    rAttrib.add(NumParams,   "size");
    rAttrib.add(radius,      "rcut");
    rAttrib.add(radius,      "cutoff");
    rAttrib.put(cur);
    if (radius < 0.0)
      if (periodic)
        app_log() << "  Jastrow cutoff unspecified.  Setting to Wigner-Seitz radius = " << cutoff_radius << ".\n";
      else
      {
        APP_ABORT("  Jastrow cutoff unspecified.  Cutoff must be given when using open boundary conditions");
      }
    else
      if (periodic && radius > cutoff_radius)
      {
        if (radius - cutoff_radius > 1e-4)
        {
          APP_ABORT( "  The Jastrow cutoff specified should not be larger than Wigner-Seitz radius.");
        }
        else
        {
          app_log() << "  The Jastrow cutoff specified is slightly larger than the Wigner-Seitz radius.";
          app_log() << "  Setting to Wigner-Seitz radius = " << cutoff_radius << ".\n";
        }
      }
      else
        cutoff_radius = radius;
    if (NumParams == 0)
    {
      PRE.error("You must specify a positive number of parameters for the Bspline jastrow function.",true);
    }
    app_log() << " size = " << NumParams << " parameters " << std::endl;
    app_log() << " cusp = " << CuspValue << std::endl;
    app_log() << " rcut = " << cutoff_radius << std::endl;
    resize(NumParams);
    // Now read coefficents
    xmlNodePtr xmlCoefs = cur->xmlChildrenNode;
    while (xmlCoefs != NULL)
    {
      std::string cname((const char*)xmlCoefs->name);
      if (cname == "coefficients")
      {
        std::string type("0"), id("0");
        std::string optimize("yes");
        OhmmsAttributeSet cAttrib;
        cAttrib.add(id, "id");
        cAttrib.add(type, "type");
        cAttrib.add(optimize, "optimize");
        cAttrib.put(xmlCoefs);
        if (type != "Array")
        {
          PRE.error("Unknown correlation type " + type + " in BsplineFunctor." + "Resetting to \"Array\"");
          xmlNewProp(xmlCoefs, (const xmlChar*) "type", (const xmlChar*) "Array");
        }
        std::vector<real_type> params;
        putContent(params, xmlCoefs);
        if (params.size() == NumParams)
          Parameters = params;
        else
        {
          app_log() << "Changing number of Bspline parameters from "
                    << params.size() << " to " << NumParams << ".  Performing fit:\n";
          // Fit function to new number of parameters
          const int numPoints = 500;
          BsplineFunctor<T> tmp_func(CuspValue);
          tmp_func.cutoff_radius = cutoff_radius;
          tmp_func.resize(params.size());
          tmp_func.Parameters = params;
          tmp_func.reset();
          std::vector<real_type> y(numPoints);
          Matrix<real_type> basis(numPoints,NumParams);
          std::vector<TinyVector<real_type,3> > derivs(NumParams);
          for (int i=0; i<numPoints; i++)
          {
            real_type r = (real_type)i / (real_type)numPoints * cutoff_radius;
            y[i] = tmp_func.evaluate(r);
            evaluateDerivatives(r, derivs);
            for (int j=0; j<NumParams; j++)
              basis(i,j) = derivs[j][0];
          }
          resize(NumParams);
          LinearFit(y, basis, Parameters);
          app_log() << "New parameters are:\n";
          for (int i=0; i < Parameters.size(); i++)
            app_log() << "   " << Parameters[i] << std::endl;
        }
        if(optimize == "yes")
        {
          notOpt=false;
        }
        else
        {
          notOpt=true;
        }
        for (int i=0; i< NumParams; i++)
        {
          std::stringstream sstr;
          sstr << id << "_" << i;
          myVars.insert(sstr.str(),Parameters[i],!notOpt,optimize::LOGLINEAR_P);
        }
        app_log() << "Parameter     Name      Value\n";
        myVars.print(app_log());
      }
      xmlCoefs = xmlCoefs->next;
    }
    reset();
    real_type zeros=0;
    for (int i=0; i< NumParams; i++)
      zeros+=Parameters[i]*Parameters[i];
    return zeros>1.0e-12; //true if Parameters are not zero
  }

  void initialize(int numPoints, std::vector<real_type>& x, std::vector<real_type>& y
                  , real_type cusp, real_type rcut, std::string& id, std::string& optimize )
  {
    ReportEngine PRE("BsplineFunctor","initialize");
    NumParams = numPoints;
    cutoff_radius = rcut;
    CuspValue = cusp;
    if (NumParams == 0)
    {
      PRE.error("You must specify a positive number of parameters for the Bspline jastrow function.",true);
    }
    app_log() << "Initializing BsplineFunctor from array. \n";
    app_log() << " size = " << NumParams << " parameters " << std::endl;
    app_log() << " cusp = " << CuspValue << std::endl;
    app_log() << " rcut = " << cutoff_radius << std::endl;
    resize(NumParams);
    int npts = x.size();
    Matrix<real_type> basis(npts,NumParams);
    std::vector<TinyVector<real_type,3> > derivs(NumParams);
    for (int i=0; i<npts; i++)
    {
      real_type r = x[i];
      if (r > cutoff_radius)
      {
        PRE.error("Error in BsplineFunctor::initialize: r > cutoff_radius.",true);
      }
      evaluateDerivatives(r, derivs);
      for (int j=0; j<NumParams; j++)
        basis(i,j) = derivs[j][0];
    }
    resize(NumParams);
    LinearFit(y, basis, Parameters);
    app_log() << "New parameters are:\n";
    for (int i=0; i < Parameters.size(); i++)
      app_log() << "   " << Parameters[i] << std::endl;
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
      app_log() << "Parameters of BsplineFunctor id:"
                <<id <<" are not being optimized.\n";
    }
    reset();
  }

  void reportStatus(std::ostream& os)
  {
    if (notOpt)
      return;
    myVars.print(os);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    if (notOpt)
      return;
    myVars.getIndex(active);
  }

  void checkInVariables(opt_variables_type& active)
  {
    if (notOpt)
      return;
    active.insertFrom(myVars);
  }

  void resetParameters(const opt_variables_type& active)
  {
    if (notOpt)
      return;
    for (int i=0; i<Parameters.size(); ++i)
    {
      int loc=myVars.where(i);
      if (loc>=0)
        Parameters[i]=myVars[i]=active[loc];
    }
//         if (ResetCount++ == 100)
//         {
//           ResetCount = 0;
//           if(ReportLevel) print();
//         }
    reset();
  }

  // check if this object has active optimizable parameters
  bool isOptimizable()
  {
    if (notOpt)
      return false;
    for (int i=0; i<Parameters.size(); ++i)
    {
      int loc=myVars.where(i);
      if (loc>=0)
        return true;
    }
    return false;
  }

  void print()
  {
    if(ReportLevel)
    {
      std::ofstream fout(fileName.c_str());
      this->print(fout);
    }
  }


  void print(std::ostream& os)
  {
    int n=100;
    T d=cutoff_radius/100.,r=0;
    T u,du,d2du;
    for (int i=0; i<n; ++i)
    {
      u=evaluate(r,du,d2du);
      os << std::setw(22) << r << std::setw(22) << u << std::setw(22) << du
         << std::setw(22) << d2du << std::endl;
      r+=d;
    }
  }
};

template<typename T>
inline T 
BsplineFunctor<T>::evaluateV(const int iat, const int iStart, const int iEnd,
    const T* restrict _distArray, T* restrict distArrayCompressed ) const
{
  const real_type* restrict distArray = _distArray + iStart;

  ASSUME_ALIGNED(distArrayCompressed);
  int iCount = 0;
  const int iLimit = iEnd-iStart;

#pragma vector always 
  for ( int jat = 0; jat < iLimit; jat++ ) {
    real_type r = distArray[jat];
    // pick the distances smaller than the cutoff and avoid the reference atom
    if ( r < cutoff_radius && iStart+jat != iat )
      distArrayCompressed[iCount++] = distArray[jat];
  }

  real_type d = 0.0;
#pragma omp simd reduction (+:d)
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
inline void BsplineFunctor<T>::evaluateVGL(const int iat, const int iStart, const int iEnd,
    const T* _distArray,  T* restrict _valArray, 
    T* restrict _gradArray, T* restrict _laplArray, 
    T* restrict distArrayCompressed, int* restrict distIndices ) const
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
    if ( r < cutoff_radius && iStart+jat != iat ) {
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

}
}
#endif
