//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_BSPLINE3D_FUNCTOR_H
#define QMCPLUSPLUS_BSPLINE3D_FUNCTOR_H
#include "Numerics/OptimizableFunctorBase.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "Numerics/LinearFit.h"
#include <cstdio>

namespace qmcplusplus
{

struct BsplineFunctor3D: public OptimizableFunctorBase
{

  typedef real_type value_type;
  int NumParams_eI, NumParams_ee;
  int Dummy;
  const TinyVector<real_type,16> A, dA, d2A, d3A;
  //static const real_type A[16], dA[16], d2A[16];
  real_type Deltax, DeltaxInv;
  real_type DeltaR_eI, DeltaRInv_eI;
  real_type Y, dY, d2Y;
  Array<real_type,3> SplineCoefs;
  // Stores the derivatives w.r.t. SplineCoefs
  // of the u, du/dr, and d2u/dr2
  std::vector<TinyVector<real_type,3> > SplineDerivs;
  std::vector<real_type> Parameters;
  Array<real_type,3> ParamArray;
  std::vector<std::string> ParameterNames;
  std::string iSpecies, eSpecies1, eSpecies2;
  int ResetCount;

  ///constructor
  BsplineFunctor3D(real_type ecusp=0.0, real_type icusp=0.0) :
    NumParams_eI(0), NumParams_ee(0),
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
    ResetCount(0)
  {
    cutoff_radius = 0.0;
  }

  OptimizableFunctorBase* makeClone() const
  {
    return new BsplineFunctor3D(*this);
  }

  void resize(int neI, int nee)
  {
    NumParams_eI = neI;
    NumParams_ee = nee;
    int numCoefs_eI = NumParams_eI + 4;
    int numCoefs_ee = NumParams_ee + 4;
    SplineCoefs.resize(numCoefs_ee, numCoefs_eI, numCoefs_eI);
    ParamArray.resize(NumParams_ee, NumParams_eI, NumParams_eI);
    int numParams = NumParams_eI*(NumParams_eI+1)/2 * NumParams_ee;
    Parameters.resize(numParams);
    int numKnots_eI = numCoefs_eI - 2;
    DeltaR_eI = 0.5*cutoff_radius / (real_type)(numKnots_eI - 1);
    DeltaRInv_eI = 1.0/DeltaR_eI;
    int numKnots_ee = numCoefs_ee - 2;
    Deltax = 1.0 / (real_type)(numKnots_ee - 1);
    DeltaxInv = 1.0/Deltax;
  }

  inline int getNumParameters()
  {
    return Parameters.size();
  }

  void reset()
  {
    int numCoefs_eI = NumParams_eI + 4;
    int numKnots_eI =  numCoefs_eI - 2;
    int numCoefs_ee = NumParams_ee + 4;
    int numKnots_ee =  numCoefs_ee - 2;
    DeltaR_eI = 0.5*cutoff_radius / (real_type)(numKnots_eI - 1);
    Deltax = 1.0 / (real_type)(numKnots_ee - 1);
    DeltaRInv_eI = 1.0/DeltaR_eI;
    DeltaxInv = 1.0/Deltax;
    // Zero out all coefficients
    for (int i=0; i<SplineCoefs.size(0); i++)
      for (int j=0; j<SplineCoefs.size(1); j++)
        for (int k=0; k<SplineCoefs.size(2); k++)
          SplineCoefs(i,j,k) = 0.0;
    // Set unconstrained coefficients
    for (int i=2; i<NumParams_ee; i++)
      for (int j=2; j<NumParams_eI; j++)
        for (int k=2; k<NumParams_eI; k++)
          SplineCoefs(i+1,j+1,k+1) = ParamArray(i,j,k);
    // j-k plane
    // Set e-I cusp parameters
    for (int j=2; j<NumParams_eI; j++)
      for (int k=2; k<NumParams_eI; k++)
      {
        SplineCoefs(1,j+1,k+1) = ParamArray(0,j,k);
        SplineCoefs(2,j+1,k+1) = ParamArray(1,j,k);
        SplineCoefs(0,k+1,j+1) = ParamArray(1,j,k);
      }
    // i-j plane
    // Set e-e cusp parameters
    for (int i=2; i<NumParams_ee; i++)
      for (int j=2; j<NumParams_eI; j++)
      {
        SplineCoefs(i+1,j+1,1) = ParamArray(i,j,0);
        SplineCoefs(i+1,j+1,2) = ParamArray(i,j,1);
        SplineCoefs(i+1,j+1,0) = ParamArray(i,j,1);
      }
    // i-k plane
    // Set e-e cusp parameters
    for (int i=2; i<NumParams_ee; i++)
      for (int k=2; k<NumParams_eI; k++)
      {
        SplineCoefs(i+1,1,k+1) = ParamArray(i,0,k);
        SplineCoefs(i+1,2,k+1) = ParamArray(i,1,k);
        SplineCoefs(i+1,0,k+1) = ParamArray(i,1,k);
      }
    // i edge
    for (int i=2; i<NumParams_ee; i++)
    {
      SplineCoefs(i+1,1,1) = ParamArray(i,0,0);
      SplineCoefs(i+1,2,1) = ParamArray(i,1,0);
      SplineCoefs(i+1,0,1) = ParamArray(i,1,0);
      SplineCoefs(i+1,1,2) = ParamArray(i,0,1);
      SplineCoefs(i+1,2,2) = ParamArray(i,1,1);
      SplineCoefs(i+1,0,2) = ParamArray(i,1,1);
      SplineCoefs(i+1,1,0) = ParamArray(i,0,1);
      SplineCoefs(i+1,2,0) = ParamArray(i,1,1);
      SplineCoefs(i+1,0,0) = ParamArray(i,1,1);
    }
    // j edge
    for (int j=2; j<NumParams_eI; j++)
    {
      SplineCoefs(1,j+1,1) = ParamArray(0,j,0);
      SplineCoefs(2,j+1,1) = ParamArray(1,j,0);
      SplineCoefs(0,j+1,1) = ParamArray(1,j,0);
      SplineCoefs(1,j+1,2) = ParamArray(0,j,1);
      SplineCoefs(2,j+1,2) = ParamArray(1,j,1);
      SplineCoefs(0,j+1,2) = ParamArray(1,j,1);
      SplineCoefs(1,j+1,0) = ParamArray(0,j,1);
      SplineCoefs(2,j+1,0) = ParamArray(1,j,1);
      SplineCoefs(0,j+1,0) = ParamArray(1,j,1);
    }
    // k edge
    for (int k=2; k<NumParams_eI; k++)
    {
      SplineCoefs(1,1,k+1) = ParamArray(0,0,k);
      SplineCoefs(2,1,k+1) = ParamArray(1,0,k);
      SplineCoefs(0,1,k+1) = ParamArray(1,0,k);
      SplineCoefs(1,2,k+1) = ParamArray(0,1,k);
      SplineCoefs(2,2,k+1) = ParamArray(1,1,k);
      SplineCoefs(0,2,k+1) = ParamArray(1,1,k);
      SplineCoefs(1,0,k+1) = ParamArray(0,1,k);
      SplineCoefs(2,0,k+1) = ParamArray(1,1,k);
      SplineCoefs(0,0,k+1) = ParamArray(1,1,k);
    }
    // Copy the 8 uniquely determined values
    SplineCoefs(1,1,1) = ParamArray(0,0,0);
    SplineCoefs(1,1,2) = ParamArray(0,0,1);
    SplineCoefs(1,2,1) = ParamArray(0,1,0);
    SplineCoefs(1,2,2) = ParamArray(0,1,1);
    SplineCoefs(2,1,1) = ParamArray(1,0,0);
    SplineCoefs(2,1,2) = ParamArray(1,0,1);
    SplineCoefs(2,2,1) = ParamArray(1,1,0);
    SplineCoefs(2,2,2) = ParamArray(1,1,1);
    // Now satisfy cusp constraints
    // ee
    SplineCoefs(1,1,0) = ParamArray(0,0,1);
    SplineCoefs(1,2,0) = ParamArray(0,1,1);
    SplineCoefs(2,1,0) = ParamArray(1,0,1);
    SplineCoefs(2,2,0) = ParamArray(1,1,1);
    SplineCoefs(1,0,1) = ParamArray(0,1,0);
    SplineCoefs(1,0,2) = ParamArray(0,1,1);
    SplineCoefs(2,0,1) = ParamArray(1,1,0);
    SplineCoefs(2,0,2) = ParamArray(1,1,1);
    // eI
    SplineCoefs(0,1,1) = ParamArray(1,0,0);
    SplineCoefs(0,1,2) = ParamArray(1,0,1);
    SplineCoefs(0,2,1) = ParamArray(1,1,0);
    SplineCoefs(0,2,2) = ParamArray(1,1,1);
    // More than one cusp constraint
    SplineCoefs(0,0,1) = ParamArray(1,1,0);
    SplineCoefs(0,1,0) = ParamArray(1,0,1);
    SplineCoefs(1,0,0) = ParamArray(0,1,1);
    SplineCoefs(0,0,2) = ParamArray(1,1,1);
    SplineCoefs(0,2,0) = ParamArray(1,1,1);
    SplineCoefs(2,0,0) = ParamArray(1,1,1);
    SplineCoefs(0,0,0) = ParamArray(1,1,1);
  }

  inline real_type evaluate(real_type r_12,
                            real_type r_1I,
                            real_type r_2I) const
  {
    if (r_12 >= cutoff_radius || r_1I >= 0.5*cutoff_radius ||
        r_2I >= 0.5*cutoff_radius)
      return 0.0;
    real_type x = r_12 / (r_1I + r_2I);
    x    *= DeltaxInv;
    r_1I *= DeltaRInv_eI;
    r_2I *= DeltaRInv_eI;
    real_type ipart, t, u, v;
    int i, j, k;
    t = std::modf(x, &ipart);
    i = (int) ipart;
    u = std::modf(r_1I, &ipart);
    j = (int) ipart;
    v = std::modf(r_2I, &ipart);
    k = (int) ipart;
    real_type tp[4], up[4], vp[4], a[4], b[4], c[4];
    tp[0] = t*t*t;
    tp[1] = t*t;
    tp[2] = t;
    tp[3] = 1.0;
    up[0] = u*u*u;
    up[1] = u*u;
    up[2] = u;
    up[3] = 1.0;
    vp[0] = v*v*v;
    vp[1] = v*v;
    vp[2] = v;
    vp[3] = 1.0;
    int index=0;
    for (int m=0; m<4; m++)
    {
      a[m] = b[m] = c[m] = 0.0;
      for (int n=0; n<4; n++)
      {
        a[m] += A[index] * tp[n];
        b[m] += A[index] * up[n];
        c[m] += A[index] * vp[n];
        index++;
      }
    }
    real_type val = 0.0;
    for (int ia=0; ia<4; ia++)
      for (int ib=0; ib<4; ib++)
        for (int ic=0; ic<4; ic++)
          val += SplineCoefs(i+ia, j+ib, k+ic)*a[ia]*b[ib]*c[ic];
    return val;
  }

  inline real_type evaluateV(int Nptcl,
                             const real_type* restrict r_12_array,
                             const real_type* restrict r_1I_array,
                             const real_type* restrict r_2I_array) const
  {
    real_type val_tot(0);
    for(int ptcl=0; ptcl<Nptcl; ptcl++)
      val_tot+=evaluate(r_12_array[ptcl],r_1I_array[ptcl],r_2I_array[ptcl]);
    return val_tot;
  }

  inline real_type evaluate(real_type r_12, real_type r_1I, real_type r_2I,
                            TinyVector<real_type,3> &grad,
                            Tensor<real_type,3> &hess) const
  {
    if (r_12 >= cutoff_radius || r_1I >= 0.5*cutoff_radius ||
        r_2I >= 0.5*cutoff_radius)
    {
      grad = 0.0;
      hess = 0.0;
      return 0.0;
    }
    // double eps = 1.0e-6;
    // grad[0] = (evaluate (r_12+eps, r_1I, r_2I) -evaluate (r_12-eps, r_1I, r_2I))/(2.0*eps);
    // grad[1] = (evaluate (r_12, r_1I+eps, r_2I) -evaluate (r_12, r_1I-eps, r_2I))/(2.0*eps);
    // grad[2] = (evaluate (r_12, r_1I, r_2I+eps) -evaluate (r_12, r_1I, r_2I-eps))/(2.0*eps);
    real_type qInv = 1.0/(r_1I + r_2I);
    real_type x = r_12 * qInv;
    x    *= DeltaxInv;
    r_1I *= DeltaRInv_eI;
    r_2I *= DeltaRInv_eI;
    real_type ipart, t, u, v;
    int i, j, k;
    t = std::modf(x, &ipart);
    i = (int) ipart;
    u = std::modf(r_1I, &ipart);
    j = (int) ipart;
    v = std::modf(r_2I, &ipart);
    k = (int) ipart;
    real_type tp[4], up[4], vp[4], a[4], b[4], c[4],
              da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
    tp[0] = t*t*t;
    tp[1] = t*t;
    tp[2] = t;
    tp[3] = 1.0;
    up[0] = u*u*u;
    up[1] = u*u;
    up[2] = u;
    up[3] = 1.0;
    vp[0] = v*v*v;
    vp[1] = v*v;
    vp[2] = v;
    vp[3] = 1.0;
    int index=0;
    for (int m=0; m<4; m++)
    {
      a[m]=b[m]=c[m]=da[m]=db[m]=dc[m]=d2a[m]=d2b[m]=d2c[m]=0.0;
      for (int n=0; n<4; n++)
      {
        a[m]+=A[index]*tp[n];
        da[m]+=dA[index]*tp[n];
        d2a[m]+=d2A[index]*tp[n];
        b[m]+=A[index]*up[n];
        db[m]+=dA[index]*up[n];
        d2b[m]+=d2A[index]*up[n];
        c[m]+=A[index]*vp[n];
        dc[m]+=dA[index]*vp[n];
        d2c[m]+=d2A[index]*vp[n];
        index++;
      }
    }
    real_type val = 0.0;
    TinyVector<real_type,3> gradF;
    Tensor<real_type,3> hessF;
    gradF = 0.0;
    hessF = 0.0;
    for (int ia=0; ia<4; ia++)
      for (int ib=0; ib<4; ib++)
        for (int ic=0; ic<4; ic++)
        {
          real_type coef = SplineCoefs(i+ia, j+ib, k+ic);
          val        += coef *  a[ia] *  b[ib] *  c[ic];
          gradF[0]   += coef * da[ia] *  b[ib] *  c[ic];
          gradF[1]   += coef *  a[ia] * db[ib] *  c[ic];
          gradF[2]   += coef *  a[ia] *  b[ib] * dc[ic];
          hessF(0,0) += coef *d2a[ia] *  b[ib] *  c[ic];
          hessF(0,1) += coef * da[ia] * db[ib] *  c[ic];
          hessF(0,2) += coef * da[ia] *  b[ib] * dc[ic];
          hessF(1,1) += coef *  a[ia] *d2b[ib] *  c[ic];
          hessF(1,2) += coef *  a[ia] * db[ib] * dc[ic];
          hessF(2,2) += coef *  a[ia] *  b[ib] *d2c[ic];
        }
    gradF[0] *= DeltaxInv;
    gradF[1] *= DeltaRInv_eI;
    gradF[2] *= DeltaRInv_eI;
    grad[0] = qInv*gradF[0];
    grad[1] = gradF[1] - r_12*qInv*qInv*gradF[0];
    grad[2] = gradF[2] - r_12*qInv*qInv*gradF[0];
    hessF(0,0) *= DeltaxInv * DeltaxInv;
    hessF(0,1) *= DeltaxInv * DeltaRInv_eI;
    hessF(0,2) *= DeltaxInv * DeltaRInv_eI;
    hessF(1,1) *= DeltaRInv_eI * DeltaRInv_eI;
    hessF(1,2) *= DeltaRInv_eI * DeltaRInv_eI;
    hessF(2,2) *= DeltaRInv_eI * DeltaRInv_eI;
    hess(0,0) = qInv*qInv*hessF(0,0);
    hess(0,1) = qInv*hessF(0,1) - gradF[0]*qInv*qInv - hessF(0,0)*r_12*qInv*qInv*qInv;
    hess(0,2) = qInv*hessF(0,2) - gradF[0]*qInv*qInv - hessF(0,0)*r_12*qInv*qInv*qInv;
    hess(1,1) = hessF(1,1) + 2.0*gradF[0]*r_12 * qInv*qInv*qInv -         2.0*hessF(0,1) *r_12*qInv*qInv + hessF(0,0)*r_12*r_12*qInv*qInv*qInv*qInv;
    hess(1,2) = hessF(1,2) + 2.0*gradF[0]*r_12 * qInv*qInv*qInv - (hessF(0,2)+hessF(0,1))*r_12*qInv*qInv + hessF(0,0)*r_12*r_12*qInv*qInv*qInv*qInv;
    hess(2,2) = hessF(2,2) + 2.0*gradF[0]*r_12 * qInv*qInv*qInv -         2.0*hessF(0,2) *r_12*qInv*qInv + hessF(0,0)*r_12*r_12*qInv*qInv*qInv*qInv;
    hess(1,0) = hess(0,1);
    hess(2,0) = hess(0,2);
    hess(2,1) = hess(1,2);
    return val;
  }

  // assume r_1I < L && r_2I < L, compression and screening is handled outside
  inline void evaluateVGL(int Nptcl, const real_type* restrict r_12_array,
                          const real_type* restrict r_1I_array,
                          const real_type* restrict r_2I_array,
                          real_type* restrict val_array,
                          real_type* restrict grad0_array,
                          real_type* restrict grad1_array,
                          real_type* restrict grad2_array,
                          real_type* restrict hess00_array,
                          real_type* restrict hess11_array,
                          real_type* restrict hess22_array,
                          real_type* restrict hess01_array,
                          real_type* restrict hess02_array) const
  {
    APP_ABORT("BsplineFunctor3D::evaluateVGL not implemented yet!");
  }

  inline real_type evaluate(real_type r_12, real_type r_1I, real_type r_2I,
                            TinyVector<real_type,3> &grad,
                            Tensor<real_type,3> &hess,
                            TinyVector<Tensor<real_type,3>,3> &d3)
  {
    return 0.0;
  }


  inline real_type evaluate(real_type r, real_type rinv)
  {
    return 0.0;
  }


  inline bool
  evaluateDerivatives(real_type r, std::vector<TinyVector<real_type,3> >& derivs)
  {
    //what is this?
    return false;
  }

  inline bool
  evaluateDerivatives (real_type r_12, real_type r_1I, real_type r_2I,
                       std::vector<real_type> &d_vals,
                       std::vector<TinyVector<real_type,3> >& d_grads,
                       std::vector<Tensor<real_type,3> > &d_hess)
  {
    return false;
  }

  inline real_type f(real_type r)
  {
    return 0.0;
  }
  inline real_type df(real_type r)
  {
    return 0.0;
  }

  bool put(xmlNodePtr cur)
  {
    ReportEngine PRE("BsplineFunctor3D","put(xmlNodePtr)");
    //CuspValue = -1.0e10;
    NumParams_eI = NumParams_ee = 0;
    cutoff_radius = 0.0;
    OhmmsAttributeSet rAttrib;
    rAttrib.add(NumParams_ee,   "esize");
    rAttrib.add(NumParams_eI,   "isize");
    rAttrib.add(cutoff_radius,  "rcut");
    rAttrib.put(cur);
    if (NumParams_eI == 0)
      PRE.error("You must specify a positive number for \"isize\"",true);
    if (NumParams_ee == 0)
      PRE.error("You must specify a positive number for \"esize\"",true);
    app_log() << " esize = " << NumParams_ee << " parameters " << std::endl;
    app_log() << " isize = " << NumParams_eI << " parameters " << std::endl;
    app_log() << " rcut = " << cutoff_radius << std::endl;
    resize(NumParams_eI, NumParams_ee);
    // Now read coefficents
    xmlNodePtr xmlCoefs = cur->xmlChildrenNode;
    while (xmlCoefs != NULL)
    {
      std::string cname((const char*)xmlCoefs->name);
      if (cname == "coefficients")
      {
        std::string type("0"), id("0");
        OhmmsAttributeSet cAttrib;
        cAttrib.add(id, "id");
        cAttrib.add(type, "type");
        cAttrib.put(xmlCoefs);
        if (type != "Array")
        {
          PRE.error("Unknown correlation type " + type +
                    " in BsplineFunctor3D." + "Resetting to \"Array\"");
          xmlNewProp(xmlCoefs, (const xmlChar*) "type", (const xmlChar*) "Array");
        }
        std::vector<real_type> params;
        putContent(params, xmlCoefs);
        if (params.size() == Parameters.size())
          Parameters = params;
        else
          if (params.size() == 0)
          {
            app_log()<<" Initializing all parameters to zero"<< std::endl;
          }
          else
          {
            app_error() << "Expected " << Parameters.size() << " parameters,"
                        << " but found only " << params.size()
                        << " in BsplineFunctor3D.\n";
            abort();
          }
        // Setup parameter names
        int index=0;
        for (int i=0; i< NumParams_ee; i++)
          for (int j=0; j < NumParams_eI; j++)
            for (int k=0; k<=j; k++)
            {
              std::stringstream sstr;
              sstr << id << "_" << i << "_" << j << "_" << k;
              myVars.insert(sstr.str(),Parameters[index],true,optimize::LOGLINEAR_P);
              ParamArray(i,j,k) = ParamArray(i,k,j) = Parameters[index];
              index++;
            }
        app_log() << "Parameter     Name      Value\n";
        myVars.print(app_log());
      }
      xmlCoefs = xmlCoefs->next;
    }
    reset();
    print();
    return true;
  }

  void checkInVariables(opt_variables_type& active)
  {
    active.insertFrom(myVars);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active);
  }

  void resetParameters(const opt_variables_type& active)
  {
    int iparam = 0;
    for (int i=0; i<NumParams_ee; i++)
      for (int j=0; j<NumParams_eI; j++)
        for (int k=0; k<=j; k++)
        {
          int loc = myVars.where(iparam);
          if (loc >=0)
            Parameters[iparam] = myVars[iparam] = active[loc];
          ParamArray(i,j,k) = Parameters[iparam];
          ParamArray(i,k,j) = Parameters[iparam];
          iparam++;
        }
    reset();
    // for(int i=0; i<Parameters.size(); ++i) {
    //   int loc=myVars.where(i);
    //   if(loc>=0) Parameters[i]=myVars[i]=active[loc];
    // }
    if (ResetCount++ == 100)
    {
      ResetCount = 0;
      //print();
    }
    reset();
  }


  void print()
  {
    const int N = 100;
    std::string fname = iSpecies + ".J3.h5";
    hid_t hid = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    Array<real_type,3> val(N,N,N);
    for (int i=0; i<N; i++)
    {
      double r_12 = (real_type)i/(real_type)(N-1);
      for (int j=0; j<N; j++)
      {
        double r_1I = (real_type)j/(real_type)(N-1) * 0.5*cutoff_radius;
        for (int k=0; k<N; k++)
        {
          double r_2I = (real_type)k/(real_type)(N-1) * 0.5*cutoff_radius;
          val(i,j,k) = evaluate(r_12*(r_1I+r_2I), r_1I, r_2I);
        }
      }
    }
    Array<double,3> SplineCoefs_DP(NumParams_ee + 4, NumParams_eI + 4, NumParams_eI + 4);
    Array<double,3> ParamArray_DP(NumParams_ee, NumParams_eI, NumParams_eI);
    Array<double,3> val_DP(N,N,N);
    HDFAttribIO<Array<double,3> > coefs_attrib(SplineCoefs_DP);
    HDFAttribIO<Array<double,3> > param_attrib(ParamArray_DP);
    HDFAttribIO<Array<double,3> > val_attrib(val_DP);
    SplineCoefs_DP = SplineCoefs;
    ParamArray_DP = ParamArray;
    val_DP = val;
    val_attrib.write(hid, "val");
    coefs_attrib.write(hid, "coefs");
    param_attrib.write(hid, "params");
    H5Fclose(hid);
    // std::string fname = (elementType != "") ? elementType : pairType;
    // fname = fname + ".dat";
    // //cerr << "Writing " << fname << " file.\n";
    // FILE *fout = fopen (fname.c_str(), "w");
    // for (double r=0.0; r<cutoff_radius; r+=0.001)
    //  fprintf (fout, "%8.3f %16.10f\n", r, evaluate(r));
    // fclose(fout);
  }


  void print(std::ostream& os)
  {
    /* no longer correct. Ye Luo
    int n=100;
    real_type d=cutoff_radius/100.,r=0;
    real_type u,du,d2du;
    for (int i=0; i<n; ++i)
    {
      u=evaluate(r,du,d2du);
      os << std::setw(22) << r << std::setw(22) << u << std::setw(22) << du
         << std::setw(22) << d2du << std::endl;
      r+=d;
    }
    */
  }
};
}
#endif
