//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_QUADRATURE_H
#define QMCPLUSPLUS_QUADRATURE_H

#include <assert.h>
#include "Numerics/Ylm.h"
#include "type_traits/scalar_traits.h"
#include "QMCWaveFunctions/lcao/SoaSphericalTensor.h"

namespace qmcplusplus
{

template<class T>
struct Quadrature3D
{
  typedef T                    RealType;
  typedef TinyVector<T,3>      PosType;
  typedef OHMMS_PRECISION_FULL mRealType;

  int nk;
  typedef enum {SINGLE, TETRA, OCTA, ICOSA} SymmType;
  SymmType symmetry;
  int lexact;
  RealType A, B, C, D;
  std::vector<PosType>  xyz_m;
  std::vector<RealType> weight_m;
  bool quad_ok;
  const bool fail_abort;

  Quadrature3D(int rule, bool request_abort=true): quad_ok(true), fail_abort(request_abort)
  {
    A = B = C = D = 0;
    switch (rule)
    {
    case 1:
      nk = 1;
      symmetry = SINGLE;
      lexact = 0;
      A = 1.0;
      break;
    case 2:
      nk = 4;
      symmetry = TETRA;
      lexact = 2;
      A=0.25;
      break;
    case 3:
      nk = 6;
      symmetry = OCTA;
      lexact = 3;
      A=1.0/6.0;
      break;
    case 4:
      nk = 12;
      symmetry = ICOSA;
      lexact = 5;
      A = 1.0/12.0;
      B = 1.0/12.0;
      break;
    case 5:
      nk = 18;
      symmetry = OCTA;
      lexact = 5;
      A = 1.0/30.0;
      B = 1.0/15.0;
      break;
    case 6:
      nk = 26;
      symmetry = OCTA;
      lexact = 7;
      A = 1.0  / 21.0;
      B = 4.0  / 105.0;
      C = 27.0 / 840.0;
      break;
    case 7:
      nk = 50;
      symmetry = OCTA;
      lexact = 11;
      A = 4.0/315.0;
      B = 64.0/2835.0;
      C = 27.0/1280.0;
      D = 14641.0/725760.0;
      break;
    default:
      ERRORMSG("Unrecognized spherical quadrature rule " << rule << ".");
      abort();
    }
    // First, build a_i, b_i, and c_i points
    std::vector<PosType> a, b, c, d;
    RealType p = 1.0/std::sqrt(2.0);
    RealType q = 1.0/std::sqrt(3.0);
    RealType r = 1.0/std::sqrt(11.0);
    RealType s = 3.0/std::sqrt(11.0);
    if (symmetry == SINGLE)
    {
      a.push_back (PosType(1.0, 0.0, 0.0));
    }
    else if (symmetry == TETRA)
    {
      a.push_back(PosType( q, q, q));
      a.push_back(PosType( q,-q,-q));
      a.push_back(PosType(-q, q,-q));
      a.push_back(PosType(-q,-q, q));
    }
    else if (symmetry == OCTA)
    {
      a.push_back(PosType( 1.0, 0.0, 0.0));
      a.push_back(PosType(-1.0, 0.0, 0.0));
      a.push_back(PosType( 0.0, 1.0, 0.0));
      a.push_back(PosType( 0.0,-1.0, 0.0));
      a.push_back(PosType( 0.0, 0.0, 1.0));
      a.push_back(PosType( 0.0, 0.0,-1.0));
      b.push_back(PosType(   p,   p, 0.0));
      b.push_back(PosType(   p,  -p, 0.0));
      b.push_back(PosType(  -p,   p, 0.0));
      b.push_back(PosType(  -p,  -p, 0.0));
      b.push_back(PosType(   p, 0.0,   p));
      b.push_back(PosType(   p, 0.0,  -p));
      b.push_back(PosType(  -p, 0.0,   p));
      b.push_back(PosType(  -p, 0.0,  -p));
      b.push_back(PosType( 0.0,   p,   p));
      b.push_back(PosType( 0.0,   p,  -p));
      b.push_back(PosType( 0.0,  -p,   p));
      b.push_back(PosType( 0.0,  -p,  -p));
      c.push_back(PosType(   q,   q,   q));
      c.push_back(PosType(   q,   q,  -q));
      c.push_back(PosType(   q,  -q,   q));
      c.push_back(PosType(   q,  -q,  -q));
      c.push_back(PosType(  -q,   q,   q));
      c.push_back(PosType(  -q,   q,  -q));
      c.push_back(PosType(  -q,  -q,   q));
      c.push_back(PosType(  -q,  -q,  -q));
      d.push_back(PosType(   r,   r,   s));
      d.push_back(PosType(   r,   r,  -s));
      d.push_back(PosType(   r,  -r,   s));
      d.push_back(PosType(   r,  -r,  -s));
      d.push_back(PosType(  -r,   r,   s));
      d.push_back(PosType(  -r,   r,  -s));
      d.push_back(PosType(  -r,  -r,   s));
      d.push_back(PosType(  -r,  -r,  -s));
      d.push_back(PosType(   r,   s,   r));
      d.push_back(PosType(   r,   s,  -r));
      d.push_back(PosType(   r,  -s,   r));
      d.push_back(PosType(   r,  -s,  -r));
      d.push_back(PosType(  -r,   s,   r));
      d.push_back(PosType(  -r,   s,  -r));
      d.push_back(PosType(  -r,  -s,   r));
      d.push_back(PosType(  -r,  -s,  -r));
      d.push_back(PosType(   s,   r,   r));
      d.push_back(PosType(   s,   r,  -r));
      d.push_back(PosType(   s,  -r,   r));
      d.push_back(PosType(   s,  -r,  -r));
      d.push_back(PosType(  -s,   r,   r));
      d.push_back(PosType(  -s,   r,  -r));
      d.push_back(PosType(  -s,  -r,   r));
      d.push_back(PosType(  -s,  -r,  -r));
    }
    else if (symmetry == ICOSA)
    {
      mRealType t, p;  // theta and phi
      // a points
      t = 0.0;
      p=0.0;
      a.push_back(PosType(std::cos(t),std::sin(t)*std::cos(p),std::sin(t)*std::sin(p)));
      t = M_PI;
      p=0.0;
      a.push_back(PosType (std::cos(t),std::sin(t)*std::cos(p),std::sin(t)*std::sin(p)));
      // b points
      for (int k=0; k<5; k++)
      {
        t = std::atan(2.0);
        p = (mRealType)(2*k+0)*M_PI/5.0;
        b.push_back(PosType (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
        t = M_PI-std::atan(2.0);
        p = (mRealType)(2*k+1)*M_PI/5.0;
        b.push_back(PosType (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
      }
      // c points
      mRealType t1 = std::acos ((2.0+std::sqrt(5.0)) / std::sqrt(15.0+6.0*std::sqrt(5.0)));
      mRealType t2 = std::acos (      1.0            / std::sqrt(15.0+6.0*std::sqrt(5.0)));
      for (int k=0; k<5; k++)
      {
        t = t1;
        p = (mRealType)(2*k+1)*M_PI/5.0;
        c.push_back(PosType (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
        t = t2;
        p = (mRealType)(2*k+1)*M_PI/5.0;
        c.push_back(PosType (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
        t = M_PI - t1;
        p = (mRealType)(2*k+0)*M_PI/5.0;
        c.push_back(PosType (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
        t = M_PI - t2;
        p = (mRealType)(2*k+0)*M_PI/5.0;
        c.push_back(PosType (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
      }
    }
    // Now, construct rule
    if (std::abs(A) > 1.0e-10)
      for (int i=0; i<a.size(); i++)
      {
        xyz_m.push_back(a[i]);
        weight_m.push_back(A);
      }
    if (std::abs(B) > 1.0e-10)
      for (int i=0; i<b.size(); i++)
      {
        xyz_m.push_back(b[i]);
        weight_m.push_back(B);
      }
    if (std::abs(C) > 1.0e-10)
      for (int i=0; i<c.size(); i++)
      {
        xyz_m.push_back(c[i]);
        weight_m.push_back(C);
      }
    if (std::abs(D) > 1.0e-10)
      for (int i=0; i<d.size(); i++)
      {
        xyz_m.push_back(d[i]);
        weight_m.push_back(D);
      }
    // Finally, check the rule for correctness
    assert (xyz_m.size() == nk);
    assert (weight_m.size() == nk);
    double wSum = 0.0;
    const RealType delta=2*std::numeric_limits<float>::epsilon();
    for (int k=0; k < nk; k++)
    {
      PosType r = xyz_m[k];
      double nrm = dot(r,r);
      assert (std::abs(nrm-1.0) < delta);
      wSum += weight_m[k];
      //cout << pp_nonloc->xyz_m[k] << " " << pp_nonloc->weight_m[k] << std::endl;
    }
    assert (std::abs(wSum - 1.0) < delta);
    // Check the quadrature rule
    // using complex spherical harmonics
    CheckQuadratureRule(lexact);
    // using real spherical harmonics
    CheckQuadratureRuleReal(lexact);
  }

  void CheckQuadratureRule(int lexact)
  {
    std::vector<PosType> &grid = xyz_m;
    std::vector<RealType> &w = weight_m;
    for (int l1=0; l1<=lexact; l1++)
      for (int l2=0; l2 <= (lexact-l1); l2++)
        for (int m1=-l1; m1<=l1; m1++)
          for (int m2=-l2; m2<=l2; m2++)
          {
            std::complex<mRealType> sum(0.0, 0.0);
            for (int k=0; k<grid.size(); k++)
            {
              std::complex<mRealType> v1 = Ylm(l1, m1, grid[k]);
              std::complex<mRealType> v2 = Ylm(l2, m2, grid[k]);
              sum += 4.0*M_PI*w[k] * qmcplusplus::conj(v1)*v2;
            }
            mRealType re = real (sum);
            mRealType im = imag (sum);
            if ((l1==l2) && (m1==m2))
              re -= 1.0;
            if ((std::abs(im) > 7*std::numeric_limits<float>::epsilon()) || (std::abs(re) > 7*std::numeric_limits<float>::epsilon()))
            {
              app_error() << "Broken spherical quadrature for " << grid.size() << "-point rule.\n" << std::endl;
              app_error() << "  Should be zero:  Real part = " << re << " Imaginary part = " << im << std::endl;
              quad_ok = false;
              if(fail_abort) APP_ABORT("Give up");
            }
//   	    fprintf (stderr, "(l1,m1,l2m,m2) = (%2d,%2d,%2d,%2d)  sum = (%20.16f %20.16f)\n",
//   	     l1, m1, l2, m2, real(sum), imag(sum));
          }
  }

  void CheckQuadratureRuleReal(int lexact)
  {
    std::vector<PosType> &grid = xyz_m;
    std::vector<RealType> &w = weight_m;
    SoaSphericalTensor<RealType> Ylm(lexact);
    const RealType* restrict Ylm_v=Ylm[0];
    for (int l1=0; l1<=lexact; l1++)
      for (int l2=0; l2 <= (lexact-l1); l2++)
        for (int m1=-l1; m1<=l1; m1++)
          for (int m2=-l2; m2<=l2; m2++)
          {
            mRealType sum(0.0);
            for (int k=0; k<grid.size(); k++)
            {
              Ylm.evaluateV(grid[k][0],grid[k][1],grid[k][2]);
              RealType v1 = Ylm_v[Ylm.index(l1, m1)];
              RealType v2 = Ylm_v[Ylm.index(l2, m2)];
              sum += 4.0*M_PI*w[k] * v1*v2;
            }
            if ((l1==l2) && (m1==m2))
              sum -= 1.0;
            if (std::abs(sum) > 15*std::numeric_limits<float>::epsilon())
            {
              app_error() << "Broken real spherical quadrature for " << grid.size() << "-point rule.\n" << std::endl;
              app_error() << "  Should be zero:  " << sum << std::endl;
              quad_ok = false;
              if(fail_abort) APP_ABORT("Give up");
            }
          }
  }
};

} // namespace qmcPlusPlus

#endif
