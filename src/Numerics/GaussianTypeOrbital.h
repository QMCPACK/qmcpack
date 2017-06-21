//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_SLATERTYPEORBITAL_H
#define QMCPLUSPLUS_SLATERTYPEORBITAL_H
#include <math.h>

/**class for Gaussian-type orbitals
 *@f[
 *\Psi_{n,l,m}({\bf R}) =N \frac{\exp{-\sigma r^2}}{r^l} \times r^l Y_{lm}(\theta,\phi)
 *@f]
 */
template<class T>
struct RadialGaussian
{
  int L;
  T Sigma;
  T Norm;
  RadialGaussian(): L(0), Sigma(1.0), Norm(1.0) { }
  RadialGaussian(int l, T sig, T norm=1.0): L(l), Sigma(sig),Norm(norm) { }

  inline void setgrid(T r) { }

  //inline void evaluate(T r, T rinv) {
  //  //Y = evaluate(r,rinv,dY,d2Y);
  //}

  inline T evaluate(T r, T rinv, T& drnl, T& d2rnl)
  {
    if(L)
    {
      T oneoverrl = pow(rinv,L);
      T rnl = Norm*exp(-Sigma*r*r)*oneoverrl;
      T x = 2.0*Sigma*r+L*rinv;
      drnl = -rnl*x;
      d2rnl = rnl*(x*x-2.0*Sigma+L*rinv*rinv);
      return rnl;
    }
    else
    {
      T rnl = Norm*exp(-Sigma*r*r);
      T x = 2.0*Sigma*r;
      drnl = -rnl*x;
      d2rnl = rnl*(x*x-2.0*Sigma);
      return rnl;
    }
  }
};

#endif
