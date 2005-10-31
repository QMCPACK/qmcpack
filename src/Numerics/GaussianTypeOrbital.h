//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_SLATERTYPEORBITAL_H
#define QMCPLUSPLUS_SLATERTYPEORBITAL_H
#include <math.h>

/**class for Gaussian-type orbitals
 *@f[
 *\Psi_{n,l,m}({\bf R}) =N \frac{\exp{-\sigma r^2}}{r^l} \times r^l Y_{lm}(\theta,\phi)
 *@f]
 */
template<class T>
struct RadialGaussian {
  int L;
  T Sigma;
  T Norm;
  RadialGaussian(): L(0), Sigma(1.0), Norm(1.0) { } 
  RadialGaussian(int l, T sig, T norm=1.0): L(l), Sigma(sig),Norm(norm) { } 

  inline void setgrid(T r) { }

  //inline void evaluate(T r, T rinv) {
  //  //Y = evaluate(r,rinv,dY,d2Y);
  //}

  inline T evaluate(T r, T rinv, T& drnl, T& d2rnl) {
    if(L) {	
      T oneoverrl = pow(rinv,L);
      T rnl = Norm*exp(-Sigma*r*r)*oneoverrl;
      T x = 2.0*Sigma*r+L*rinv;
      drnl = -rnl*x;
      d2rnl = rnl*(x*x-2.0*Sigma+L*rinv*rinv);
      return rnl;
    } else {
      T rnl = Norm*exp(-Sigma*r*r);
      T x = 2.0*Sigma*r;
      drnl = -rnl*x;
      d2rnl = rnl*(x*x-2.0*Sigma);
      return rnl;
    }
  }
};

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
