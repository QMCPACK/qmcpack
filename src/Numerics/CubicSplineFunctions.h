//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//
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

/***************************************************************************
 * This class implements the orthogonal tight-binding model
 * for Si proposed by Lenosky, Kress, Kwon, Voter, Edwards, Richards
 * Yang, and Adams in PRB 55 (1997) p.1528  using cubic splines
 * to parameterize hopping terms and pair-potential.
 * Implemented by D. Trinkle July '98 on TBMD
 * Implemented by Jeongnim Kim Sep '99
 ***************************************************************************/
#ifndef CUBICSPLINEFUNCTIONS_H
#define CUBICSPLINEFUNCTIONS_H

#include <vector>
using std::vector;

template<class T>
class CubicSpline
{
public:
  int Npt;
  T Yp1, Ypn;
  vector<T> X, Y, Y2;

  CubicSpline() { }
  inline CubicSpline(const int n)
  {
    resize(n);
  }

  inline void resize(const int n)
  {
    Npt = n;
    X = vector<T>(n+1);
    Y = vector<T>(n+1);
    Y2= vector<T>(n+1);
  }
  T operator()(T x0);
  T operator()(T x0, T& yval);
  T operator()(T x0, T& yval, T& yp);
  void spline();
  void spline(const int n, T* x, T* y, const T yp1, const T ypn)
  {
    resize(n);
    X[0] = 0.0e0;
    Y[0] = 0.0e0; // not needed
    for(int i=1; i<=n; i++)
      X[i] = x[i-1];
    for(int i=1; i<=n; i++)
      Y[i] = y[i-1];
    Yp1 = yp1;
    Ypn = ypn;
    spline();
  }
};

template<class T>
class RegCubicSpline
{
  int Npt;
  T Xmin, Xmax, Dx, DxInv;
  T Yp1, Ypn, h2over6;
  vector<T> Y, Y2;
public:

  RegCubicSpline() { }
  inline RegCubicSpline(const int n)
  {
    resize(n);
  }

  inline void resize(const int n)
  {
    Npt = n;
    Y = vector<T>(n+1);
    Y2= vector<T>(n+1);
  }
  T operator()(T x0);
  T operator()(T x0, T& yval);
  T operator()(T x0, T& yval, T& yp);

  void spline();
  void spline(const int n, const T x0, const T dx, const T yp1, const T ypn, T* y)
  {
    resize(n);
    Dx = dx;
    DxInv = 1/dx;
    h2over6 = dx*dx*0.16666666666666666667e0;
    Xmin = x0;
    Xmax = static_cast<T>(n-1)*dx+x0;
    Y[0] = 0;
    for(int i=1; i<=n; i++)
      Y[i] = y[i-1];
    Yp1 = yp1;
    Ypn = ypn;
    spline();
  }
};

#include "ScaleFunctors/CubicSplineFunctions.cpp"
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * POOMA_VERSION_ID: $Id$
 ***************************************************************************/
