//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim and Simone Chiesa
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file OneDimNumGridFunctor.h
 * @brief Definition of OneDimNumGridFunctor
 */
#ifndef QMCPLUSPLUS_ONEDIMNUMGRIDFUNCTOR_H
#define QMCPLUSPLUS_ONEDIMNUMGRIDFUNCTOR_H
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimCubicSpline.h"

namespace qmcplusplus
{

/** adaptor class to handle a temporary OneDimCubicSpline on a numerical grid.
 */
template<class T>
struct OneDimNumGridFunctor
{

  typedef NumericalGrid<T> GridType;
  typedef OneDimCubicSpline<T> FuncType;

  GridType myGrid;
  FuncType myFunc;

  inline OneDimNumGridFunctor()
  {
    myFunc.m_grid=&myGrid;
  }

  inline T r(int i) const
  {
    return myGrid[i];
  }
  inline T operator()(int i) const
  {
    return myFunc(i);
  }

  inline T rmax() const
  {
    return myGrid.rmax();
  }

  inline T rmin() const
  {
    return myGrid.rmin();
  }

  inline T splint(T r)
  {
    return myFunc.splint(r);
  }

  void put(int npoints, std::istream& fin)
  {
    myGrid.resize(npoints);
    myFunc.resize(npoints);
    for (int j=0; j<npoints; j++)
    {
      fin >> myGrid(j) >> myFunc(j);
    }
    myGrid.set(myGrid(0),myGrid(npoints-1),npoints);
    T yprime_i = (myFunc(1)-myFunc(0))/(myGrid(1)-myGrid(0));
    myFunc.spline(0,yprime_i,npoints-1,0.0);
  }

  bool put(xmlNodePtr cur)
  {
  }
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
