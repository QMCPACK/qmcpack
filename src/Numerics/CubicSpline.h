//////////////////////////////////////////////////////////////////
// (c) Copyright 2007-  Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Modified by Jeongnim Kim for qmcpack
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_CUBIC_SPLINE_GRIDTEMPLATED_H
#define QMCPLUSPLUS_CUBIC_SPLINE_GRIDTEMPLATED_H

#include "Numerics/CubicSplineGrid.h"

template<class T, unsigned GRIDTYPE, bool PBC>
struct CubicSpline: public CubicSplineGrid<T,GRIDTYPE,PBC>
{
  typedef typename CubicSplineGrid<T,GRIDTYPE,PBC>::point_type point_type;
  typedef typename CubicSplineGrid<T,GRIDTYPE,PBC>::value_type value_type;
  typedef typename CubicSplineGrid<T,GRIDTYPE,PBC>::container_type container_type;

  using CubicSplineGrid<T,GRIDTYPE,PBC>::GridDelta;
  using CubicSplineGrid<T,GRIDTYPE,PBC>::GridDeltaInv;
  point_type p1,p2,q1,q2;
  point_type dp1,dq1,dq2;
  point_type d2p1,d2q1,d2q2;
  int i0;

  /// data
  container_type P;
  /// first derivatives
  container_type dP;

  inline CubicSpline() { }

  void Init(T start, T end, const container_type& datain, bool closed)
  {
    this->spline(start,end,datain,P,dP,closed);
  }

  inline void updateFirstOrder0(point_type x) {
    point_type cL = this->getDelta(x,i0); 
    point_type cR = cL-1.0;
    p1=cR*cR*(1.0+2.0*cL);
    p2=cL*cL*(3.0-2.0*cL);
    q1=cL*cR*cR;
    q2=cL*cL*cR;
  }

  inline void updateFirstOrder1(point_type x) {
    point_type cL = this->getDelta(x,i0); 
    point_type cR = cL-1.0;
    dp1=6.0*cL*cR*GridDeltaInv;
    dq1=1.0-4.0*cL+3.0*cL*cL;
    dq2=cL*(3.0*cL-2.0);
  }

  inline void updateFirstOrder2(point_type x) {
    point_type cL = this->getDelta(x,i0); 
    point_type cR = cL-1.0;
    d2p1=(12.0*cL-6.0)*GridDeltaInv*GridDeltaInv;
    d2q1=(6.0*cL-4.0)*GridDeltaInv;
    d2q2=(6.0*cL-2.0)*GridDeltaInv;
  }

  /** evaluate coefficients upto second order */
  inline void updateFirstOrder(point_type x) {
    point_type cL = this->getDelta(x,i0); 
    point_type cR = cL-1.0;
    p1=cR*cR*(1.0+2.0*cL);
    p2=cL*cL*(3.0-2.0*cL);
    q1=cL*cR*cR;
    q2=cL*cL*cR;

    dp1=6.0*cL*cR*GridDeltaInv;
    dq1=1.0-4.0*cL+3.0*cL*cL;
    dq2=cL*(3.0*cL-2.0);

    d2p1=(12.0*cL-6.0)*GridDeltaInv*GridDeltaInv;
    d2q1=(6.0*cL-4.0)*GridDeltaInv;
    d2q2=(6.0*cL-2.0)*GridDeltaInv;
  }

  inline value_type getValue(point_type x)
  {
    updateFirstOrder0(x);
    return  p1*P[i0]+p2*P[i0+1]+GridDelta*(q1*dP[i0]+q2*dP[i0+1]);
  }

  inline value_type getDeriv(point_type x)
  {
    updateFirstOrder1(x);
    return dp1*(P[i0]-P[i0+1])+dq1*dP[i0]+dq2*dP[i0+1];
  }

  inline value_type getDeriv2(point_type x)
  {
    updateFirstOrder2(x);
    return d2p1*(P[i0]-P[i0+1])+d2q1*dP[i0]+d2q2*dP[i0+1];
  }

  inline value_type operator()(T x) 
  {
    return getValue(x);
  }

  inline value_type splint(point_type x)
  {
    updateFirstOrder0(x);
    return  p1*P[i0]+p2*P[i0+1]+GridDelta*(q1*dP[i0]+q2*dP[i0+1]);
  }

  inline value_type splint(point_type x, value_type& dy, value_type& d2y)
  {
    updateFirstOrder(x);
    return interpolate(P[i0],P[i0+1],dP[i0],dP[i0+1],dy,d2y);
    //dy = dp1*(P[i0]-P[i0+1])+dq1*dP[i0]+dq2*dP[i0+1];
    //d2y = d2p1*(P[i0]-P[i0+1])+d2q1*dP[i0]+d2q2*dP[i0+1];
    //return p1*P[i0]+p2*P[i0+1]+GridDelta*(q1*dP[i0]+q2*dP[i0+1]);
  }

  inline value_type interpolate(value_type a, value_type b,
      value_type a1, value_type b1,
      value_type& dy, value_type& d2y)
  {
    dy = dp1*(a-b)+dq1*a1+dq2*b1;
    d2y = d2p1*(a-b)+d2q1*a1+d2q2*b1;
    return p1*a+p2*b+GridDelta*(q1*a1+q2*b1);
  }
};

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1580 $   $Date: 2007-01-04 10:00:43 -0600 (Thu, 04 Jan 2007) $
 * $Id: TricubicBsplineSet.h 1580 2007-01-04 16:00:43Z jnkim $
 ***************************************************************************/
