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

///dummy class for CubicSpline function
template<class T, unsigned GRIDTYPE, unsigned PBC>
struct CubicSpline { };

/** specialization for a cubic spline on any grid with PBC
 */
template<class T, unsigned GRIDTYPE>
struct CubicSpline<T,GRIDTYPE,PBC_CONSTRAINTS>: public CubicSplineGrid<T,GRIDTYPE,PBC_CONSTRAINTS>
{
  typedef typename CubicSplineGrid<T,GRIDTYPE,PBC_CONSTRAINTS>::point_type point_type;
  typedef typename CubicSplineGrid<T,GRIDTYPE,PBC_CONSTRAINTS>::value_type value_type;
  typedef typename CubicSplineGrid<T,GRIDTYPE,PBC_CONSTRAINTS>::container_type container_type;

  using CubicSplineGrid<T,GRIDTYPE,PBC_CONSTRAINTS>::GridDelta;
  using CubicSplineGrid<T,GRIDTYPE,PBC_CONSTRAINTS>::GridDeltaInv;
  point_type cL, cR;
  point_type p1,p2,q1,q2;
  point_type dp1,dq1,dq2;
  point_type d2p1,d2q1,d2q2;
  int i0;

  /// data
  container_type P;
  /// first derivatives
  container_type dP;

  inline CubicSpline() { }

  void Init(point_type start, point_type end, 
      const container_type& datain, bool closed)
  {
    this->spline(start,end,datain,P,dP,closed);
  }

  inline value_type splint(point_type x)
  {
    cL=getDelta(x,i0);
    cR = cL-1.0;
    p1=cR*cR*(1.0+2.0*cL);
    p2=cL*cL*(3.0-2.0*cL);
    q1=cL*cR*cR;
    q2=cL*cL*cR;
    return  p1*P[i0]+p2*P[i0+1]+GridDelta*(q1*dP[i0]+q2*dP[i0+1]);
  }

  inline value_type splint(point_type x, value_type& dy, value_type& d2y)
  {
    cL = this->getDelta(x,i0); 
    cR = cL-1.0;
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

  inline value_type getValue(point_type x)
  {
    point_type cL=getDelta(x,i0);
    point_type cR = cL-1.0;
    p1=cR*cR*(1.0+2.0*cL);
    p2=cL*cL*(3.0-2.0*cL);
    q1=cL*cR*cR;
    q2=cL*cL*cR;
    return  p1*P[i0]+p2*P[i0+1]+GridDelta*(q1*dP[i0]+q2*dP[i0+1]);
  }

  inline value_type getDeriv(point_type x)
  {
    point_type cL=getDelta(x,i0);
    point_type cR = cL-1.0;
    dp1=6.0*cL*cR*GridDeltaInv;
    dq1=1.0-4.0*cL+3.0*cL*cL;
    dq2=cL*(3.0*cL-2.0);
    return dp1*(P[i0]-P[i0+1])+dq1*dP[i0]+dq2*dP[i0+1];
  }

  inline value_type getDeriv2(point_type x)
  {
    point_type cL=getDelta(x,i0);
    point_type cR = cL-1.0;
    d2p1=(12.0*cL-6.0)*GridDeltaInv*GridDeltaInv;
    d2q1=(6.0*cL-4.0)*GridDeltaInv;
    d2q2=(6.0*cL-2.0)*GridDeltaInv;
    return d2p1*(P[i0]-P[i0+1])+d2q1*dP[i0]+d2q2*dP[i0+1];
  }

  inline value_type operator()(T x) 
  {
    return getValue(x);
  }
};

/** specialization for a cubic spline on any grid with PBC
 */
template<class T, unsigned GRIDTYPE>
struct CubicSpline<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>: 
public CubicSplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>
{
  typedef typename CubicSplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>::point_type point_type;
  typedef typename CubicSplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>::value_type value_type;
  typedef typename CubicSplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>::container_type container_type;

  using CubicSplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>::GridStart;
  using CubicSplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>::GridEnd;
  using CubicSplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>::GridDelta;
  using CubicSplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>::GridDeltaInv;
  using CubicSplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>::StartDeriv;
  using CubicSplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>::EndDeriv;

  int EndIndex;
  int i0;
  value_type ConstValue;
  point_type cL, cR;
  point_type q1,q2;
  point_type dq1,dq2;
  
  /// data
  container_type P;
  /// first derivatives
  container_type d2P;

  inline CubicSpline():ConstValue(0.0) { }

  void Init(point_type end, const container_type& datain, bool closed)
  {
    P=datain;
    d2P.resize(datain.size(),0.0);
    this->spline(0.0,end,P,d2P,closed);
    ConstValue=datain.back();
  }

  void Init(point_type end, const container_type& datain,  bool closed,
      value_type yp1, value_type ypn)
  {
    P=datain;
    d2P.resize(datain.size(),0.0);
    this->spline(0.0,end,yp1,ypn,P,d2P,closed);
    ConstValue=datain.back();
  }

  inline value_type splint(point_type x)
  {
    if (x>GridEnd) return ConstValue;
    cL = this->getDelta(x,i0); 
    cR = 1.0-cL;
    const point_type onesixth = 1.0/6.0;
    point_type h6(GridDelta*GridDelta*onesixth);
    q1 = cR*(cR*cR-1.0)*h6; //C
    q2 = cL*(cL*cL-1.0)*h6; //D
    return cR*P[i0]+cL*P[i0+1]+q1*d2P[i0]+q2*d2P[i0+1];
  }

  inline value_type splint(point_type x, value_type& dy, value_type& d2y)
  {
    if (x>GridEnd)
    {
      dy=0.0;d2y=0.0;
      return ConstValue;
    }
    cL = this->getDelta(x,i0); 
    cR = 1.0-cL;
    const point_type onesixth = 1.0/6.0;
    point_type h6(GridDelta*onesixth);
    q1 = cR*(cR*cR-1.0)*h6*GridDelta; //C
    q2 = cL*(cL*cL-1.0)*h6*GridDelta; //D
    dq1 = h6*(1.0-3.0*cR*cR);
    dq2 = h6*(3.0*cL*cL-1.0);

    return interpolate(P[i0],P[i0+1],d2P[i0],d2P[i0+1],dy,d2y);
  }

  inline value_type interpolate(value_type y1, value_type y2,
      value_type d2y1, value_type d2y2,
      value_type& du, value_type& d2u)
  {
    du = GridDeltaInv*(y2-y1)+dq1*d2y1+dq2*d2y2;
    d2u = cR*d2y1+cL*d2y2; 
    return cR*y1+cL*y2+q1*d2y1+q2*d2y2;
  }

};
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1580 $   $Date: 2007-01-04 10:00:43 -0600 (Thu, 04 Jan 2007) $
 * $Id: TricubicBsplineSet.h 1580 2007-01-04 16:00:43Z jnkim $
 ***************************************************************************/
