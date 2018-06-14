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
    
    



#ifndef QMCPLUSPLUS_CUBIC_SPLINE_ENGINE_H
#define QMCPLUSPLUS_CUBIC_SPLINE_ENGINE_H

#include "Numerics/CubicSplineGrid.h"

template<class T, unsigned GRIDTYPE, unsigned BC>
struct CubicSplineEngine
{
};

template<class T, unsigned GRIDTYPE>
struct CubicSplineEngine<T,GRIDTYPE,PBC_CONSTRAINTS>
{
  typedef CubicSplineGrid<T,GRIDTYPE,PBC_CONSTRAINTS> GridType;
  typedef typename GridType::point_type point_type;
  typedef typename GridType::value_type value_type;
  typedef typename GridType::container_type container_type;

  GridType* myGrid;
  point_type cL, cR;
  point_type p1,p2,q1,q2;
  point_type dp1,dq1,dq2;
  point_type d2p1,d2q1,d2q2;

  inline CubicSplineEngine(GridType* agrid=0):myGrid(0)
  {
    if(agrid)
      setGrid(agrid);
  }

  inline void setGrid(GridType* agrid)
  {
    myGrid=agrid;
  }

  inline int getLowerGridBound(point_type x)
  {
    int i0;
    cL=myGrid->getDelta(x,i0);
    cR = cL-1.0;
    p1=cR*cR*(1.0+2.0*cL);
    p2=cL*cL*(3.0-2.0*cL);
    q1=cL*cR*cR;
    q2=cL*cL*cR;
    point_type GridDeltaInv(myGrid->GridDeltaInv);
    dp1=6.0*cL*cR*GridDeltaInv;
    dq1=1.0-4.0*cL+3.0*cL*cL;
    dq2=cL*(3.0*cL-2.0);
    d2p1=(12.0*cL-6.0)*GridDeltaInv*GridDeltaInv;
    d2q1=(6.0*cL-4.0)*GridDeltaInv;
    d2q2=(6.0*cL-2.0)*GridDeltaInv;
    return i0;
  }

  inline int getLowerGridBound0(point_type x)
  {
    int i0;
    cL=myGrid->getDelta(x,i0);
    cR = cL-1.0;
    p1=cR*cR*(1.0+2.0*cL);
    p2=cL*cL*(3.0-2.0*cL);
    q1=cL*cR*cR;
    q2=cL*cL*cR;
    return i0;
  }

  inline value_type interpolate(value_type a, value_type b, value_type a1, value_type b1,
                                value_type& du, value_type& d2u)
  {
    du = dp1*(a-b)+dq1*a1+dq2*b1;
    d2u = d2p1*(a-b)+d2q1*a1+d2q2*b1;
    return p1*a+p2*b+myGrid->GridDelta*(q1*a1+q2*b1);
  }

  inline value_type interpolate(value_type a, value_type b, value_type a1, value_type b1)
  {
    return p1*a+p2*b+myGrid->GridDelta*(q1*a1+q2*b1);
  }

  void spline(value_type yp1, value_type ypn, container_type& p, container_type& auxp)
  {
    int n=p.size();
    std::vector<point_type> gr(n);
    container_type d2p(n);
    for(int i=0; i<n; i++)
      gr[i]=myGrid->getGridValue(i);
    NRCubicSplinePBC(&(gr[0]), &(p[0]),n,&(auxp[0]),&(d2p[0]));
  }

};

template<class T, unsigned GRIDTYPE>
struct CubicSplineEngine<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>
{
  typedef CubicSplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS> GridType;
  typedef typename GridType::point_type point_type;
  typedef typename GridType::value_type value_type;
  typedef typename GridType::container_type container_type;

  GridType* myGrid;
  point_type cL, cR;
  point_type q1,q2;
  point_type dq1,dq2;

  inline CubicSplineEngine(GridType* agrid=0):myGrid(0)
  {
    if(agrid)
      setGrid(agrid);
  }

  inline void setGrid(GridType* agrid)
  {
    myGrid=agrid;
  }

  inline int getLowerGridBound(point_type x)
  {
    int i0;
    cL=myGrid->getDelta(x,i0);
    cR = 1.0-cL;
    const point_type onesixth = 0.166666666666666666666667;
    point_type GridDelta(myGrid->GridDelta);
    point_type h6(GridDelta*onesixth);
    q1 = cR*(cR*cR-1.0)*h6*GridDelta; //C
    q2 = cL*(cL*cL-1.0)*h6*GridDelta; //D
    dq1 = h6*(1.0-3.0*cR*cR);
    dq2 = h6*(3.0*cL*cL-1.0);
    return i0;
  }

  inline int getLowerGridBound0(point_type x)
  {
    int i0;
    cL=myGrid->getDelta(x,i0);
    cR = 1.0-cL;
    const point_type onesixth = 0.166666666666666666666667;
    point_type GridDelta(myGrid->GridDelta);
    point_type h6(GridDelta*GridDelta*onesixth);
    q1 = cR*(cR*cR-1.0)*h6; //C
    q2 = cL*(cL*cL-1.0)*h6; //D
    return i0;
  }

  inline value_type interpolate(value_type y1, value_type y2, value_type d2y1, value_type d2y2,
                                value_type& du, value_type& d2u)
  {
    du = myGrid->GridDeltaInv*(y2-y1)+dq1*d2y1+dq2*d2y2;
    d2u = cR*d2y1+cL*d2y2;
    return cR*y1+cL*y2+q1*d2y1+q2*d2y2;
  }

  inline value_type interpolate(value_type y1, value_type y2, value_type d2y1, value_type d2y2)
  {
    return cR*y1+cL*y2+q1*d2y1+q2*d2y2;
  }

  void spline(value_type yp1, value_type ypn, container_type& p, container_type& d2p)
  {
    int n=p.size();
    std::vector<point_type> gr(n);
    for(int i=0; i<n; i++)
      gr[i]=myGrid->getGridValue(i);
    NRCubicSpline(&gr[0],&p[0],n,yp1,ypn,&d2p[0]);
  }
};


#endif
