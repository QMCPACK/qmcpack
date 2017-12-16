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
    
    
    
    



/** @file CubicBspline.h
 * @brief declaration of CubicBspline class
 */
#ifndef QMCPLUSPLUS_CUBIC_B_SPLINE_H
#define QMCPLUSPLUS_CUBIC_B_SPLINE_H

#include "Numerics/CubicBsplineGrid.h"

template<class T, unsigned GRIDTYPE, unsigned BC>
struct CubicBspline: public CubicBsplineGrid<T,GRIDTYPE,BC>
{
  typedef typename CubicBsplineGrid<T,GRIDTYPE,BC>::point_type point_type;
  typedef typename CubicBsplineGrid<T,GRIDTYPE,BC>::value_type value_type;
  typedef typename CubicBsplineGrid<T,GRIDTYPE,BC>::container_type container_type;

  using CubicBsplineGrid<T,GRIDTYPE,BC>::tp;
  using CubicBsplineGrid<T,GRIDTYPE,BC>::GridDeltaInv;
  using CubicBsplineGrid<T,GRIDTYPE,BC>::GridDeltaInv2;

  ///index of current grid point
  int i0,i1,i2,i3;
  ///constant shift
  value_type OffSet;
  ///coefficients
  point_type A[16], dA[12], d2A[8], d3A[4];
  /// The control points
  container_type P;

  /** default constructor
   *
   * Initialize linear coefficients
   */
  inline CubicBspline():OffSet(0.0)
  {
    A[ 0] = -1.0/6.0;
    A[ 1] =  3.0/6.0;
    A[ 2] = -3.0/6.0;
    A[ 3] = 1.0/6.0;
    A[ 4] =  3.0/6.0;
    A[ 5] = -6.0/6.0;
    A[ 6] =  3.0/6.0;
    A[ 7] = 0.0/6.0;
    A[ 8] = -3.0/6.0;
    A[ 9] =  0.0/6.0;
    A[10] =  3.0/6.0;
    A[11] = 0.0/6.0;
    A[12] =  1.0/6.0;
    A[13] =  4.0/6.0;
    A[14] =  1.0/6.0;
    A[15] = 0.0/6.0;
    dA[0]=-0.5;
    dA[1]= 1.5;
    dA[ 2]=-1.5;
    dA[ 3]= 0.5;
    dA[4]= 1.0;
    dA[5]=-2.0;
    dA[ 6]= 1.0;
    dA[ 7]= 0.0;
    dA[8]=-0.5;
    dA[9]= 0.0;
    dA[10]= 0.5;
    dA[11]= 0.0;
    d2A[0]=-1.0;
    d2A[1]= 3.0;
    d2A[2]=-3.0;
    d2A[3]= 1.0;
    d2A[4]= 1.0;
    d2A[5]=-2.0;
    d2A[6]= 1.0;
    d2A[7]= 0.0;
    d3A[0]=-1.0;
    d3A[1]= 3.0;
    d3A[2]=-3.0;
    d3A[3]= 1.0;
  }

  void Init(point_type start, point_type end, const container_type& datain, bool closed)
  {
    this->spline(start,end,datain,P,closed);
  }

  void Init(point_type start, point_type end, const container_type& datain, bool closed,
            T yp1, T ypn)
  {
    this->spline(start,end,yp1,ypn,datain,P);
  }

  inline value_type getValue(point_type x)
  {
    if(this->getGridPoint(x,i0))
      return interpolate0(P[i0],P[i0+1],P[i0+2],P[i0+3]);
    else
      return OffSet;
  }

  inline value_type getDeriv(point_type x)
  {
    if(this->getGridPoint(x,i0))
      return interpolate1(P[i0],P[i0+1],P[i0+2],P[i0+3]);
    else
      return OffSet;
  }

  inline value_type getDeriv2(point_type x)
  {
    if(this->getGridPoint(x,i0))
      return interpolate2(P[i0],P[i0+1],P[i0+2],P[i0+3]);
    else
      return OffSet;
  }

  inline value_type getDeriv3(point_type x)
  {
    if(this->getGridPoint(x,i0))
      return GridDeltaInv * GridDeltaInv* GridDeltaInv*
             (tp[3]*(d2A[0]*P[i0]+d2A[1]*P[i1]+d2A[2]*P[i2]+d2A[3]*P[i3]));
    else
      return OffSet;
  }

  inline value_type operator()(T x)
  {
    return getValue(x);
  }

  inline value_type splint(T x)
  {
    if(this->getGridPoint(x,i0))
      return interpolate0(P[i0],P[i0+1],P[i0+2],P[i0+3]);
    else
      return OffSet;
  }

  inline value_type splint(point_type x, value_type& dy, value_type& d2y)
  {
    if(this->getGridPoint(x,i0))
    {
      return interpolate(P[i0],P[i0+1],P[i0+2],P[i0+3],dy,d2y);
    }
    else
    {
      dy=0.0;
      d2y=0.0;
      return OffSet;
    }
    //Too slow
    //dy= GridDeltaInv *
    //  (tp[1]*(dA[0]*P[i0]+dA[1]*P[i1]+dA[ 2]*P[i2]+dA[ 3]*P[i3])+
    //   tp[2]*(dA[4]*P[i0]+dA[5]*P[i1]+dA[ 6]*P[i2]+dA[ 7]*P[i3])+
    //   tp[3]*(dA[8]*P[i0]+dA[9]*P[i1]+dA[10]*P[i2]+dA[11]*P[i3]));
    //d2y=GridDeltaInv2 *
    //  (tp[2]*(d2A[0]*P[i0]+d2A[1]*P[i1]+d2A[2]*P[i2]+d2A[3]*P[i3])+
    //   tp[3]*(d2A[4]*P[i0]+d2A[5]*P[i1]+d2A[6]*P[i2]+d2A[7]*P[i3]));
    //return
    //  tp[0]*(A[ 0]*P[i0]+A[ 1]*P[i1]+A[ 2]*P[i2]+A[ 3]*P[i3])+
    //  tp[1]*(A[ 4]*P[i0]+A[ 5]*P[i1]+A[ 6]*P[i2]+A[ 7]*P[i3])+
    //  tp[2]*(A[ 8]*P[i0]+A[ 9]*P[i1]+A[10]*P[i2]+A[11]*P[i3])+
    //  tp[3]*(A[12]*P[i0]+A[13]*P[i1]+A[14]*P[i2]+A[15]*P[i3]);
  }
  inline value_type interpolate(value_type p0, value_type p1, value_type p2, value_type p3,
                                value_type& dy, value_type& d2y)
  {
    dy= GridDeltaInv*
        (tp[1]*(-0.5*p0+1.5*p1-1.5*p2+0.5*p3)+
         tp[2]*(     p0-2.0*p1+    p2)+
         tp[3]*(-0.5*p0       +0.5*p2));
    d2y=GridDeltaInv2*
        (tp[2]*(-p0+3.0*p1-3.0*p2+p3)+ tp[3]*(p0-2.0*p1+p2));
    const point_type onesixth=1.0/6.0;
    return onesixth*
           (tp[0]*(    -p0+3.0*p1-3.0*p2+p3)+
            tp[1]*( 3.0*p0-6.0*p1+3.0*p2)+
            tp[2]*(-3.0*p0+3.0*p2)+
            tp[3]*(     p0+4.0*p1+p2));
  }

  inline value_type interpolate0(value_type p0, value_type p1, value_type p2, value_type p3)
  {
    const point_type onesixth=1.0/6.0;
    return onesixth*
           (tp[0]*(    -p0+3.0*p1-3.0*p2+p3)+
            tp[1]*( 3.0*p0-6.0*p1+3.0*p2)+
            tp[2]*(-3.0*p0+3.0*p2)+
            tp[3]*(     p0+4.0*p1+p2));
  }

  inline value_type interpolate1(value_type p0, value_type p1, value_type p2, value_type p3)
  {
    return GridDeltaInv*
           (tp[1]*(-0.5*p0+1.5*p1-1.5*p2+0.5*p3)+
            tp[2]*(     p0-2.0*p1+    p2)+
            tp[3]*(-0.5*p0       +0.5*p2));
  }

  inline value_type interpolate2(value_type p0, value_type p1, value_type p2, value_type p3)
  {
    return GridDeltaInv2*
           (tp[2]*(-p0+3.0*p1-3.0*p2+p3)+ tp[3]*(p0-2.0*p1+p2));
  }
};


#endif
