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
    
    
    
    



#ifndef QMCPLUSPLUS_CUBIC_SPLINE_GRIDTEMPLATED_H
#define QMCPLUSPLUS_CUBIC_SPLINE_GRIDTEMPLATED_H

#include "Numerics/CubicSplineGrid.h"
#include "Numerics/CubicSplineEngine.h"

//template<class T>
template<class T, unsigned GRIDTYPE, unsigned BC>
struct CubicSpline
{
  typedef CubicSplineGrid<T,GRIDTYPE,BC> GridType;
  typedef CubicSplineEngine<T,GRIDTYPE,BC> SplineEngineType;
  typedef typename GridType::point_type point_type;
  typedef typename GridType::value_type value_type;
  typedef typename GridType::container_type container_type;

  int EndIndex;
  ///vconst value
  value_type ConstValue;
  ///grid
  GridType myGrid;
  ///spline engine
  //CubicSplineFirst<GridType> myEngine;
  SplineEngineType myEngine;
  /// data
  container_type P;
  /// supporting coeffs
  container_type auxP;

  inline CubicSpline():ConstValue(0.0) { }

  void Init(point_type start, point_type end, const container_type& datain,  bool closed,
            value_type yp1, value_type ypn)
  {
    myEngine.setGrid(&myGrid);
    int n=datain.size();
    myGrid.setGrid(start,end,(closed)?n-1:n);
    ConstValue=datain.back();
    EndIndex = (closed)?n-1:n-2;
    P=datain;
    auxP.resize(datain.size(),0.0);
    myEngine.spline(yp1,ypn,P,auxP);
  }

  void Init(point_type start, point_type end, const container_type& datain,  bool closed)
  {
    Init(start,end,datain,closed,0.0,0.0);
  }

  inline value_type splint(point_type x)
  {
    int i0=myEngine.getLowerGridBound0(x);
    if(i0>EndIndex)
      return ConstValue;
    else
      return myEngine.interpolate(P[i0],P[i0+1],auxP[i0],auxP[i0+1]);
  }

  inline value_type splint(point_type x, value_type& dy, value_type& d2y)
  {
    int i0=myEngine.getLowerGridBound(x);
    if(i0>EndIndex)
    {
      dy=0.0;
      d2y=0.0;
      return ConstValue;
    }
    else
    {
      return myEngine.interpolate(P[i0],P[i0+1],auxP[i0],auxP[i0+1],dy,d2y);
    }
  }

};
#endif
