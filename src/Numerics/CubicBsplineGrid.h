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
    
    



#ifndef QMCPLUSPLUS_CUBIC_B_SPLINE_GRID_H
#define QMCPLUSPLUS_CUBIC_B_SPLINE_GRID_H
#include "Numerics/GridTraits.h"
#include "Numerics/BsplineOneDimSolvers.h"
#include <limits>

/** CubicBsplineGrid
 *
 * Empty declaration to be specialized. Three template parameters are
 * - T data type
 * - GRIDTYPE enumeration of the grid type
 * - PBC true for periodic boundary condition
 */
template<class T, unsigned GRIDTYPE, unsigned BC>
struct CubicBsplineGrid { };

/** specialization for linear grid with PBC */
template<class T>
struct CubicBsplineGrid<T,LINEAR_1DGRID,PBC_CONSTRAINTS>
{
  typedef typename GridTraits<T>::point_type point_type;
  typedef typename GridTraits<T>::value_type value_type;
  typedef std::vector<T>                     container_type;
  int i0,i1,i2,i3;
  point_type GridStart, GridEnd, GridDelta, GridDeltaInv, GridDeltaInv2, L, Linv;
  point_type curPoint;
  point_type tp[4];

  inline CubicBsplineGrid():curPoint(-10000) {}

  inline bool getGridPoint(point_type x, int &i)
  {
    //point_type delta = x - GridStart;
    //delta -= std::floor(delta*Linv)*L;
    //point_type ipart;
    //point_type t = modf (delta*GridDeltaInv, &ipart);
    //int i = (int) ipart;
    point_type delta = x - GridStart;
    delta -= std::floor(delta*Linv)*L;
    i = static_cast<int>(delta*GridDeltaInv);
    point_type t =delta*GridDeltaInv-i;
    tp[0] = t*t*t;
    tp[1] = t*t;
    tp[2] = t;
    tp[3] = 1.0;
    return true;
  }

  void spline(point_type start, point_type end, const container_type& data, container_type& p, bool closed)
  {
    GridStart=start;
    GridEnd=end;
    int N =data.size();
    if(closed)
      N--;
    p.resize(N+3);
    L=end-start;
    Linv=1.0/L;
    GridDelta=L/static_cast<T>(N);
    GridDeltaInv=1.0/GridDelta;
    GridDeltaInv2=1.0/GridDelta/GridDelta;
    SolvePeriodicInterp1D<T>::apply(data,p,N);
    p[0]=p[N];
    p[N+1]=p[1];
    p[N+2]=p[2];
  }
};

/** specialization for linear grid with PBC */
template<class T>
struct CubicBsplineGrid<T,LINEAR_1DGRID,FIRSTDERIV_CONSTRAINTS>
{
  typedef typename GridTraits<T>::point_type point_type;
  typedef typename GridTraits<T>::value_type value_type;
  typedef std::vector<T>                     container_type;
  typedef std::size_t    size_t;
  ///number of points
  int Npts;
  point_type GridStart, GridEnd, GridDelta, GridDeltaInv, GridDeltaInv2, L, Linv;
  point_type tp[4];

  inline CubicBsplineGrid() {}

  inline bool getGridPoint(point_type x, int &i)
  {
    if(x<GridStart || x>GridEnd)
      return false;
    point_type delta = x - GridStart;
    delta -= std::floor(delta*Linv)*L;
    i = static_cast<int>(delta*GridDeltaInv);
    point_type t =delta*GridDeltaInv-i;
    tp[0] = t*t*t;
    tp[1] = t*t;
    tp[2] = t;
    tp[3] = 1.0;
    return true;
  }

  /** set linear grid
   * @param start starting grid
   * @param end ending grid
   * @param n size of data
   */
  void setGrid(point_type start, point_type end, size_t n)
  {
    Npts=n;
    L=end-start;
    Linv=1.0/L;
    GridDelta=L/static_cast<point_type>(Npts-1);
    GridDeltaInv=1.0/GridDelta;
    GridDeltaInv2=1.0/GridDelta/GridDelta;
    GridStart=start;
    GridEnd=end;
  }

  void spline(point_type start, point_type end, const container_type& data, container_type& p,
              bool closed)
  {
    setGrid(start,end,data.size());
    p.resize(Npts+2);
    point_type bcLower[]= {-3.0,0.0,3.0,0.0};
    point_type bcUpper[]= {-3.0,0.0,3.0,0.0};
    bcLower[3]=data[1]-data[0];
    bcUpper[3]=data[Npts-1]-data[Npts-2];
    SolveFirstDerivInterp1D<point_type>::apply(data,p,Npts,bcLower,bcUpper);
  }

  void spline(point_type start, point_type end,
              value_type startDeriv, value_type endDeriv,
              const container_type& data, container_type& p)
  {
    setGrid(start,end,data.size());
    p.resize(Npts+2);
    point_type bcLower[]= {-3.0,0.0,3.0,0.0};
    point_type bcUpper[]= {-3.0,0.0,3.0,0.0};
    bcLower[3]=startDeriv*GridDelta;
    bcUpper[3]=endDeriv*GridDelta;
    SolveFirstDerivInterp1D<point_type>::apply(data,p,Npts,bcLower,bcUpper);
  }

};
#endif
