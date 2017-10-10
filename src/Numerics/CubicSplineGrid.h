//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_CUBIC_SPLINE_GRID_H
#define QMCPLUSPLUS_CUBIC_SPLINE_GRID_H

#include "Numerics/GridTraits.h"
#include "Numerics/NRSplineFunctions.h"
/** CubicSplineGrid
 *
 * Empty declaration to be specialized. Three template parameters are
 * - T data type
 * - GRIDTYPE enumeration of the grid type
 * - PBC true for periodic boundary condition
 */
template<class T, unsigned GRIDTYPE, unsigned BC>
struct CubicSplineGrid { };

/** specialization for linear grid with PBC */
template<class T>
struct CubicSplineGrid<T,LINEAR_1DGRID,PBC_CONSTRAINTS>
{
  typedef typename GridTraits<T>::point_type point_type;
  typedef typename GridTraits<T>::value_type value_type;
  typedef std::vector<T>                     container_type;
  point_type GridStart, GridEnd, GridDelta, GridDeltaInv, GridDeltaInv2, L, Linv;

  inline CubicSplineGrid() {}

  inline point_type getGridValue(int i)
  {
    return GridStart+i*GridDelta;
  }

  /** evaluate i , where \f$x[i] <= xIn < x[i+1]\f$
   * @param xinx input grid point
   * @param i return index
   * @return (xIn-x[i])/(x[i+1]-x[i])
   */
  inline point_type getDelta(point_type xIn, int& i)
  {
    point_type delta = xIn - GridStart;
    delta -= std::floor(delta*Linv)*L;
    i = static_cast<int>(delta*GridDeltaInv);
    return delta*GridDeltaInv-i;
  }

  void setGrid(point_type start, point_type end, int n)
  {
    GridStart=start;
    GridEnd=end;
    L=end-start;
    Linv=1.0/L;
    GridDelta=L/static_cast<point_type>(n);
    GridDeltaInv=1.0/GridDelta;
  }
  void spline(point_type start, point_type end,
              const container_type& datain, container_type& p, container_type& dp, bool closed)
  {
    int n(datain.size());
    setGrid(start,end,(closed)?n-1:n);
    p.resize(n);
    copy(datain.begin(),datain.end(),p.begin());
    dp.resize(n);
    container_type gr(n),d2p(n);
    //don't forget to include the end point
    for(int i=0; i<n; i++)
      gr[i]=i*GridDelta+GridStart;
    p.back()=p.front();
    NRCubicSplinePBC(&(gr[0]), &(p[0]),n,&(dp[0]),&(d2p[0]));
  }
};

/** specialization for linear grid with First-Deriv */
template<class T>
struct CubicSplineGrid<T,LINEAR_1DGRID,FIRSTDERIV_CONSTRAINTS>
{
  typedef typename GridTraits<T>::point_type point_type;
  typedef typename GridTraits<T>::value_type value_type;
  typedef std::vector<T>                     container_type;
  point_type L, Linv;
  point_type GridStart, GridEnd;
  point_type GridDelta, GridDeltaInv, GridDeltaInv2;
  value_type StartDeriv;
  value_type EndDeriv;

  inline CubicSplineGrid():StartDeriv(0.0),EndDeriv(0.0) {}

  inline point_type getGridValue(int i)
  {
    return GridStart+i*GridDelta;
  }

  /** evaluate i , where \f$x[i] <= xIn < x[i+1]\f$
   * @param xinx input grid point
   * @param i return index
   * @return (xIn-x[i])/(x[i+1]-x[i])
   */
  inline point_type getDelta(point_type xIn, int& i)
  {
    point_type delta = xIn - GridStart;
    delta -= std::floor(delta*Linv)*L;
    i = static_cast<int>(delta*GridDeltaInv);
    return delta*GridDeltaInv-i;
  }

  inline int checkGrid(point_type xIn, int& i, point_type& dl)
  {
    if(xIn>GridStart && xIn<GridEnd)
    {
      point_type delta = xIn - GridStart;
      delta -= std::floor(delta*Linv)*L;
      i = static_cast<int>(delta*GridDeltaInv);
      dl=delta*GridDeltaInv-i;
      return true;
    }
    else
      return false;
  }

  /** set linear grid
   * @param start starting grid point
   * @param end ending grid point
   * @param n number of bins
   */
  void setGrid(point_type start, point_type end, int n)
  {
    GridStart=start;
    GridEnd=end;
    GridDelta=(end-start)/static_cast<point_type>(n);
    GridDeltaInv=1.0/GridDelta;
    L=end-start;
    Linv=1.0/L;
  }

  /** spline
   * @param start starting grid point
   * @param end ending grid point
   * @param p data on the grid
   * @param d2p coefficients for 2nd derivate
   * @param closed
   *
   * \if closed == true
   * p is valid in [start,end].
   * \else
   * p is valid in [start,end)
   */
  void spline(point_type start, point_type end, container_type& p, container_type& d2p, bool closed)
  {
    int n=p.size();
    setGrid(start,end,(closed)?n-1:n);
    spline((p[1]-p[0])*GridDeltaInv,0.0,p,d2p);
  }

  /** spline
   * @param start starting grid point
   * @param end ending grid point
   * @param yp1 first derivative at start
   * @param ypn first derivative at end
   * @param p data on the grid
   * @param d2p coefficients for 2nd derivate
   * @param closed
   *
   * \if closed == true
   * p is valid in [start,end].
   * \else
   * p is valid in [start,end)
   */
  void spline(point_type start, point_type end, value_type yp1, value_type ypn,
              container_type& p, container_type& d2p, bool closed)
  {
    int n=p.size();
    setGrid(start,end,(closed)?n-1:n);
    spline(yp1,ypn,p,d2p);
  }

  /** spline
   * @param yp1 first derivative at start
   * @param ypn first derivative at end
   * @param p data on the grid
   * @param d2p coefficients for 2nd derivate
   */
  void spline(value_type yp1, value_type ypn, container_type& p, container_type& d2p)
  {
    StartDeriv=yp1;
    EndDeriv=ypn;
    int n=p.size();
    std::vector<point_type> gr(n);
    for(int i=0; i<n; i++)
      gr[i]=GridStart+static_cast<point_type>(i)*GridDelta;
    NRCubicSpline(&gr[0],&p[0],n,StartDeriv,EndDeriv,&d2p[0]);
  }
};

#endif
