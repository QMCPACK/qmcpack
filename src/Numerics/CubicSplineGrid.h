//////////////////////////////////////////////////////////////////
// (c) Copyright 2007-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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

  inline CubicSplineGrid(){}

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

  void spline(point_type start, point_type end, 
      const container_type& datain, container_type& p, container_type& dp, bool closed) 
  {
    int N(datain.size());
    if(closed) N--;
    GridStart=start;
    GridEnd=end;
    L=end-start;
    Linv=1.0/L;
    GridDelta=L/static_cast<point_type>(N);
    GridDeltaInv=1.0/GridDelta;

    p.resize(N+1);
    std::copy(datain.begin(),datain.end(),p.begin());
    dp.resize(N+1);
    container_type gr(N+1),d2p(N+1);
    //don't forget to include the end point
    for(int i=0; i<=N; i++) gr[i]=i*GridDelta+GridStart;

    p[N]=p[0];
    NRCubicSplinePBC(&(gr[0]), &(p[0]),N+1,&(dp[0]),&(d2p[0]));
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

  inline CubicSplineGrid():StartDeriv(0.0),EndDeriv(0.0){}

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

  inline bool checkGrid(point_type xIn, int& i, point_type& dl)
  {
    if(xIn<GridStart||xIn> GridEnd) return false;
    point_type delta = xIn - GridStart;
    delta -= std::floor(delta*Linv)*L;
    i = static_cast<int>(delta*GridDeltaInv);
    dl=delta*GridDeltaInv-i;
    return true;
  }

  template<typename GRIDCONTAINER>
  void setGrid(point_type start, point_type end, int n, GRIDCONTAINER& gr)
  {
    GridStart=start;
    GridEnd=end;
    GridDelta=(end-start)/static_cast<point_type>(n);
    GridDeltaInv=1.0/GridDelta;
    L=end-start;
    Linv=1.0/L;
    for(int i=0; i<=n; i++) 
      gr[i]=start+static_cast<point_type>(i)*GridDelta;
  }

  void spline(point_type start, point_type end, const container_type& datain, 
      container_type& p, container_type& d2p, bool closed) 
  {
    int n=datain.size();
    if(closed) n--;
    std::vector<point_type> gr(n+1);
    setGrid(start,end,n,gr);

    p.resize(n+1);
    d2p.resize(n+1,0.0);
    std::copy(datain.begin(), datain.end(),p.begin());
    StartDeriv=(datain[1]-datain[0])*GridDeltaInv;
    EndDeriv=0.0;//cheating
    NRCubicSpline(&gr[0],&datain[0],n+1,StartDeriv,EndDeriv,&d2p[0]);
  }

  void spline(point_type start, point_type end, const container_type& datain,  
      value_type yp1, value_type ypn,
      container_type& p, container_type& d2p, bool closed) 
  {
    int n=datain.size();
    if(closed) n--;//odd thing about the grid. 
    std::vector<point_type> gr(n+1);
    setGrid(start,end,n,gr);

    p.resize(n+1);
    d2p.resize(n+1,0.0);
    std::copy(datain.begin(), datain.end(),p.begin());
    NRCubicSpline(&gr[0],&datain[0],n+1,StartDeriv,EndDeriv,&d2p[0]);
    std::cout << "Last point " << p[n] << " " << d2p[n] << std::endl;
  }
};

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1528 $   $Date: 2006-11-23 16:57:19 -0600 (Thu, 23 Nov 2006) $
 * $Id: OneDimCubicSpline.h 1528 2006-11-23 22:57:19Z jnkim $ 
 ***************************************************************************/
