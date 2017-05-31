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
    
    



/* @file LinearSpline.h
 * @brief New implementation of LinearSpline with templated grid
 *
 * Cleaner design but the performance is not good compared to OneDimLinearSpline
 */
#ifndef QMCPLUSPLUS_LINEAR_SPLINE_GRIDTEMPLATED_H
#define QMCPLUSPLUS_LINEAR_SPLINE_GRIDTEMPLATED_H

#include "Numerics/GridTraits.h"

template<class T, unsigned GRIDTYPE, unsigned BC>
struct LinearSplineGrid { };

template<class T>
struct LinearSplineGrid<T,LINEAR_1DGRID,PBC_CONSTRAINTS>
{
  typedef typename GridTraits<T>::point_type point_type;
  typedef typename GridTraits<T>::value_type value_type;
  typedef std::vector<T>                     container_type;
  point_type GridStart, GridEnd, GridDelta, GridDeltaInv, L, Linv;

  inline int getGridPoint(point_type x, point_type& dh)
  {
    point_type delta = x - GridStart;
    delta -= std::floor(delta*Linv)*L;
    int i = static_cast<int>(delta*GridDeltaInv);
    dh =delta*GridDeltaInv-i;
    return i;
  }

  void spline(point_type start, point_type end,
              const container_type& data, container_type& p, bool closed)
  {
    GridStart=start;
    GridEnd=end;
    int N =data.size();
    if(closed)
      N--;
    p.resize(N+1);
    L=end-start;
    Linv=1.0/L;
    GridDelta=L/static_cast<T>(N);
    GridDeltaInv=1.0/GridDelta;
    copy(data.begin(),data.end(),p.data());
    p[N]=p[0];
  }

  void spline(const container_type& ng,
              const container_type& data, container_type& p, bool closed)
  {
    //empty
  }
};

template<class T>
struct LinearSplineGrid<T,LINEAR_1DGRID,FIRSTDERIV_CONSTRAINTS>
{
  typedef typename GridTraits<T>::point_type point_type;
  typedef typename GridTraits<T>::value_type value_type;
  typedef std::vector<T>                     container_type;
  point_type GridStart, GridEnd, GridDelta, GridDeltaInv, L, Linv;

  inline int getGridPoint(point_type x, point_type& dh)
  {
    if(x<GridStart||x>=GridEnd)
      return -1;
    point_type delta = x - GridStart;
    delta -= std::floor(delta*Linv)*L;
    int i = static_cast<int>(delta*GridDeltaInv);
    dh =delta*GridDeltaInv-i;
    return i;
  }

  void spline(point_type start, point_type end,
              const container_type& data, container_type& p, bool closed)
  {
    GridStart=start;
    GridEnd=end;
    int N =data.size();
    if(closed)
      N--;
    p.resize(N+1);
    L=end-start;
    Linv=1.0/L;
    GridDelta=L/static_cast<T>(N);
    GridDeltaInv=1.0/GridDelta;
    copy(data.begin(),data.end(),p.begin());
    //only if the last is missing
    if(!closed)
      p[N]=p[N-1];
  }

  void spline(point_type start, point_type end,
              const container_type& data, container_type& p, container_type& dp,
              bool closed)
  {
    spline(start,end,data,p,closed);
    dp.resize(p.size());
    for(int i=0; i<p.size()-1; i++)
      dp[i]=(p[i+1]-p[i])*GridDeltaInv;
  }

  void spline(const container_type& ng,
              const container_type& data, container_type& p, bool closed)
  {
    //empty
  }
};

template<class T>
struct LinearSplineGrid<T,CUSTOM_1DGRID,FIRSTDERIV_CONSTRAINTS>
{
  typedef typename GridTraits<T>::point_type point_type;
  typedef typename GridTraits<T>::value_type value_type;
  typedef std::vector<T>                     container_type;
  point_type GridStart, GridEnd, GridDelta, GridDeltaInv;

  container_type X;

  inline int getGridPoint(point_type xIn, point_type& dh)
  {
    if(xIn<GridStart||xIn>=GridEnd)
      return -1;
    int k;
    int klo=0;
    int khi=this->size()-1;
    while(khi-klo > 1)
    {
      k=(khi+klo) >> 1;
      if(X[k] > xIn)
        khi=k;
      else
        klo=k;
    }
    dh =(xIn-X[klo])/(X[klo+1]-X[klo]);
    return klo;
  }

  void spline(point_type start, point_type end,
              const container_type& data, container_type& p, bool closed)
  {
    GridStart=start;
    GridEnd=end;
    int N =data.size();
    if(closed)
      N--;
    X.resize(N+1);
    p.resize(N+1);
    point_type L=end-start;
    GridDelta=L/static_cast<T>(N);
    GridDeltaInv=1.0/GridDelta;
    for(int i=0; i<=N; i++)
      X[i]=i*GridDelta+start;
    copy(data.begin(),data.end(),p.begin());
  }

  void spline(const container_type& ng,
              const container_type& data, container_type& p, bool closed)
  {
    int N=ng.size();
    GridStart=ng[0];
    GridEnd=ng[N-1];
    X.resize(N);
    copy(ng.begin(), ng.end(), X.begin());
    p.resize(N);
    copy(data.begin(),data.end(),p.begin());
  }
};


template<class T, unsigned GRIDTYPE, unsigned BC>
struct LinearSpline: public LinearSplineGrid<T,GRIDTYPE,BC>
{
  typedef typename LinearSplineGrid<T,GRIDTYPE,BC>::point_type point_type;
  typedef typename LinearSplineGrid<T,GRIDTYPE,BC>::value_type value_type;
  typedef typename LinearSplineGrid<T,GRIDTYPE,BC>::container_type container_type;

  ///const value
  value_type ConstValue;
  ///control point data
  container_type P;

  LinearSpline():ConstValue(0.0) {}

  void Init(T start, T end, const container_type& datain, bool closed)
  {
    this->spline(start,end,datain,P,closed);
  }

  void Init(const container_type& ng, const container_type& datain, bool closed)
  {
    this->spline(ng,datain,P,closed);
  }

  inline value_type splint(point_type x)
  {
    point_type dh;
    int i=this->getGridPoint(x,dh);
    if(i<0)
      return ConstValue;
    else
      return P[i]+(P[i+1]-P[i])*dh;
  }
};
#endif
