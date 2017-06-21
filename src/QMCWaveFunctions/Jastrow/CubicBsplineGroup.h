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
    
    
#ifndef QMCPLUSPLUS_CUBICBSPLINEGROUP_H
#define QMCPLUSPLUS_CUBICBSPLINEGROUP_H
#include "Numerics/OneDimGridBase.h"
#include "Numerics/CubicBspline.h"
#include "Numerics/OptimizableFunctorBase.h"

namespace qmcplusplus
{

template<class T, unsigned GRIDTYPE>
class CubicBsplineGroup: public CubicBsplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>
{

public:

  typedef OptimizableFunctorBase<T>                                                FNIN;
  typedef CubicBsplineGroup<T,GRIDTYPE>                                            ThisType;
  typedef typename CubicBsplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>::point_type point_type;
  typedef typename CubicBsplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>::value_type value_type;
  typedef typename CubicBsplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>::container_type container_type;

  ///current value
  value_type Y;
  ///current first derivative
  value_type dY;
  ///current second derivative
  value_type d2Y;

  /** default constructor
   *
   * Initialize linear coefficients
   */
  inline CubicBsplineGroup(): GridManager(true),numSiblings(1),OffSet(0.0), InFunc(0)
  {
    Siblings.push_back(this);//add itself to the Siblings
  }

  /** evaluate value at x to perform SphericalBasisSet::evaluateForWalkerMove
   */
  inline value_type evaluate(point_type x, point_type xinv)
  {
    if(GridManager)
      getValuesOnly(x);
    return Y;
  }

  /** evaluate value, first and second derivates at x to perform SphericalBasisSet::evaluateForWalkerMove
   */
  inline value_type evaluateAll(point_type x, point_type xin)
  {
    if(GridManager)
      getValues(x);
    return Y;
  }

  /** Initialize the spline function
   */
  void initialize(point_type start, point_type end, const container_type& datain, bool closed,
                  T yp1, T ypn)
  {
    this->spline(start,end,yp1,ypn,datain,P);
    OffSet=datain.back();
  }

  /** Initialize the spline function with an input functor
   */
  void  initialize(FNIN* in_, point_type rmax, int npts)
  {
    if(in_==0)
      APP_ABORT("Cannot initialize with null functor");
    InFunc=in_;
    this->setGrid(0.0,rmax,npts);
    reset();
  }

  /** set GridManager
   */
  void setGridManager(bool manage)
  {
    GridManager=manage;
  }

  /** add a sibling
   */
  void addSibling(ThisType* a)
  {
    numSiblings++;
    Siblings.push_back(a);
  }

  /** reset the spline function if an inFunc is set
   */
  inline void reset()
  {
    if(InFunc)
    {
      std::cout << "Getting thru " << std::endl;
      InFunc->reset();
      container_type datain(Npts);
      point_type r=GridStart;
      for(int i=0; i<Npts; i++, r+=GridDelta)
        datain[i] = InFunc->f(r);
      this->spline(GridStart,GridEnd,InFunc->df(GridStart),InFunc->df(GridEnd),datain,P);
      OffSet=datain.back();
    }
  }


private:
  using CubicBsplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>::Npts;
  using CubicBsplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>::GridStart;
  using CubicBsplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>::GridEnd;
  using CubicBsplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>::GridDelta;
  using CubicBsplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>::GridDeltaInv;
  using CubicBsplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>::GridDeltaInv2;
  using CubicBsplineGrid<T,GRIDTYPE,FIRSTDERIV_CONSTRAINTS>::tp;

  ///boolean to indicate
  bool GridManager;
  ///the number of siblings
  int numSiblings;
  ///constant shift
  value_type OffSet;
  ///input functor to be interpolated
  FNIN* InFunc;
  /// The control points
  container_type P;
  ///a list of cubic splines which are managed by this object
  std::vector<ThisType*> Siblings;

//  /** evaluate the value of this object
//  */
//  inline void getValueOnly(int i0)
//  {
//    Y=interpolate(P[i0],P[i0+1],P[i0+2],P[i0+3]);
//  }
//
//  /** assign the default value of this object for the out-of-bound point
//  */
//  inline void getDefaultValueOnly()
//  {
//    Y=OffSet;
//  }
//
//  /** evaluate the value, first and second derivatives of this object
//  */
//  inline void getValue(int i0)
//  {
//    Y=interpolate(P[i0],P[i0+1],P[i0+2],P[i0+3],dY,d2Y);
//  }
//
//  /** assign the default value, first and second derivatives of this object
//  */
//  inline void getDefaultValue()
//  {
//    Y=OffSet; dY=0.0; d2Y=0.0;
//  }
//
  /** evaluate values of all the Siblings
  */
  inline void getValuesOnly(point_type x)
  {
    int i0;
    if(this->getGridPoint(x,i0))
    {
      for(int j=0; j<numSiblings; j++)
      {
        Siblings[j]->Y = interpolate(Sinblings[j]->P.data()+i0);
      }
    }
    else
    {
      for(int j=0; j<numSiblings; j++)
        Siblings[j]->Y=Siblings[j]->OffSet;
    }
  }

  /** evaluate values, first and second derivatives of all the Siblings
  */
  inline void getValues(point_type x)
  {
    int i0;
    if(this->getGridPoint(x,i0))
    {
      for(int j=0; j<numSiblings; j++)
        Siblings[j]->getValue(i0);
      Siblings[j]->Y = interpolate(Sinblings[j]->P.data()+i0);
    }
    else
    {
      for(int j=0; j<numSiblings; j++)
        Siblings[j]->getDefaultValue();
    }
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

  inline value_type interpolate(value_type p0, value_type p1, value_type p2, value_type p3)
  {
    const point_type onesixth=1.0/6.0;
    return onesixth*
           (tp[0]*(    -p0+3.0*p1-3.0*p2+p3)+
            tp[1]*( 3.0*p0-6.0*p1+3.0*p2)+
            tp[2]*(-3.0*p0+3.0*p2)+
            tp[3]*(     p0+4.0*p1+p2));
  }

};

}
#endif
