//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_GRID_FUNCTOR_CUBIC_SPLINE_H
#define QMCPLUSPLUS_GRID_FUNCTOR_CUBIC_SPLINE_H

#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/NRSplineFunctions.h"

namespace qmcplusplus
{

/**Perform One-Dimensional Cubic Spline Interpolation using M-relation.
 *
 Given a function evaluated on a grid \f$ \{x_i\},
 i=1\ldots N, \f$ such that \f$ y_i = y(x_i), \f$ we would like to
 interpolate for a point \f$ x \f$ in the interval \f$ [x_j,x_{j+1}]. \f$
 The linear interpolation formula
 \f[
 y = Ay_j + By_{j+1}
 \f]
 where
 \f[
 A = \frac{x_{j+1}-x}{x_{j+1}-x_j} \;\;\;\;\;\;\;\;\;\;\;\;\;
 B = 1-A = \frac{x-x_{j+1}}{x_{j+1}-x_j}
 \f]
 Satisfies the conditions at the endpoints \f$ x_j \mbox{ and } x_{j+1},\f$
 but suffers from some major drawbacks.  The problem with this approach is
 that over the range of the function \f$ [x_1,x_N] \f$ we have a series of
 piecewise linear equations with a zero second derivative within each interval
 and an undefined or infinite second derivative at the interval boundaries,
 the grid points \f$ \{x_i\}. \f$  Ideally we would like to construct an
 interpolation function with a smooth first derivative and a continuous second
 derivative both within the intervals and at the the grid points.

 By adding a cubic polynomial to the linear interpolation equation within
 each interval, we can construct an interpolation function that varies
 linearly in the second derivative.  Assume for a moment that we have the
 values of the second derivative evaluated at each grid point,
 \f$ y_i'' = d^2y(x_i)/dx^2, i=1\ldots N. \f$  Now we can construct a cubic
 polynomial that has the correct second derivatives \f$y_j'' \mbox{ and }
 y_{j+1}''\f$ at the endpoints and also evaluates to zero at the endpoints.
 The reason the polynomial must be zero at the endpoints is to not spoil
 the agreement that is already built into the linear function.  A function
 constructed from these principals is given by the equation
 \f[
 y = Ay_j + By_{j+1} + Cy_j'' + Dy_{j+1}''
 \f]
 where
 \f[
 C = \frac{1}{6}(A^3-A)(x_{j+1}-x_j)^2 \;\;\;\;\;\;\;
 D = \frac{1}{6}(B^3-B)(x_{j+1}-x_j)^2.
 \f]


 To explicitly check that this function does indeed satisfy the conditions
 at the endpoints take the derivatives
 \f[
 \frac{dy}{dx} = \frac{y_{j+1}-y_j}{x_{j+1}-x_j}
 - \frac{3A^2-1}{6}(x_{j+1}-x_j)y_j''
 + \frac{3B^2-1}{6}(x_{j+1}-x_j)y_{j+1}''
 \f]
 and
 \f[
 \frac{d^2y}{dx^2} = Ay_j'' + By_{j+1}''.
 \f]
 The second derivative is continuous across the boundary between two
 intervals, e.g. \f$ [x_{j-1},x_j] \f$ and \f$ [x_j,x_{j+1}], \f$ and
 obeys the conditions at the endpoints since at \f$ x=x_j, (A=1,B=0) \f$
 and at  \f$ x=x_{j+1}, (A=0,B=1). \f$


 We had made the assumption that the values of the second derivative are
 known at the grid points, which they are not.  By imposing the condition
 that the first derivative is smooth and continuous across the boundary
 between two intervals it is possible to derive a set of equations to
 generate the \f$ y_i''\f$'s.  Evaluate the equation for the first
 derivative at \f$x=x_j\f$ in the inverval \f$ [x_{j-1},x_j] \f$ and set
 it equal to the same equation evaluated at \f$x=x_j\f$ in the inverval
 \f$ [x_j,x_{j+1}]; \f$ rearranging the terms

 \f[
 \frac{x_j-x_{j+1}}{6}y_{j+1}'' + \frac{x_{j+1}-x_{j-1}}{3}y_j''
 + \frac{x_{j+1}-x_j}{6}y_{j+1}'' = \frac{y_{j+1}-y_j}{x_{j+1}-x_j}
 -  \frac{y_j-y_{j+1}}{x_j-x_{j+1}},
 \f]
 where \f$ j=2\ldots N-1.\f$  To generate a unique solution for the system
 of \f$N-2\f$ equations we have to impose boundary conditions at \f$x_1
 \mbox{ and } x_N,\f$ the possibilities being either to set \f$y_1''
 \mbox{ and } y_N''\f$ to zero, the natural cubic spline, or, if you want
 to make the first derivative at the boundaries to have a specified value,
 use \f$y_1' \mbox{ and } y_N'\f$ to calculate the second derivatives at
 the endpoints using equation.
 *
*/

template <class Td,
         class Tg = Td,
         class CTd= Vector<Td>,
         class CTg= Vector<Tg> >
class OneDimCubicSpline: public OneDimGridFunctor<Td,Tg,CTd,CTg>
{

public:

  typedef OneDimGridFunctor<Td,Tg,CTd,CTg> base_type;
  typedef typename base_type::value_type  value_type;
  typedef typename base_type::point_type  point_type;
  typedef typename base_type::data_type data_type;
  typedef typename base_type::grid_type grid_type;

  using base_type::GridManager;
  using base_type::m_grid;
  using base_type::Y;
  using base_type::dY;
  using base_type::d2Y;
  using base_type::m_Y;
  //using base_type::FirstAddress;

  data_type m_Y2;

  using base_type::NumNodes;

  point_type r_min;
  point_type r_max;
  value_type first_deriv;
  value_type last_deriv;
  value_type ConstValue;

  //OneDimCubicSpline(const OneDimCubicSpline<Td,Tg,CTd,CTg>& rhs):
  //  base_type(rhs), m_Y2(rhs.m_Y2)
  //  { }

  OneDimCubicSpline(grid_type* gt = 0):base_type(gt) { }

  template<class VV>
  OneDimCubicSpline(grid_type* gt, const VV& nv):
    base_type(gt),first_deriv(0.0),last_deriv(0.0)
  {
    m_Y.resize(nv.size());
    m_Y2.resize(nv.size());
    copy(nv.begin(), nv.end(), m_Y.data());
  }

  OneDimCubicSpline<Td,Tg,CTd,CTg>* makeClone() const
  {
    return new OneDimCubicSpline<Td,Tg,CTd,CTg>(*this);
  }

  OneDimCubicSpline<Td,Tg,CTd,CTg>(const OneDimCubicSpline<Td,Tg,CTd,CTg>& a)
    : OneDimGridFunctor<Td,Tg,CTd,CTg>(a)
  {
    m_Y2.resize(a.m_Y2.size());
    m_Y2        = a.m_Y2;
    ConstValue  = a.ConstValue;
    r_min       = a.r_min;
    r_max       = a.r_max;
    first_deriv = a.first_deriv;
    last_deriv  = a.last_deriv;
  }

  const OneDimCubicSpline<Td,Tg,CTd,CTg>&
  operator=(const OneDimCubicSpline<Td,Tg,CTd,CTg>& a)
  {
    shallow_copy(a);
    return *this;
  }

  void shallow_copy(const OneDimCubicSpline<Td,Tg,CTd,CTg>& a)
  {
    this->GridManager = a.GridManager;
    this->OwnGrid=false;
    m_grid = a.m_grid;
    m_Y.resize(a.m_Y.size());
    m_Y2.resize(a.m_Y2.size());
    m_Y=a.m_Y;
    m_Y2=a.m_Y2;
    ConstValue = a.ConstValue;
    r_min = a.r_min;
    r_max = a.r_max;
    first_deriv = a.first_deriv;
    last_deriv = a.last_deriv;
  }

  //void setgrid(point_type r) {
  //  m_grid->locate(r);
  //}

  inline value_type splint(point_type r)
  {
    if(r<r_min)
    {
      return m_Y[0]+first_deriv*(r-r_min);
    }
    else
      if(r>=r_max)
      {
        return ConstValue;
      }
    if(GridManager)
    {
      m_grid->updateSecondOrder(r,false);
    }
    int Loc(m_grid->currentIndex());
    return m_grid->cubicInterpolateSecond(m_Y[Loc],m_Y[Loc+1],m_Y2[Loc],m_Y2[Loc+1]);
  }

  /** Interpolation to evaluate the function and itsderivatives.
   *@param r the radial distance
   *@param du return the derivative
   *@param d2u return the 2nd derivative
   *@return the value of the function
  */
  inline value_type
  splint(point_type r, value_type& du, value_type& d2u)
  {
    if(r<r_min)
      //linear-extrapolation returns y[0]+y'*(r-r[0])
    {
      du = first_deriv;
      d2u = 0.0;
      return m_Y[0]+first_deriv*(r-r_min);
    }
    else
      if(r>=r_max)
      {
        du = 0.0;
        d2u = 0.0;
        return ConstValue;
      }
    if(GridManager)
    {
      m_grid->updateSecondOrder(r,true);
    }
    int Loc(m_grid->currentIndex());
    return
      m_grid->cubicInterpolateSecond(m_Y[Loc],m_Y[Loc+1],m_Y2[Loc],m_Y2[Loc+1],du,d2u);
  }

  /** Interpolation to evaluate the function and itsderivatives.
   *@param r the radial distance
   *@param du return the derivative
   *@param d2u return the 2nd derivative
   *@param d3u return the 3nd derivative
   *@return the value of the function
  */
  inline value_type
  splint(point_type r, value_type& du, value_type& d2u, value_type& d3u)
  {
    if(r<r_min)
      //linear-extrapolation returns y[0]+y'*(r-r[0])
    {
      du = first_deriv;
      d2u = 0.0;
      return m_Y[0]+first_deriv*(r-r_min);
    }
    else
      if(r>=r_max)
      {
        du = 0.0;
        d2u = 0.0;
        return ConstValue;
      }
    if(GridManager)
    {
      m_grid->updateSecondOrder(r,true);
    }
    int Loc(m_grid->currentIndex());
    // no third derivatives yet, only for templating purposes
    d3u = 0.0;
    return
      m_grid->cubicInterpolateSecond(m_Y[Loc],m_Y[Loc+1],m_Y2[Loc],m_Y2[Loc+1],du,d2u);
  }
  /** Evaluate the 2nd derivative on the grid points
   *\param imin the index of the first valid data point
   *\param yp1 the derivative at the imin-th grid point
   *\param imax the index of the last valid data point
   *\param ypn the derivative at the imax-th grid point
   *
   *In general, a grid is shared by several OneDimCubicSpline objects
   *and each object can have its own range of valid grid points.
   *r_min and r_max are used to specify the range.
   */
  inline
  void spline(int imin, value_type yp1, int imax, value_type ypn)
  {
    first_deriv = yp1;
    last_deriv = ypn;
    r_min = m_grid->r(imin);
    r_max = m_grid->r(imax);
    int npts(this->size());
    m_Y2.resize(npts);
    m_Y2 = 0.0;
    NRCubicSpline(m_grid->data()+imin, m_Y.data()+imin,
                  npts-imin, yp1, ypn, m_Y2.data()+imin);
    ConstValue=m_Y[imax];
    //FirstAddress[0]=m_Y.data()+imin;
    //FirstAddress[2]=m_Y2.data()+imin;
  }

  inline
  void spline()
  {
    spline(0,0.0,m_grid->size()-1,0.0);
  }
};


}
#endif
