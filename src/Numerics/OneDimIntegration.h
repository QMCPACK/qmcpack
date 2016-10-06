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
    
    



/** @file OneDimIntegration.h
 * @brief Inline functions to perform 1-D integration
 */
#ifndef QMCPLUSPLUS_ONEDIMINTEGRATION_H
#define QMCPLUSPLUS_ONEDIMINTEGRATION_H

/**@file OneDimIntegration.h
   @brief Functions to perform integration in one-dimension
   @note  The original Prim was written in F90 by Tim Wilkens.
*/

/**
 *@param f the integrand
 *@param g the integral of f
 *@return the numerical integral of f
 *@brief Performs the second order Runge-Kutta
 algorithm to evaluate the integral of a radial
 grid function
 \f$ f(r) \f$: \f[ g(r) = \int_a^r dr' f(r') \f]
*/

template<class GF>
inline
typename GF::value_type
integrate_RK2(const GF& f, GF& g)
{
  typedef typename GF::value_type value_type;
  //value_type ysum = 0.0;
  value_type yold=0.0;
  value_type ynew=0.0;
  g(0) = 0.0;
  for(int i=0; i < f.size()-1; i++)
  {
    //ysum += 0.5*f.dr(i)*(f(i)+f(i+1));
    ynew=yold+0.5*f.dr(i)*(f(i)+f(i+1));
    //g(i+1) = ysum;
    g(i+1) = ynew;
    yold = ynew;
  }
  //return ysum;
  return yold;
}

/**
 *@param f the integrand
 *@param g the integral of f
 *@return the numerical integral of f
 *@brief Performs the second order Runge-Kutta
 algorithm to evaluate the integral (in the
 forward direction) of a radial grid function
 \f$ f(r) \f$: \f[ g(r) = \int_a^r dr' f(r') \f]
*/

template<class GF>
inline
typename GF::value_type
integrate_RK2_forward(const GF& f, GF& g)
{
  return integrate_RK2(f,g);
}

/**
 *@param f the integrand
 *@param g the integral of f
 *@return the numerical integral of f
 *@brief Performs the Runge-Kutta algorithm to
 evaluate the integral (in the backwards direction)
 of a radial grid function
 \f$ f(r) \f$: \f[ g(x) = \int_x^b dx' f(x'), \f]
 where \f$ b \f$ is the upper limit of \f$ x. \f$
*/

template<class GF>
inline
typename GF::value_type
integrate_RK2_backward(const GF& f, GF& g)
{
  // typedef typename GF::value_type value_type;
  //   int last = std::min(f.size(),g.size())-1;
  //   value_type ysum=0.0;
  //   g(last) = 0.0;
  //   for(int i=last; i > 0; i--){
  //     ysum += 0.5*f.dr(i)*(f(i)+f(i-1));
  //     g(i)=ysum;
  //   }
  //   return ysum;
  typedef typename GF::value_type value_type;
  int last = std::min(f.size(),g.size())-1;
  value_type yold = 0.0;
  value_type ynew = 0.0;
  g(last) = 0.0;
  for(int i=last; i > 0; i--)
  {
    ynew = yold+0.5*f.dr(i-1)*(f(i)+f(i-1));
    g(i-1)=ynew;
    yold = ynew;
  }
  return yold;
}

/**
 *@param f the integrand
 *@return the numerical integral of f
 *@brief Performs the second order Runge-Kutta
 algorithm to evaluate the integral of a radial
 grid function
 \f$ f(r) \f$: \f[ y = \int_a^b dr' f(r'), \f]
 where \f$ (a,b) \f$ are the lower and upper
 limits of \f$ x.\f$
*/

template<class GF>
inline
typename GF::value_type
integrate_RK2(const GF& f)
{
  typedef typename GF::value_type value_type;
  value_type sum = 0.0;
  for(int i=0; i < f.size()-1; i++)
  {
    sum += f.dr(i)*(f(i)+f(i+1));
  }
  return 0.5*sum;
}

/**
 *@param f the integrand
 *@brief Normalizes the function \f$ f(r) \f$:
 \f[ f(r) = \frac{1}{\sqrt{C}} f(r) \f] where
 \f[ C = \int_a^b dr f^2(r), \f]
 where \f$ (a,b) \f$ are the lower and upper
 limits of \f$ x.\f$
*/

template<class GF>
inline
void normalize_RK2(GF& f)
{
  typedef typename GF::value_type value_type;
  value_type sum = 0.0;
  for(int i=0; i < f.size()-1; i++)
  {
    sum += f.dr(i)*(pow(f(i),2)+pow(f(i+1),2));
  }
  value_type norm = 1.0/sqrt(0.5*sum);
  for(int i=0; i < f.size(); i++)
    f(i) *= norm;
}

/**
 *@param grid the radial grid
 *@param f the integrand
 *@return the numerical integral of f
 *@brief Performs the second order Runge-Kutta
 algorithm to evaluate the integral of a radial
 grid function
 \f$ f(r) \f$: \f[ y = \int_a^b dr' f(r'), \f]
 where \f$ (a,b) \f$ are the lower and upper
 limits of \f$ x.\f$
*/

template<class GT, class Fn>
inline
typename Fn::value_type
integrate_RK2(const GT& grid, const Fn& f)
{
  typedef typename GT::value_type value_type;
  value_type sum = 0.0;
  for(int i=0; i < f.size()-1; i++)
  {
    sum += grid.dr(i)*(f(grid.r(i)) + f(grid.r(i+1)));
  }
  return 0.5*sum;
}

template<class GF>
typename GF::value_type
integrate(const GF& f)
{
  typedef typename GF::value_type value_type;
  const value_type one_third = 0.33333333333333333333333333333333;
  const value_type three_eighths = 0.375;
  const value_type BODES_FACTOR = 0.0444444444444444444444444444444;
  value_type sum = 0.0;
  int NumIntervals = f.size() - 1;
  int rem = NumIntervals % 4;
  //f.func(int i) return dr(i)*f(i)
  switch ( rem )
  {
  case 0:
    sum += 0.0;
    break;
  case 1:
    sum += 0.5 * (f.func(0)+ f.func(1));
    break;
  case 2:
    sum += one_third * (f.func(0)+  4.0 * f.func(1) + f.func(2));
    break;
  case 3:
    sum += three_eighths * (f.func(0)+3.0*f.func(1)
                            + 3.0*f.func(2) + f.func(3));
    break;
  }
  for(int i = 0; i < f.size()-1; i+=4)
  {
    //std::cout << i << '\n';
    int pt0=(i+0);
    int pt1=(i+1);
    int pt2=(i+2);
    int pt3=(i+3);
    int pt4=(i+4);
    sum +=  7.0 * f.func(pt0)  +
            32.0 * f.func(pt1) +
            12.0 * f.func(pt2) +
            32.0 * f.func(pt3) +
            7.0 * f.func(pt4);
  }
  return BODES_FACTOR * sum;
}

#endif
