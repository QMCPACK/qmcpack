//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include "Grid.h"

inline double mag (double x)
{
  return (std::abs(x));
}

inline double mag (const Vec2 &x)
{
  return (std::abs(x[0]) + std::abs(x[1]));
}

inline double mag (const Vec3 &x)
{
  return (std::abs(x[0])+std::abs(x[1])+std::abs(x[2]));
}

template <class IntegrandClass, class T>
class RungeKutta
{
private:
  IntegrandClass &Integrand;

  // Separate functions are just a little faster...
  inline void IntegrateForw (const Grid &grid, int startPoint, int endPoint,
			     Array<T,1> &result, bool scale=true)
  {
    T k1, k2, k3, k4;
    T y, yplus;
    double h, x;
    
    const static double OneSixth = 1.0/6.0;
    const static double OneThird = 1.0/3.0;
    
    for (int i=startPoint; i<endPoint; i++) {
      x = grid(i);
      h = grid(i+1) - x;
      y = result (i);
      k1 = h * Integrand(x,       y);
      yplus = y + 0.5*k1;
      k2 = h * Integrand(x+0.5*h, yplus);
      yplus = y + 0.5*k2;
      k3 = h * Integrand(x+0.5*h, yplus);
      yplus = y + k3;
      k4 = h * Integrand(x+h, yplus);
      result(i+1) = 
	y + (OneSixth*(k1+k4) + OneThird*(k2+k3));
      if (scale && mag(result(i+1)) > 1.0e10)
	for (int j=startPoint; j<=endPoint; j++)
	  result(j) *= 1.0e-10;
    }
  }

  inline void IntegrateRev (const Grid &grid, int startPoint, int endPoint,
			    Array<T,1> &result, bool scale=true)
  {
    T k1, k2, k3, k4;
    T y, yplus;
    double h, x;
    
    const static double OneSixth = 1.0/6.0;
    const static double OneThird = 1.0/3.0;
    
    for (int i=startPoint; i>endPoint; i--) {
      x = grid(i);
      h = grid(i-1) - x;
      y = result (i);
      k1 = h * Integrand(x,       y);
      yplus = y + 0.5*k1;
      k2 = h * Integrand(x+0.5*h, yplus);
      yplus = y + 0.5*k2;
      k3 = h * Integrand(x+0.5*h, yplus);
      yplus = y + k3;
      k4 = h * Integrand(x+h, yplus);
      result(i-1) = 
	y + (OneSixth*(k1+k4) + OneThird*(k2+k3));
      if (scale && mag(result(i-1)) > 1.0e10)
	for (int j=startPoint; j>=endPoint; j--)
	  result(j) *= 1.0e-10;
    }
  }


public:
  inline void Integrate (const Grid &grid, int startPoint, int endPoint,
			 Array<T,1> &result, bool scale=true)
  {
    if (endPoint > startPoint)
      IntegrateForw(grid, startPoint, endPoint, result, scale);
    else
      IntegrateRev (grid, startPoint, endPoint, result, scale);
  }

  RungeKutta(IntegrandClass &integrand) : Integrand(integrand)
  {
    // Do nothing 
  }
};



template <class IntegrandClass>
class RungeKutta2
{
private:
  IntegrandClass &Integrand;

  // Separate functions are just a little faster...
  inline void IntegrateForw (const Grid &grid, int startPoint, int endPoint,
			     Array<double,2> &result)
  {
    int numVars = result.cols();
    Array<double,1> k1(numVars),    k2(numVars), k3(numVars), k4(numVars);
    Array<double,1>  y(numVars), yplus(numVars);
    double h, x;
    
    const static double OneSixth = 1.0/6.0;
    const static double OneThird = 1.0/3.0;
    
    for (int i=startPoint; i<endPoint; i++) {
      x = grid(i);
      h = grid(i+1) - x;
      y = result (i,Range::all());
      k1 = h * Integrand(x,       y);
      yplus = y + 0.5*k1;
      k2 = h * Integrand(x+0.5*h, yplus);
      yplus = y + 0.5*k2;
      k3 = h * Integrand(x+0.5*h, yplus);
      yplus = y + k3;
      k4 = h * Integrand(x+h, yplus);
      result(i+1,Range::all()) = 
	y + (OneSixth*(k1+k4) + OneThird*(k2+k3));
    }
  }

  inline void IntegrateRev (const Grid &grid, int startPoint, int endPoint,
			    Array<double,2> &result)
  {
    int numVars = result.cols();
    Array<double,1> k1(numVars),    k2(numVars), k3(numVars), k4(numVars);
    Array<double,1>  y(numVars), yplus(numVars);
    double h, x;
    
    const static double OneSixth = 1.0/6.0;
    const static double OneThird = 1.0/3.0;
    
    for (int i=startPoint; i>endPoint; i--) {
      x = grid(i);
      h = grid(i-1) - x;
      y = result (i,Range::all());
      k1 = h * Integrand(x,       y);
      yplus = y + 0.5*k1;
      k2 = h * Integrand(x+0.5*h, yplus);
      yplus = y + 0.5*k2;
      k3 = h * Integrand(x+0.5*h, yplus);
      yplus = y + k3;
      k4 = h * Integrand(x+h, yplus);
      result(i-1,Range::all()) = 
	y + (OneSixth*(k1+k4) + OneThird*(k2+k3));
    }
  }


public:
  inline void Integrate (const Grid &grid, int startPoint, int endPoint,
			 Array<double,2> &result)
  {
    if (endPoint > startPoint)
      IntegrateForw(grid, startPoint, endPoint, result);
    else
      IntegrateRev (grid, startPoint, endPoint, result);
  }

  RungeKutta2(IntegrandClass &integrand) : Integrand(integrand)
  {
    // Do nothing 
  }
};


#endif
