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
    
    



#ifndef CUB_SPLINE_H
#define CUB_SPLINE_H

#include "GeneralGrid.h"
#include <iostream>
#include <cstdlib>


/// The CubSpline class is a third-order spline representation of a
/// function.  It stores a pointer to a grid and the values of the
/// function and its second derivative at the points defined by the
/// grid. 
class CubSpline
{
private:
  /// This flag records whether or not the stored second derivatives
  /// are in sync with the function values.  It is used to determine
  /// whether the second derivatives need recomputation.
  bool UpToDate;
  /// The function values on the grid points.
  std::vector<double> y;   
  /// The second derivatives of the function
  std::vector<double> d2y;  

  /// The values of the derivative of the represented function on the
  /// boundary.  If each value is greater that 1e30, we compute
  /// bondary conditions assuming that the second derivative is zero at
  /// that boundary.
  double StartDeriv, EndDeriv;
public:
  SimpleGrid grid;
  bool Initialized;

  /// Returns the interpolated value.
  inline int size() { return grid.NumPoints(); }
  inline double operator()(double x);
  /// Returns the interpolated first derivative.
  inline double Deriv(double x);
  /// Returns the interpolated second derivative.
  inline double Deriv2(double x);
  /// Returns the interpolated third derivative.
  inline double Deriv3(double x);
  /// Recompute the second derivatives from the function values
  void Update();
  inline std::vector<double>& Data() { return y; } 
  
  /// Initialize the cubic spline.  See notes about start and end
  /// deriv above.
  inline void Init(SimpleGrid &newGrid, std::vector<double> yvals,
		   double startderiv, double endderiv)
  {
    StartDeriv = startderiv;
    EndDeriv   = endderiv;
    if (newGrid.NumPoints() != yvals.size()) {
      std::cerr << "Size mismatch in CubSpline.\n";
      std::cerr << "Grid Points = " << newGrid.NumPoints() << std::endl;
      std::cerr << "Y points    = " << yvals.size() << std::endl;
      abort();
    }
    grid = newGrid;
    y.resize(grid.NumPoints());
    d2y.resize(grid.NumPoints());
    y = yvals;
    Update();
    Initialized = true;
  }

  /// Simplified form which assumes that the second derivative at both
  /// boundaries are zero.
  inline void Init (SimpleGrid &newGrid, std::vector<double> &yvals)
  {
    Init (newGrid, yvals, 5.0e30, 5.0e30);
    Update();
  }
  
  /// Simplified constructor.
  inline CubSpline (SimpleGrid &newGrid, std::vector<double> &yvals)
  {
    StartDeriv = EndDeriv = 5.0e30;
    Init (newGrid, yvals, 5.0e30, 5.0e30);
    Update();
  }

  /// Full constructor.
  inline CubSpline (SimpleGrid &newGrid, std::vector<double> &yvals,
		      double startderiv, double endderiv)
  {
    Init (newGrid, yvals, startderiv, endderiv);
    Update();
  }

  /// Returns the value of the function at the ith grid point.
  inline double operator()(int i) const
  {
    return (y[i]);
  }
  /// Returns a reference to the value at the ith grid point.
  inline double & operator()(int i)
  {
    UpToDate = false;
    return (y[i]);
  }

  /// Trivial constructor
  CubSpline()
  {
    UpToDate = false;
    Initialized = false;
  }
};

inline double CubSpline::operator()(double x)
{
  if (!UpToDate)
    Update();

  SimpleGrid &X = grid;
#ifdef DEBUG
  if (x > X.End())
    {
      if (x < (X.End() * 1.000000001))
	x = X.End();
      else
	{
	  std::cerr << "x outside grid in CubSpline.\n";
	  std::cerr << "x = " << x << " X.End = " << X.End() << "\n";
	  abort();
	}
    }
#endif
  int hi = X.ReverseMap(x)+1;
  int low = hi-1;
  if (low<0) {
    low = 0;
    hi = 1;
  }
  if (hi>(X.NumPoints()-1)) {
    hi = (X.NumPoints()-1);
    low = hi-1;
  }

  double h = X[hi] - X[low];
  double hinv = 1.0/h;
  double a = (X[hi]-x)*hinv;
  double b = (x-X[low])*hinv;
  double sixinv = 0.1666666666666666666;
  
  return (a*y[low] + b*y[hi] +
	  ((a*a*a-a)*d2y[low]+(b*b*b-b)*d2y[hi])*(h*h*sixinv));
}

double CubSpline::Deriv(double x)
{
  if(!UpToDate)
    Update();

  SimpleGrid &X = grid;
  int hi = X.ReverseMap(x)+1;
  int low = hi-1;
  if (low<0) {
      low = 0;
      hi = 1;
  }
  if (hi>(X.NumPoints()-1)) {
    hi = (X.NumPoints()-1);
    low = hi-1;
  }
  
  double h = X[hi] - X[low];
  double hinv = 1.0/h;
  double a = (X[hi]-x)*hinv;
  double b = (x-X[low])*hinv;
  double sixinv = 0.1666666666666666666;
  
  return ((y[hi]-y[low])*hinv + (h*sixinv)*((3.0*b*b-1.0)*d2y[hi] -
				      (3.0*a*a-1.0)*d2y[low]));
}

inline double CubSpline::Deriv2(double x)
{
  if(!UpToDate)
    Update();
  SimpleGrid &X = grid;
  int hi = X.ReverseMap(x)+1;
  int low = hi-1;
  if (low<0) {
    low = 0;
    hi = 1;
  }
  if (hi>(X.NumPoints()-1)) {
    hi = (X.NumPoints()-1);
    low = hi-1;
  }

  
  double h = X[hi] - X[low];
  double hinv = 1.0/h;
  double a = (X[hi]-x)*hinv;
  double b = (x-X[low])*hinv;
  return (a*d2y[low] + b*d2y[hi]);
}


inline double CubSpline::Deriv3(double x)
{
  if(!UpToDate)
    Update();
  SimpleGrid &X = grid;
  int hi = X.ReverseMap(x)+1;
  int low = hi-1;
  if (low<0) {
    low = 0;
    hi = 1;
  }
  if (hi>(X.NumPoints()-1)) {
    hi = (X.NumPoints()-1);
    low = hi-1;
  }

  double h = X[hi]-X[low];
  
  return ((d2y[hi]-d2y[low])/h);
}





#endif
