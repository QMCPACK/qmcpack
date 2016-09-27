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

#ifndef CUBIC_SPLINE_COMMON_H
#define CUBIC_SPLINE_COMMON_H

#include "Grid.h"
#include <iostream>


/// The CubicSplineCommon class is a third-order spline representation of a
/// function.  It stores a pointer to a grid and the values of the
/// function and its second derivative at the points defined by the
/// grid. 
class CubicSplineCommon
{
private:
  /// This flag records whether or not the stored second derivatives
  /// are in sync with the function values.  It is used to determine
  /// whether the second derivatives need recomputation.
  int UpToDate;
  /// The function values on the grid points.
  Array<double, 1> y;   
  /// The second derivatives of the function
  Array<double, 1> d2y;  
public:
  int NumParams;
  Grid *grid;
  /// The values of the derivative of the represented function on the
  /// boundary.  If each value is greater that 1e30, we compute
  /// bondary conditions assuming that the second derivative is zero at
  /// that boundary.
  double StartDeriv, EndDeriv;

  inline Array<double,1>& Data();
  /// Returns the interpolated value.
  inline int size() { return grid->NumPoints; }
  inline double operator()(double x);
  /// Returns the interpolated first derivative.
  inline double Deriv(double x);
  /// Returns the interpolated second derivative.
  inline double Deriv2(double x);
  /// Returns the interpolated third derivative.
  inline double Deriv3(double x);
  /// Recompute the second derivatives from the function values
  void Update();
  
  /// Initialize the cubic spline.  See notes about start and end
  /// deriv above.
  inline void Init(Grid *NewGrid, Array<double,1> NewYs,
		   double startderiv, double endderiv)
  {
    StartDeriv = startderiv;
    EndDeriv = endderiv;
    if (NewGrid->NumPoints != NewYs.rows())
      {
	cerr << "Size mismatch in CubicSplineCommon.\n";
	cerr << "Grid Points = " << NewGrid->NumPoints << std::endl;
	cerr << "Y points    = " << NewYs.rows() << std::endl;
	exit(1);
      }
    grid = NewGrid;
    NumParams = grid->NumPoints;
    y.resize(grid->NumPoints);
    d2y.resize(grid->NumPoints);
    y = NewYs;
    Update();
  }

  /// Simplified form which assumes that the second derivative at both
  /// boundaries are zero.
  inline void Init (Grid *NewGrid, Array<double,1> NewYs)
  {
    Init (NewGrid, NewYs, 5.0e30, 5.0e30);
  }
  
  /// Simplified constructor.
  inline CubicSplineCommon (Grid *NewGrid, Array<double,1> NewYs)
  {
    StartDeriv = EndDeriv = 5.0e30;
    Init (NewGrid, NewYs, 5.0e30, 5.0e30);
  }

  /// Full constructor.
  inline CubicSplineCommon (Grid *NewGrid, Array<double,1> NewYs,
		      double startderiv, double endderiv)
  {
    Init (NewGrid, NewYs, startderiv, endderiv);
    Update();
  }

  /// Returns the value of the function at the ith grid point.
  inline double operator()(int i) const
  {
    return (y(i));
  }
  /// Returns a reference to the value at the ith grid point.
  inline double & operator()(int i)
  {
    UpToDate = 0;
    return (y(i));
  }

  /// Returns the value of the function at the ith grid point.
  inline double Params(int i) const
  {  return (y(i));  }
  /// Returns a reference to the value at the ith grid point.
  inline double & Params(int i)
  {  return ((*this)(i));   }
  void Write(IOSectionClass &outSection)
  {
    outSection.WriteVar("StartDeriv", StartDeriv);
    outSection.WriteVar("EndDeriv", EndDeriv);
    outSection.WriteVar("y", y);

    outSection.NewSection("Grid");
    grid->Write(outSection);
    outSection.CloseSection();
  }
  void Read(IOSectionClass &inSection)
  {
    assert(inSection.ReadVar("StartDeriv", StartDeriv));
    assert(inSection.ReadVar("EndDeriv", EndDeriv));
    assert(inSection.ReadVar("y", y));
    NumParams = y.size();
    d2y.resize(NumParams);
    assert(inSection.OpenSection("Grid"));
    grid = ReadGrid(inSection);
    inSection.CloseSection();
    Update();
  }

  CubicSplineCommon& operator= (const CubicSplineCommon& spline);

  /// Trivial constructor
  CubicSplineCommon()
  {
    UpToDate = 0;
    grid = NULL;
  }
};


inline Array<double,1>& CubicSplineCommon::Data()
{
  UpToDate = false;
  return y;
}



inline double CubicSplineCommon::operator()(double x)
{
  if (!UpToDate)
    Update();


  Grid &X = *grid;
#ifdef BZ_DEBUG
  if (x > X.End)
    {
      if (x < (X.End * 1.000000001))
	x = X.End;
      else
	{
	  std::cerr << "x outside grid in CubicSplineCommon.\n";
	  std::cerr << "x = " << x << " X.End = " << X.End << "\n";
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
  if (hi>(X.NumPoints-1)) {
    hi = (X.NumPoints-1);
    low = hi-1;
  }

  double h = X(hi) - X(low);
  double hinv = 1.0/h;
  double a = (X(hi)-x)*hinv;
  double b = (x-X(low))*hinv;
  double sixinv = 0.1666666666666666666;
  
  return (a*y(low) + b*y(hi) +
	  ((a*a*a-a)*d2y(low)+(b*b*b-b)*d2y(hi))*(h*h*sixinv));
}

double CubicSplineCommon::Deriv(double x)
{
  if(!UpToDate)
    Update();

  Grid &X = *grid;
  int hi = X.ReverseMap(x)+1;
  int low = hi-1;
  if (low<0) {
      low = 0;
      hi = 1;
  }
  if (hi>(X.NumPoints-1)) {
    hi = (X.NumPoints-1);
    low = hi-1;
  }
  
  double h = X(hi) - X(low);
  double hinv = 1.0/h;
  double a = (X(hi)-x)*hinv;
  double b = (x-X(low))*hinv;
  double sixinv = 0.1666666666666666666;
  
  return ((y(hi)-y(low))*hinv + (h*sixinv)*((3.0*b*b-1.0)*d2y(hi) -
				      (3.0*a*a-1.0)*d2y(low)));
}

inline double CubicSplineCommon::Deriv2(double x)
{
  if(!UpToDate)
    Update();
  Grid &X = *grid;
  int hi = X.ReverseMap(x)+1;
  int low = hi-1;
  if (low<0) {
    low = 0;
    hi = 1;
  }
  if (hi>(X.NumPoints-1)) {
    hi = (X.NumPoints-1);
    low = hi-1;
  }

  
  double h = X(hi) - X(low);
  double hinv = 1.0/h;
  double a = (X(hi)-x)*hinv;
  double b = (x-X(low))*hinv;
  return (a*d2y(low) + b*d2y(hi));
}


inline double CubicSplineCommon::Deriv3(double x)
{
  if(!UpToDate)
    Update();
  Grid &X = *grid;
  int hi = X.ReverseMap(x)+1;
  int low = hi-1;
  if (low<0) {
    low = 0;
    hi = 1;
  }
  if (hi>(X.NumPoints-1)) {
    hi = (X.NumPoints-1);
    low = hi-1;
  }

  double h = X(hi)-X(low);
  
  return ((d2y(hi)-d2y(low))/h);
}





/// The MulitCubicSpline class is nearly identical to the CubicSplineCommon
/// class, except that it stores the values of several functions at
/// the same grid point.  This can be used to simultaneously
/// interpolate several functions at the same point without redundant
/// computation.  
class MultiCubicSplineCommon
{
private:
  /// Stores whether each functions second derivatives are in sync
  /// with the stored function values.
  Array<bool,1> UpToDate;
  Array<double, 2> y;   //< The values of the function at the data points
  Array<double, 2> d2y;  //< The second derivatives of the function
public:
  int NumGridPoints, NumSplines;
  Grid *grid;
  Array<double, 1> StartDeriv, EndDeriv;

  inline double operator()(int, double x);
  inline void operator()(double x, Array<double,1> &yVec);
  inline double Deriv(int i, double x);
  inline void   Deriv (double x, Array<double,1> &deriv);
  inline double Deriv2(int i, double x);
  inline void   Deriv2 (double x, Array<double,1> &deriv2);
  inline double Deriv3(int i, double x);
  inline void   Deriv3 (double x, Array<double,1> &deriv3);
  void Update(int i);
  
  inline void Init(Grid *NewGrid, const Array<double,2> &NewYs,
		   const Array<double,1> &startderiv, 
		   const Array<double,1> &endderiv)
  {

    grid = NewGrid;
    NumGridPoints = grid->NumPoints;
    NumSplines = NewYs.cols(); 
    UpToDate.resize(NumSplines);
    StartDeriv.resize(NumSplines);

    StartDeriv = startderiv;

    EndDeriv.resize(NumSplines);

    EndDeriv = endderiv;

    if (NewGrid->NumPoints != NewYs.rows())
      {
	cerr << "Size mismatch in CubicSplineCommon.\n";
	cerr << "Grid Points = " << NewGrid->NumPoints << std::endl;
	cerr << "Y points    = " << NewYs.rows() << std::endl;
	exit(1);
      }

    y.resize(NumGridPoints,NumSplines);
    d2y.resize(NumGridPoints,NumSplines);

    y = NewYs;

    for (int i=0; i<NumSplines; i++)
      Update(i);

  }

  inline void Init (Grid *NewGrid, const Array<double,2> &NewYs)
  {

    grid = NewGrid;
    NumGridPoints = grid->NumPoints;
    NumSplines = NewYs.cols();
    UpToDate.resize(NumSplines);
    StartDeriv.resize(NumSplines);
    EndDeriv.resize(NumSplines);
    StartDeriv = 5.0e30;
    EndDeriv = 5.0e30;

    if (NewGrid->NumPoints != NewYs.rows())
      {
	cerr << "Size mismatch in CubicSplineCommon.\n";
	exit(1);
      }

    y.resize(NumGridPoints,NumSplines);
    d2y.resize(NumGridPoints,NumSplines);

    y = NewYs;

    for (int i=0; i<NumSplines; i++)
      Update(i);

  }

  inline MultiCubicSplineCommon (Grid *NewGrid, Array<double,2> NewYs,
			   Array<double,1> startderiv, 
			   Array<double,1> endderiv)
  {
    Init (NewGrid, NewYs, startderiv, endderiv);
  }

  inline MultiCubicSplineCommon (Grid *NewGrid, Array<double,2> NewYs)
  {
    Init (NewGrid, NewYs);
  }

  inline double Params(int i, int j) const
  { return (y(i,j));  }
  
  inline double & Params(int i, int j)
  {
    UpToDate(j) = false;
    return (y(i,j));
  }

  inline double operator()(int i, int j) const
  { return (y(i,j)); }

  inline double & operator()(int i, int j)
  { 
    UpToDate(j) = false; 
    return (y(i,j)); 
  }
  /// Copy constructor
  MultiCubicSplineCommon (const MultiCubicSplineCommon &a)
  {
    UpToDate.resize(a.UpToDate.rows());
    UpToDate = a.UpToDate;
    y.resize(a.y.rows(), a.y.cols());
    y = a.y;
    d2y.resize(a.d2y.rows(), a.d2y.cols());
    d2y = a.d2y;
    NumGridPoints = a.NumGridPoints;
    NumSplines = a.NumSplines;
    grid = a.grid;
    StartDeriv.resize(a.StartDeriv.rows());
    StartDeriv = a.StartDeriv;
    EndDeriv.resize(a.EndDeriv.rows());
    EndDeriv = a.EndDeriv;
  }

  /// Assignment operator -- necessary for array resizeAndPreserve
  MultiCubicSplineCommon & operator= (MultiCubicSplineCommon &a)
  {
    UpToDate.resize(a.UpToDate.rows());
    UpToDate = a.UpToDate;
    y.resize(a.y.rows(), a.y.cols());
    y = a.y;
    d2y.resize(a.d2y.rows(), a.d2y.cols());
    d2y = a.d2y;
    NumGridPoints = a.NumGridPoints;
    NumSplines = a.NumSplines;
    grid = a.grid;
    StartDeriv.resize(a.StartDeriv.rows());
    StartDeriv = a.StartDeriv;
    EndDeriv.resize(a.EndDeriv.rows());
    EndDeriv = a.EndDeriv;
    return *this;
  }

  /// Assignment operator -- necessary for array resizeAndPreserve
  MultiCubicSplineCommon & operator= (MultiCubicSplineCommon a)
  {
    UpToDate.resize(a.UpToDate.rows());
    UpToDate = a.UpToDate;
    y.resize(a.y.rows(), a.y.cols());
    y = a.y;
    d2y.resize(a.d2y.rows(), a.d2y.cols());
    d2y = a.d2y;
    NumGridPoints = a.NumGridPoints;
    NumSplines = a.NumSplines;
    grid = a.grid;
    StartDeriv.resize(a.StartDeriv.rows());
    StartDeriv = a.StartDeriv;
    EndDeriv.resize(a.EndDeriv.rows());
    EndDeriv = a.EndDeriv;
    return *this;
  }


  MultiCubicSplineCommon()
  {
  }
};



//////////////////////////////////////////////////////////////////////
// MultiCubicSplineCommon operator
//////////////////////////////////////////////////////////////////////
inline double MultiCubicSplineCommon::operator()(int i, double x)
{
  if (!UpToDate(i))
    Update(i);

  Grid &X = *grid;
#ifdef BZ_DEBUG
  if (x > X.End)
    {
      if (x < (X.End * 1.000000001))
	x = X.End;
      else
	{
	  std::cerr << "x outside grid in CubicSplineCommon.\n";
	  std::cerr << "x = " << x << " X.End = " << X.End << "\n";
	  exit(1);
	}
    }
#endif

  int hi = X.ReverseMap(x)+1;
  int low = hi-1;
  if (low<0)
    {
      low = 0;
      hi = 1;
    }
  if (hi>(X.NumPoints-1))
    {
      hi = (X.NumPoints-1);
      low = hi-1;
    }

  double h = X(hi) - X(low);
  double hinv = 1.0/h;
  double a = (X(hi)-x)*hinv;
  double b = (x-X(low))*hinv;
  double sixinv = 0.1666666666666666666;
  
  return (a*y(low,i) + b*y(hi,i) +
	  ((a*a*a-a)*d2y(low,i)+(b*b*b-b)*d2y(hi,i))*(h*h*sixinv));
}


inline void MultiCubicSplineCommon::operator()(double x, Array<double,1> &yVec)
{
  for (int i=0; i<NumSplines; i++)
    if (!UpToDate(i))
      Update(i);

  Grid &X = *grid;
#ifdef BZ_DEBUG
  if (x > X.End)
    {
      if (x < (X.End * 1.000000001))
	x = X.End;
      else
	{
	  std::cerr << "x outside grid in CubicSplineCommon.\n";
	  std::cerr << "x = " << x << " X.End = " << X.End << "\n";
	  exit(1);
	}
    }
#endif
  int hi = X.ReverseMap(x)+1;
  int low = hi-1;
  if (low<0)
    {
      low = 0;
      hi = 1;
    }
  if (hi>(X.NumPoints-1))
    {
      hi = (X.NumPoints-1);
      low = hi-1;
    }

  double h = X(hi) - X(low);
  double hinv = 1.0/h;
  double a = (X(hi)-x)*hinv;
  double b = (x-X(low))*hinv;
  
  double a3minusa = a*a*a-a;
  double b3minusb = b*b*b-b;
  double h2over6 = h*h*0.1666666666666666666;
  
  for (int i=0; i<NumSplines; i++)
    {
      yVec(i) = a*y(low,i) + b*y(hi,i) +
	(a3minusa*d2y(low,i)+b3minusb*d2y(hi,i))*h2over6;
    }
}





inline double MultiCubicSplineCommon::Deriv(int i, double x)
{
  if(!UpToDate(i))
    Update(i);

  Grid &X = *grid;
  int hi = X.ReverseMap(x)+1;
  int low = hi-1;
  if (low<0)
    {
      low = 0;
      hi = 1;
    }
  
  double h = X(hi) - X(low);
  double hinv = 1.0/h;
  double a = (X(hi)-x)*hinv;
  double b = (x-X(low))*hinv;
  double sixinv = 0.166666666666666666666;
  
  return ((y(hi,i)-y(low,i))*hinv + (h*sixinv)*((3.0*b*b-1.0)*d2y(hi,i) -
						(3.0*a*a-1.0)*d2y(low,i)));
}


inline void MultiCubicSplineCommon::Deriv(double x, Array<double,1> &deriv)
{
  for (int i=0; i<NumSplines; i++)
    if(!UpToDate(i))
      Update(i);

  Grid &X = *grid;
  int hi = X.ReverseMap(x)+1;
  int low = hi-1;
  if (low<0)
    {
      low = 0;
      hi = 1;
    }
  
  double h = X(hi) - X(low);
  double hinv = 1.0/h;
  double a = (X(hi)-x)*hinv;
  double b = (x-X(low))*hinv;
  
  double a2times3 = 3.0*a*a;
  double b2times3 = 3.0*b*b;
  double hover6 = h*0.1666666666666666666666;

  for (int i=0; i<NumSplines; i++)
    {
      deriv(i) = (y(hi,i)-y(low,i))*hinv + 
	hover6 * ((b2times3-1.0)*d2y(hi,i) -
		  (a2times3-1.0)*d2y(low,i));
    }
}



inline double MultiCubicSplineCommon::Deriv2(int i, double x)
{
  if(!UpToDate(i))
    Update(i);
  Grid &X = *grid;
  int hi = X.ReverseMap(x)+1;
  int low = hi-1;
  if (low<0)
    {
      low = 0;
      hi = 1;
    }
  
  double h = X(hi) - X(low);
  double hinv = 1.0/h;
  double a = (X(hi)-x)*hinv;
  double b = (x-X(low))*hinv;
  return (a*d2y(low, i) + b*d2y(hi, i));
}


inline void MultiCubicSplineCommon::Deriv2(double x, Array<double,1> &deriv2)
{
  for (int i=0; i<NumSplines; i++)  
    if(!UpToDate(i))
      Update(i);
  Grid &X = *grid;
  int hi = X.ReverseMap(x)+1;
  int low = hi-1;
  if (low<0)
    {
      low = 0;
      hi = 1;
    }
  
  double h = X(hi) - X(low);
  double hinv = 1.0/h;
  double a = (X(hi)-x)*hinv;
  double b = (x-X(low))*hinv;

  for (int i=0; i<NumSplines; i++)  
    deriv2(i) = a*d2y(low,i) + b*d2y(hi,i);
}



inline double MultiCubicSplineCommon::Deriv3(int i, double x)
{
  if(!UpToDate(i))
    Update(i);
  Grid &X = *grid;
  int hi = X.ReverseMap(x)+1;
  int low = hi-1;
  if (low<0)
    {
      low = 0;
      hi = 1;
    }
  double h = X(hi)-X(low);
  
  return ((d2y(hi, i)-d2y(low, i))/h);
}


inline void MultiCubicSplineCommon::Deriv3(double x, Array<double,1> &deriv3)
{
  for (int i=0; i<NumSplines; i++)
    if(!UpToDate(i))
      Update(i);
  Grid &X = *grid;
  int hi = X.ReverseMap(x)+1;
  int low = hi-1;
  if (low<0)
    {
      low = 0;
      hi = 1;
    }
  double h = X(hi)-X(low);
  double hinv = 1.0/h;

  for (int i=0; i<NumSplines; i++)
    deriv3(i) = (d2y(hi,i) - d2y(low,i))*hinv;
}





#endif
