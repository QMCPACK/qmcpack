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
  std::shared_ptr<Grid> grid;
  /// The values of the derivative of the represented function on the
  /// boundary.  If each value is greater that 1e30, we compute
  /// boundary conditions assuming that the second derivative is zero at
  /// that boundary.
  double StartDeriv, EndDeriv;

  inline Array<double, 1>& Data();
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
  inline void Init(std::shared_ptr<Grid>& NewGrid, Array<double, 1> NewYs, double startderiv, double endderiv)
  {
    StartDeriv = startderiv;
    EndDeriv   = endderiv;
    if (NewGrid->NumPoints != NewYs.rows())
    {
      std::cerr << "Size mismatch in CubicSplineCommon.\n";
      std::cerr << "Grid Points = " << NewGrid->NumPoints << std::endl;
      std::cerr << "Y points    = " << NewYs.rows() << std::endl;
      exit(1);
    }
    grid      = NewGrid;
    NumParams = grid->NumPoints;
    y.resize(grid->NumPoints);
    d2y.resize(grid->NumPoints);
    y = NewYs;
    Update();
  }

  /// Simplified form which assumes that the second derivative at both
  /// boundaries are zero.
  inline void Init(std::shared_ptr<Grid>& NewGrid, Array<double, 1> NewYs) { Init(NewGrid, NewYs, 5.0e30, 5.0e30); }

  /// Simplified constructor.
  inline CubicSplineCommon(std::shared_ptr<Grid>& NewGrid, Array<double, 1> NewYs)
  {
    StartDeriv = EndDeriv = 5.0e30;
    Init(NewGrid, NewYs, 5.0e30, 5.0e30);
  }

  /// Full constructor.
  inline CubicSplineCommon(std::shared_ptr<Grid>& NewGrid, Array<double, 1> NewYs, double startderiv, double endderiv)
  {
    Init(NewGrid, NewYs, startderiv, endderiv);
    Update();
  }

  /// Returns the value of the function at the ith grid point.
  inline double operator()(int i) const { return (y(i)); }
  /// Returns a reference to the value at the ith grid point.
  inline double& operator()(int i)
  {
    UpToDate = 0;
    return (y(i));
  }

  /// Returns the value of the function at the ith grid point.
  inline double Params(int i) const { return (y(i)); }
  /// Returns a reference to the value at the ith grid point.
  inline double& Params(int i) { return (y(i)); }
  void Write(IOSectionClass& outSection)
  {
    outSection.WriteVar("StartDeriv", StartDeriv);
    outSection.WriteVar("EndDeriv", EndDeriv);
    outSection.WriteVar("y", y);

    outSection.NewSection("Grid");
    grid->Write(outSection);
    outSection.CloseSection();
  }
  void Read(IOSectionClass& inSection)
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

  CubicSplineCommon& operator=(const CubicSplineCommon& spline);

  /// Trivial constructor
  CubicSplineCommon()
  {
    UpToDate = 0;
    NumParams = 0;
    StartDeriv = 0;
    EndDeriv = 0;
    grid     = NULL;
  }
};


inline Array<double, 1>& CubicSplineCommon::Data()
{
  UpToDate = false;
  return y;
}


inline double CubicSplineCommon::operator()(double x)
{
  if (!UpToDate)
    Update();


  Grid& X = *grid;
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
  int hi  = X.ReverseMap(x) + 1;
  int low = hi - 1;
  if (low < 0)
  {
    low = 0;
    hi  = 1;
  }
  if (hi > (X.NumPoints - 1))
  {
    hi  = (X.NumPoints - 1);
    low = hi - 1;
  }

  double h      = X(hi) - X(low);
  double hinv   = 1.0 / h;
  double a      = (X(hi) - x) * hinv;
  double b      = (x - X(low)) * hinv;
  double sixinv = 0.1666666666666666666;

  return (a * y(low) + b * y(hi) + ((a * a * a - a) * d2y(low) + (b * b * b - b) * d2y(hi)) * (h * h * sixinv));
}

double CubicSplineCommon::Deriv(double x)
{
  if (!UpToDate)
    Update();

  Grid& X = *grid;
  int hi  = X.ReverseMap(x) + 1;
  int low = hi - 1;
  if (low < 0)
  {
    low = 0;
    hi  = 1;
  }
  if (hi > (X.NumPoints - 1))
  {
    hi  = (X.NumPoints - 1);
    low = hi - 1;
  }

  double h      = X(hi) - X(low);
  double hinv   = 1.0 / h;
  double a      = (X(hi) - x) * hinv;
  double b      = (x - X(low)) * hinv;
  double sixinv = 0.1666666666666666666;

  return ((y(hi) - y(low)) * hinv + (h * sixinv) * ((3.0 * b * b - 1.0) * d2y(hi) - (3.0 * a * a - 1.0) * d2y(low)));
}

inline double CubicSplineCommon::Deriv2(double x)
{
  if (!UpToDate)
    Update();
  Grid& X = *grid;
  int hi  = X.ReverseMap(x) + 1;
  int low = hi - 1;
  if (low < 0)
  {
    low = 0;
    hi  = 1;
  }
  if (hi > (X.NumPoints - 1))
  {
    hi  = (X.NumPoints - 1);
    low = hi - 1;
  }


  double h    = X(hi) - X(low);
  double hinv = 1.0 / h;
  double a    = (X(hi) - x) * hinv;
  double b    = (x - X(low)) * hinv;
  return (a * d2y(low) + b * d2y(hi));
}


inline double CubicSplineCommon::Deriv3(double x)
{
  if (!UpToDate)
    Update();
  Grid& X = *grid;
  int hi  = X.ReverseMap(x) + 1;
  int low = hi - 1;
  if (low < 0)
  {
    low = 0;
    hi  = 1;
  }
  if (hi > (X.NumPoints - 1))
  {
    hi  = (X.NumPoints - 1);
    low = hi - 1;
  }

  double h = X(hi) - X(low);

  return ((d2y(hi) - d2y(low)) / h);
}

#endif
