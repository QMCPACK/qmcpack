//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "CubicSplineCommon.h"
#include "Numerics/SplineSolvers.h"

#include <iostream>

void CubicSplineCommon::Update()
{
  Grid &x = *grid;
  int N = x.NumPoints;

	CubicSplineSolve(x.data(), y.data(), N, StartDeriv, EndDeriv, d2y.data());

  UpToDate = 1;
}

//void Test2()
//{
//  LinearGrid x(0.0, 5.0, 6);
//  Array<double,1> y(6);
//  
//  y = 0.0, 2.0, -2.0, 3.0, -3.0, 1.0;
//  CubicSplineCommon spline(&x, y, 1.0, -1.0);
//  for (int i=0; i<=5000; i++)
//    {
//      double r = ((double)i/5000.01) * 5.0;
//      fprintf (stderr, "%15.12f %15.12f %15.12f %15.12f %15.12f\n",
//	       r, spline(r), spline.Deriv(r), spline.Deriv2(r),
//	       spline.Deriv3(r));  
//    }
//  
//}


//void TestSpline()
//{
//  LinearGrid x(0.0, 20.0*M_PI, 501);
//  Array<double,1> SinVals(501);

//  for (int i=0; i<501; i++)
//    SinVals(i) = sin(x(i));
//  
//  CubicSplineCommon spline (&x, SinVals);

//  //  for (int i=0; i<500; i++)
//  //  fprintf (stderr, "%15.12f %15.12f\n", spline.y(i), spline.d2y(i));

//  for (int i=0; i<=5000; i++)
//    {
//      double r = ((double)i/5000.0)*2.0*M_PI;
//      fprintf (stderr, "%15.16f %1.16f %1.16f %15.16f %15.16f %15.16f %15.16f\n",	       r, sin(r), cos(r), spline(r), spline.Deriv(r), spline.Deriv2(r),
//	       spline.Deriv3(r));  
//    }
//  
//}


//void TestMultiCubicSplineCommon()
//{
//  LinearGrid x(0.0, 2.0*M_PI, 51);
//  Array<double,2> SinCosVals(51,2);

//  for (int i=0; i<51; i++)
//    {
//      SinCosVals(i,0) = sin(x(i));
//      SinCosVals(i,1) = cos(x(i));
//    }
//  
//  Array<double,1> StartDeriv(2);
//  Array<double,1> EndDeriv(2);
//  StartDeriv(0) = 1.0;
//  StartDeriv(1) = 0.0;
//  EndDeriv(0) = 1.0;
//  EndDeriv(1) = 0.0;
//  MultiCubicSplineCommon spline (&x, SinCosVals, StartDeriv, EndDeriv);

//  //  for (int i=0; i<500; i++)
//  //  fprintf (stderr, "%15.12f %15.12f\n", spline.y(i), spline.d2y(i));

//  Array<double,1> yVal(2), yDeriv(2), yDeriv2(2), yDeriv3(2);
//  for (int i=0; i<=5000; i++)
//    {
//      double r = ((double)i/5000.0)*2.0*M_PI;
//      spline(r, yVal);
//      spline.Deriv (r, yDeriv);
//      spline.Deriv2 (r, yDeriv2);
//      spline.Deriv3 (r, yDeriv3);

//      fprintf (stderr,"%15.16f %1.16f %1.16f %15.16f %15.16f %15.16f %15.16f ",
//	       r, sin(r), cos(r), yVal(0), yDeriv(0), yDeriv2(0), yDeriv3(0));
//      fprintf (stderr,"%15.16f %15.16f %15.16f %15.16f\n",
//	       yVal(1), yDeriv(1), yDeriv2(1), yDeriv3(1));
//    }
//  
//}


CubicSplineCommon&
CubicSplineCommon::operator=(const CubicSplineCommon &spline)
{
  UpToDate = spline.UpToDate;
  y.resize(spline.y.shape());     y   = spline.y;
  d2y.resize(spline.d2y.shape()); d2y = spline.d2y;
  NumParams = spline.NumParams;
  grid = spline.grid;
  StartDeriv = spline.StartDeriv;
  EndDeriv   = spline.EndDeriv;
  
  return *this;
}


//main()
//{
  //Test2();
  //TestSpline();
//  TestMultiCubicSplineCommon();
//}

