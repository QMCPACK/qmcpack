//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef GUARD_CUBICSPLINE_H
#define GUARD_CUBICSPLINE_H

#include "Numerics/Spline3D/Config.h"
#include "Numerics/Spline3D/uGrid1D.h"
#include "Numerics/Spline3D/SetSplinePoint.h"
#include <fstream>
#include <blitz/array.h>

/// This is Cubic Splines with natural boundary conditions, i.e.
/// the second derivative vanishes at the boundary, not suitable for
/// functions like sines and cosines.

/// Each point of F contains:
/// 0) F(x,y,z)
/// 1) dF/dx

class CubicSpline
{

  /// functions which depend on the point where the interpolated value
  /// is required. t = (x - xi)/h
  inline double p1(double t)
  {
    return ((t-1.0)*(t-1.0)*(1.0+2.0*t));
  }
  inline double p2(double t)
  {
    return (t*t*(3.0-2.0*t));
  }
  inline double q1(double t)
  {
    return (t*(t-1.0)*(t-1.0));
  }
  inline double q2(double t)
  {
    return (t*t*(t-1.0));
  }
  inline double dp1(double t)
  {
    return (6.0*t*(t-1.0));
  }
  inline double dq1(double t)
  {
    return ((t-1.0)*(3.0*t-1.0));
  }
  inline double dp2(double t)
  {
    return (-dp1(t));
  }
  inline double dq2 (double t)
  {
    return ((3.0*t - 2.0)*t);
  }
  inline double d2p1(double t)
  {
    return (12.0*t-6.0);
  }
  inline double d2q1 (double t)
  {
    return (6.0*t - 4.0);
  }
  inline double d2p2 (double t)
  {
    return (-d2p1(t));
  }
  inline double d2q2 (double t)
  {
    return (6.0*t - 2.0);
  }

  // dim:     Dimension to calculate derivative w.r.t
  // source:  Function to differentiate
  // dest:    where to put result
  void UpdateX (int source, int dest);

  /// whether the first derivatives have been calculated using the
  /// m-relations.
  bool UpToDate;



  int n_x;
  double d2i;
  double d2f;
  double h;
  double invh;


public:

  typedef double value_type;

  /// function and derivatives at each point
  blitz::Array<blitz::TinyVector<double,2>,1> F;

  uGrid1D* m_grid;

  /// constructor
  CubicSpline(uGrid1D* agrid)
  {
    m_grid = agrid;
    n_x = agrid->m_size;
    h = agrid->m_h;
    invh = 1.0/h;
    d2i = 0.0;
    d2f = 0.0;
    F.resize(n_x);
    UpToDate = false;
  }

  inline void set_bc(double ad2i, double ad2f)
  {
    d2i = ad2i;
    d2f = ad2f;
    return;
  }

  inline double operator()(int ix) const
  {
    return (F(ix)[0]);
  }

  inline double& operator()(int ix)
  {
    UpToDate = false;
    return (F(ix)[0]);
  }

  /// update the derivatives using the m-relations
  void Update();
  void Update(double,double);

  inline double evaluate(const double x,
                         double& gradf,
                         double& lapf)
  {
    double val;
    //if( !UpToDate )      Update();              /// m-relations
    int ix = m_grid->xl(x);  /// get the lowest grid-point
    double u = (x-m_grid->m_coord[ix])*invh;
    double& Y00 = F(ix)[0];
    double& Y01 = F(ix)[1];
    double& Y10 = F(ix+1)[0];
    double& Y11 = F(ix+1)[1];
    val = Y00 * p1(u) + Y10 * p2(u) + h* ( Y01 * q1(u) + Y11 * q2(u) );
    gradf = invh * ( Y00 * dp1(u) + Y10 * dp2(u) )
            + Y01 * dq1(u) + Y11 * dq2(u) ;
    lapf = invh * ( invh * ( Y00 * d2p1(u) + Y10 * d2p2(u) )
                    + Y01 * d2q1(u) + Y11 * d2q2(u) );
    return val;
  }

};
#endif
