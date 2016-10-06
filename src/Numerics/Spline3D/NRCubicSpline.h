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
    
    



#ifndef GUARD_NRCUBICSPLINE_H
#define GUARD_NRCUBICSPLINE_H

#include "Numerics/Spline3D/Config.h"
#include "Numerics/Spline3D/uGrid1D.h"
#include <vector>

class NRCubicSpline
{


  /// default value of y1max
  double y1max;

  /// initial first derivative
  double y1i;

  /// final first derivative
  double y1f;

  /// function value
  std::vector<double> y;

  /// second derivative
  std::vector<double> y2;

  void update();

public:

  typedef double value_type;

  /// uniform Grid1D
  uGrid1D* m_grid;

  /// Constructor
  CubicSpline(uGrid1D* agrid)
  {
    m_grid = agrid;
    n_x = agrid->m_size;
    h = agrid->m_h;
    hin = 1.0/h;
    y.resize(n_x);
    y2.resize(n_x);
    y1max = 1.e30;
    y1i = y1max;
    y1f = y1max;
    UpToDate = false;
  }

  /// set boundary conditions on first derivatives
  inline void set_bc(double ypi, double ypf)
  {
    y1i = ypi;
    y1f = ypf;
  }


  inline evaluate(const double, double&);




};
#endif
