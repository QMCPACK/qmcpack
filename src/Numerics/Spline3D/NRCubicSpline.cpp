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
    
    



#include "Numerics/Spline3D/NRCubicSpline.h"


void NRCubicSpline::update()
{
  double qn, un;
  std::vector<double> u(n_x-1);
  /// set lower boundary condition
  if( y1i > y1max )
  {
    y2[0] = 0.0;
    u[0] = 0.0;
  }
  else
  {
    y2[0] = -0.5;
    u[0] = 3.0 * hin * ( y[1] - y[0] ) / ( h - y1i);
  }
  /// decomposition loop of the tridiagonal matrix
  for(int i = 1; i < n_x - 1; i++)
  {
    double sig = 0.5;
    double p = sig * y2[i-1] + 2.0;
    double pin = 1.0/p;
    y2[i] = (sig - 1.0) * pin;
    u[i] = hin * ( y[i+1] - 2.0 * y[i] + y[i-1] );
    u[i] = pin * ( 3 * hin * u[i] - sig * u[i-1] );
  }
  if( y1f > y1max )
  {
    qn = 0.0;
    un = 0.0;
  }
  else
  {
    qn = 0.5;
    un = 3.0 * hin * ( y1f - hin * ( y[n_x-1] - y[n_x-2] ) );
  }
  /// backsubstitution loop
  y2[n_x-1] = ( un - qn * u[n_x-2] ) / ( qn * y2[n_x-2] + 1.0 );
  for( int i = n_x - 2; i >= 0; i--)
    y2[i] = y2[i] * y2[i+1] + u[i];
  UpToDate = true;
  return;
}

double NRCubicSpline::evaluate( const double x,
                                double& gradf)
{
  if(! UpToDate )
    update();
  double onesixth = 1.0/6.0;
  int ix = m_grid->xl(x);
  double a = hin * ( m_grid->m_coord[ix+1] - x );
  double b = hin * ( x - m_grid->m_coord[ix] );
  double yval = a * y[ix] + b * y[ix+1] + ( a * ( a * a - 1.0 ) * y2[ix] +
                b * ( b * b - 1.0 ) * y2[ix+1] ) *
                h * h * onesixth;
  gradf = hin * ( y[ix+1] - y[ix] ) + onesixth * h *
          ( ( 3.0 * b * b - 1.0 ) * y2[ix+1] - ( 3.0 * a * a - 1.0 ) * y2[ix] );
  return yval;
}
