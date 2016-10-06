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
    
    
#ifndef GUARD_INTERFACE_H
#define GUARD_INTERFACE_H

#include "Numerics/Spline3D/Config.h"
#include "Domain.h"
#include <math.h>

struct Interface
{

  /// id_m = i is between i - 1, and i th Layers
  int id_m;

  int sign_eps;
  int xsign;

  /// z coordinate of interface
  double z_val;

  /// dielectric constant
  double eps_d;

  /// normal distance to interface / domain radius
  double d_frac;

  /// | eps(i) - eps(i-1) |
  double d_epsdiff;

  /// sampling probability
  double prob_d;

  Interface(double zvalue,
            double eps_difference,
            double eps_average)
  {
    z_val = zvalue;
    eps_d = eps_average;
    d_epsdiff = std::abs( eps_difference );
    sign_eps = sign( eps_difference );
  }


  void setid( int i )
  {
    id_m = i;
  }

  inline bool find(const Domain& domain) const
  {
    return ( domain.runner[2] == z_val );
  }

  inline bool find(const posvec_t r) const
  {
    return ( r[2] == z_val );
  }

  void calc_dfrac(double z0, double radius)
  {
    double x = ( z_val - z0 ) / radius ;
    d_frac = 0.5 * ( fabs ( x + 1.0 ) - fabs ( x - 1.0 ) );
  }

  void sample_prob(const Domain& domain)
  {
    prob_d = 0.5 * d_epsdiff * d_frac * ( 1.0 - fabs ( d_frac ) ) /
             ( domain.eps_d * ( d_frac + 1.e-10 ) );
    xsign = sign_eps * sign ( d_frac );
  }


};
#endif
