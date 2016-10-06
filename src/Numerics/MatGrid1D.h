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
    
    



#ifndef GUARD_MATGRID1D_H
#define GUARD_MATGRID1D_H
#include <vector>
#include "Numerics/Spline3D/Config.h"
#include "Numerics/Spline3D/Grid1D.h"


class MatGrid1D
{

public:

  /// number of different grid density sections (regions)
  int nsecs_m;

  /// starting coordinate of grid
  int ix0_m;

  /// the axis of the grid: namely the z direction z = 2;
  int iaxis;

  /// the vector containing the intervals and properties
  std::vector<int> d_ivals;
  std::vector<int> prior_m;
  std::vector<double> prop_m;
  std::vector<double> dx_ivals;


  /// the constructor
  MatGrid1D(int nsections)
  {
    nsecs_m = nsections;
    d_ivals.resize(nsecs_m+1);
    dx_ivals.resize(nsecs_m+1);
    prop_m.resize(nsecs_m);
    prior_m.resize(nsecs_m);
    std::cout << "material Grid:" << nsecs_m << std::endl;
  }

  /// initialise with the entire data from input
  inline void init(const Grid1D& aGrid1D,
                   const std::vector<int>& ix,
                   const std::vector<int>& prior,
                   const std::vector<double>& props)
  {
    d_ivals[0] = 0;
    dx_ivals[0] = aGrid1D.d_x[d_ivals[0]];
    for( int isec = 0; isec < prop_m.size(); isec++)
    {
      d_ivals[isec+1] = ix[isec+1];
      dx_ivals[isec+1] = aGrid1D.d_x[d_ivals[isec+1]];
      prop_m[isec] = props[isec];
      prior_m[isec] = prior[isec];
    }
  }

  /// get the property of the material
  inline double prop(double x)
  {
    double property;
    bool flag=true;
    for(int isec = 0; isec < nsecs_m; isec++)
    {
      double xsec = dx_ivals[isec+1];
      if( flag && x < xsec )
      {
        property = prop_m[isec];
        flag = false;
      }
      else
        if( flag && x == xsec )
        {
          int jsec = (prior_m[isec+1] > prior_m[isec]) ? isec + 1 : isec;
          property = prop_m[jsec];
          flag = false;
        }
    }
    return property;
  }

};
#endif







