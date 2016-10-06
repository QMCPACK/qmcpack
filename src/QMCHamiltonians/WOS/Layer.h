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
    
    
#ifndef GUARD_LAYER_H
#define GUARD_LAYER_H

#include "Numerics/Spline3D/Config.h"
#include "QMCHamiltonians/WOS/Domain.h"
#include <algorithm>


struct Layer
{

  /// layer number
  int id_m;

  /// dielectric constant
  double eps_d;

  /// Conduction Band offset
  double offset_d;

  /// how probable is sampling the layer
  double prob_d;

  /// layer boundaries
  posvec_t r_min, r_max;

  /// Constructor :: initialise data members
  Layer(double dielectric,
        double BandOffset,
        const posvec_t& bot_corner,
        const posvec_t& top_corner)
  {
    eps_d = dielectric;
    offset_d = BandOffset;
    r_min = bot_corner;
    r_max = top_corner;
  }


  void setid(int i)
  {
    id_m = i;
  }

  /// if point r is in the layer
  inline bool find(const posvec_t r) const
  {
    return (r[2] > r_min[2] && r[2] < r_max[2] );
  }

  inline bool find( const Domain& domain ) const
  {
    return find( domain.runner );
  }


  /// create spherical domain within layer
  void makeSphere(Domain& domain) const
  {
    double d[6];       /// six faces of a parallelopiped
    /// setup distances to all six faces
    d[0] = domain.runner[0] - r_min[0];
    d[1] = r_max[0]-domain.runner[0];
    d[2] = domain.runner[1]-r_min[1];
    d[3] = r_max[1]-domain.runner[1];
    d[4] = domain.runner[2]-r_min[2];
    d[5] = r_max[2]-domain.runner[2];
    /// find minimum distance to side
    double* it = std::min_element(d,d+6);
    /// domain radius is distance to nearest side
    domain.radius = (*it);
  };

  /// get sampling probability
  void sample_prob(const Domain& domain,
                   double dfrac1,
                   double dfrac2)
  {
    prob_d = 0.5 * eps_d * ( dfrac1 - dfrac2 ) / domain.eps_d;
  } ;


};
#endif
