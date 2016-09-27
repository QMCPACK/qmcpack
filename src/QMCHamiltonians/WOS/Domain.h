//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    D. Das, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef OHMMS_QMC_DOMAIN_H
#define OHMMS_QMC_DOMAIN_H
#include "Numerics/Spline3D/Config.h"
#include "Utilities/RandomGenerator.h"

struct Domain
{

  /// if runner is within the Device
  bool in_device;

  /// radius of the spherical domain
  double radius;

  /// the runner at the domain center
  posvec_t runner;

  /// the edge nearest to the runner
  int edge;

  /// constructor:: initilaise to be within Device
  Domain()
  {
    in_device = true;
  }

  inline void WalkOnSphere()
  {
    double cosq = 2.0 * Random() - 1.0;
    double phi = 2.0 * M_PI * Random();
    double rsinq = radius * sqrt( 1.0 - cosq * cosq );
    runner[0] += rsinq * cos( phi );
    runner[1] += rsinq * sin( phi );
    runner[2] += radius * cosq;
    //    std::cout << "cos1/phi: " << cosq << '\t' << phi << std::endl;
    return;
  }

  inline double ImportanceWalk()
  {
    double xi = Random();
    double cosq = sqrt(1.0-xi);
    double phi = 2.0 * M_PI * Random();
    double rsinq = radius * xi;
    if(Random() > 0.5)
      cosq = -cosq;
    runner[0] += rsinq * cos( phi );
    runner[1] += rsinq * sin( phi );
    runner[2] += radius * cosq;
    return std::abs(cosq);
  }




};
#endif
