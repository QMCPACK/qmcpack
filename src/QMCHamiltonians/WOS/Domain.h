#ifndef GUARD_DOMAIN_H
#define GUARD_DOMAIN_H

#include <math.h>
#include <vector>
#include "Numerics/Spline3D/Config.h"
#include "Utilities/RandomGenerator.h"

struct Domain{

  /// if domain center is within HeteroStructure
  bool in_device;

  /// if domain center on an interface
  bool inter_face;

  /// id of Layer in which runner is
  int id_m;

  /// radius of spherical domain
  double radius;

  /// dielectric constant of ic layer
  double eps_d;

  /// runner at domain center
  posvec_t runner;

  /// normal distance to interfaces / domain radius
  std::vector<double> inter_frac;

  /// sampling probability of each section
  std::vector<double> sample_prob;

  /// constructore :: initilaise to be within HeteroStructure
  Domain(){
    in_device = true;
  }

  inline void WalkOnSphere(){
    
    double cosq = 2.0 * Random() - 1.0;
    double phi  = 2.0 * M_PI * Random();
    double rsinq = radius * sqrt( 1.0 - cosq * cosq );

    runner[0] += rsinq * cos( phi );
    runner[1] += rsinq * sin( phi );
    runner[2] += radius * cosq;

  }

  inline void WalkOnSphere(double theta, double phi){
    
    double rsq = radius * sin ( theta );

    runner[0] += rsq * cos( phi );
    runner[1] += rsq * sin( phi );
    runner[2] += radius * cos ( theta );

  }

  inline void WalkOnDisk(double rho, double phi, double z){

    runner[0] += rho * cos( phi );
    runner[1] += rho * sin( phi );
    runner[2] = z;

  }

  void resize( int n ){
    inter_frac.resize(n);
    sample_prob.resize(n);
  }

};
#endif
