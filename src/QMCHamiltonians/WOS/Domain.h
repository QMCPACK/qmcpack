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
    //    cout << "cos1/phi: " << cosq << '\t' << phi << endl;
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
    return fabs(cosq);
  }




};
#endif
