//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: D. Das, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: D. Das, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef OHMMS_QMC_DEVICE_H
#define OHMMS_QMC_DEVICE_H
#include <vector>
#include "Numerics/Spline3D/Config.h"
#include "QMCHamiltonians/WOS/Domain.h"
#include "QMCHamiltonians/WOS/WOSParticles.h"

using namespace qmcplusplus;

class Device
{

public:

  /// width of skin region
  double m_skin;

  /// dielectric of the device
  double m_eps;

  /// factor for continuous charge density
  double q_fac;

  /// continuous charge density
  double rho0;

  /// boundaries of the device : capacitor
  posvec_t r_min;
  posvec_t r_max;

  /// boundary conditions: the applied Gate voltages
  std::vector<double> Vapp;


  /// initialise and construct the Device
  Device(double delta,
         double eps,
         double qdensity,
         const std::vector<double>& appV,
         const posvec_t& min,
         const posvec_t& max)
  {
    m_skin = delta;
    m_eps = eps;
    rho0 = qdensity;
    r_min = min;
    r_max = max;
    q_fac = rho0 / ( 6.0 * m_eps );
    Vapp.resize(appV.size());
    Vapp = appV;
    //    flush();
  }

  /// make the maximum sphere, i.e. find the nearest distance to boundary
  void MaximumSphere( Domain& );

  /// check for first passage out of the boundary
  double passage( Domain& );
  double OC_passage( const double,
                     Domain&,
                     WOSParticles*);

  double contrib0( int, const Domain&, WOSParticles* );
  double contribk( const Domain&, const WOSParticles* );
  double OC_contrib0( double,
                      const posvec_t&,
                      WOSParticles*);

  void flush()
  {
    std::cout << m_skin << '\t' << m_eps << '\t' << rho0 << std::endl;
    std::cout << r_min << std::endl;
    std::cout << r_max << std::endl;
    std::cout << Vapp[4] << '\t' << Vapp[5] << std::endl;
  }


};
#endif
