//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: D. Das, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: D. Das, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCHamiltonians/WOS/Device.h"
#include <algorithm>

// void Device::MaximumSphere( Domain& domain){

//   double d[6];
//   /// setup distances to all six faces
//   d[0] = domain.runner[0] - r_min[0];
//   d[1] = r_max[0] - domain.runner[0];
//   d[2] = domain.runner[1] - r_min[1];
//   d[3] = r_max[1] - domain.runner[1];
//   d[4] = domain.runner[2] - r_min[2];
//   d[5] = r_max[2] - domain.runner[2];

//   /// find minimum distance to side
//   double* it = std::min_element(d,d+6);

//   /// domain radius is distance to nearest side
//   domain.radius = (*it);
//   domain.edge = it - d;

//   return;

// }

void Device::MaximumSphere( Domain& domain)
{
  double lower = domain.runner[2] - r_min[2];
  double upper = r_max[2] - domain.runner[2];
  if( upper < lower )
  {
    domain.radius = upper;
    domain.edge = 5;
  }
  else
  {
    domain.radius = lower;
    domain.edge = 4;
  }
  return;
}

double Device::passage( Domain& domain )
{
  /// sampled boundary Voltage
  double vapp = 0.0;
  if( domain.radius <= m_skin )
  {
    /// wieghted (different eps) and sampled potential
    vapp = Vapp[domain.edge];
    /// runner now on surface
    domain.in_device = false;
  }
  return vapp;
}

double Device::OC_passage(const double V0,
                          Domain& domain ,
                          WOSParticles* WP)
{
  double vapp = 0.0;
  if( domain.radius <= m_skin )
  {
    /// wieghted (different eps) and sampled potential
    for(int i = 0; i < WP->m_qpts; i++)
      vapp += WP->Q[i]*(WP->wt[i]-1.0);
    vapp *= Vapp[domain.edge];
    vapp += V0; //(*sumQ) shouldn't this sum to zero for charge neutral system?
    domain.in_device = false;/// runner now on surface
  }
  return vapp;
}

double Device::contrib0( int i,
                         const Domain& domain,
                         WOSParticles* WP)
{
  double contrib = -WP->Q[i] / domain.radius;
  double Qi = WP->Q[i];
  WP->Q[i] = 0.0;
  for(int j = 0; j < WP->m_qpts; j++)
  {
    posvec_t rdiff = WP->R[j] - domain.runner;
    double r = sqrt(dot(rdiff,rdiff));
    /// r should be smaller than radius but not too close to centre.
    contrib += WP->Q[j] * ( std::abs( r - domain.radius ) - r + domain.radius )
               / ( 2.0 * r * domain.radius + 1.e-20 );
  }
  WP->Q[i] = Qi;
  return contrib;
}

double Device::contribk( const Domain& domain,
                         const WOSParticles* WP)
{
  double contrib = 0.0;
  for(int j = 0; j < WP->m_qpts; j++)
  {
    posvec_t rdiff = WP->R[j] - domain.runner;
    double r = sqrt(dot(rdiff,rdiff));
    contrib +=  WP->Q[j] * ( std::abs( r - domain.radius) - r + domain.radius )
                / ( 2.0 * r * domain.radius + 1.e-20 );
  }
  return contrib;
}

double Device::OC_contrib0( double d0,
                            const posvec_t& r1,
                            WOSParticles* WP)
{
  double contrib = 0.0;
  double vimg = 0.0;
  for(int i = 0; i < WP->m_qpts; i++)
  {
    posvec_t Ri = WP->R[i] - WP->R0;
    double risq = dot(Ri,Ri);
    double ri = sqrt(risq);
    double drsq = d0 * d0 - risq;
    /// the image contribution from the zeroth domain:
    vimg = -0.5 * WP->Q[i] * d0 / ( drsq + 1.e-20 );
    /// interparticle Goc
    double Goc = 0.0;
    for(int j = 0; j < i; j++)
    {
      posvec_t vec1 = WP->R[i] - WP->R[j];
      posvec_t vec2 = WP->R[j] - (d0 * d0 / risq ) * WP->R[i];
      double rij = sqrt(dot(vec1,vec1));
      double oc = ri * sqrt(dot(vec2,vec2)) + 1.e-20;
      Goc += WP->Q[j] * (1.0/rij - d0/oc);
    }
    contrib += WP->Q[i] * ( vimg + Goc );
    posvec_t vec3 = WP->R[i] - r1;
    double denom = dot(vec3,vec3);
    denom *= sqrt(denom);
    WP->wt[i] = d0 * drsq / denom;
  }
  return contrib;
}






// double Device::contrib0( int i,
// 			 const Domain& domain,
// 			 WOSParticles* WP){

//   double contrib = 0.0; ///-WP->Q[i] / domain.radius;
//   double Qi = WP->Q[i]; WP->Q[i] = 0.0;
//   for(int j = 0; j < WP->m_qpts; j++){
//     posvec_t rdiff = WP->R[j] - domain.runner;
//     double r = sqrt(dot(rdiff,rdiff));
//     /// r should be smaller than radius but not too close to centre.
//     contrib += WP->Q[j] / (r+1.e-20);
//     std::cout << "contrib: " << contrib << '\t' << WP->Q[j] / (r+1.e-20) << std::endl;
//   }
//   WP->Q[i] = Qi;
//   return contrib;

// }

