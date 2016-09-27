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
    
    
#include "QMCHamiltonians/WOS/WOSPotential.h"

using namespace qmcplusplus;

/*
  Here we should have a number of experimental methods to implement
  WOS, and try them out.

  method0:   WOS basics
  method1:   WOS with antithetic variates
  method2:   WOS with a reference value
  method3:   WOS with importance sampled first sphere. : dipole


*/

WOSPotential::ValueType
WOSPotential::method0(ParticleSet& P)
{
  WP->setP(P);
  Domain domain;
  double pe = 0.0;
  double dpe = 0.0;
  double tau = 0.01;
  double v0 = -8.709;
  double branch = 0.0;
  for(int irun = 0; irun < m_runs; irun++)
  {
    double vrun = 0.0;
    for(int i = 0; i < WP->m_qpts; i++)
      /// start run
    {
      domain.runner = WP->R[i];
      domain.in_device = true;
      device->MaximumSphere(domain);
      double vi = device->contrib0(i,domain,WP);
      domain.WalkOnSphere();
      double vbare = device->passage(domain);
      while(domain.in_device)
      {
        device->MaximumSphere(domain);
        vi += device->contribk(domain,WP);
        domain.WalkOnSphere();
        vbare += device->passage(domain);
      }                                  /// end walk: runner on boundary
      vi = WP->Q[i] * ( 0.5 * vi + vbare );
      vrun += vi;
    }                                       /// end run/particle loop
    pe += vrun;
    dpe += vrun * vrun;         /// collect statistics
    branch += exp(-(vrun-v0)*tau);
  }                                         /// end runners-loop
  /// compute estimates
  pe *= m_norm;
  branch *= m_norm;
  dpe *= m_norm;
  double var = dpe - pe*pe;
  dpe = sqrt(m_norm * fabs ( dpe - pe * pe ));
  std::cout << "basic: " << pe << '\t' << dpe << '\t' << branch
       << '\t' << branch*exp(-0.5*var*tau*tau) << std::endl;
  exit(-1);
  return pe;
}

/// Antithetic variates
WOSPotential::ValueType
WOSPotential::method1(ParticleSet& P)
{
  WP->setP(P);
  Domain domain;
  double pe = 0.0;
  double dpe = 0.0;
  int hfruns = m_runs/2;
  double h_norm = 1.0/double(hfruns);
  for(int irun = 0; irun < hfruns; irun++)
  {
    double vrun = 0.0;
    for(int i = 0; i < WP->m_qpts; i++)
    {
      domain.runner = WP->R[i];
      domain.in_device = true;
      device->MaximumSphere(domain);
      double vi = device->contrib0(i,domain,WP);
      domain.WalkOnSphere();
      /// store the antithetic pair (refelcted about particle position)
      posvec_t r1(2*WP->R[i][0]-domain.runner[0],
                  2*WP->R[i][1]-domain.runner[1],
                  2*WP->R[i][2]-domain.runner[2]);
      double vbare = device->passage(domain);
      while(domain.in_device)
      {
        device->MaximumSphere(domain);
        vi += device->contribk(domain,WP);
        domain.WalkOnSphere();
        vbare += device->passage(domain);
      }
      vi = WP->Q[i] * ( 0.5 * vi + vbare );
      vrun += vi;  /// since there will be a pair of runners
      domain.runner = WP->R[i];
      domain.in_device = true;
      device->MaximumSphere(domain);
      vi = device->contrib0(i,domain,WP);
      //domain.WalkOnSphere();
      domain.runner = r1;
      vbare = device->passage(domain);
      while(domain.in_device)
      {
        device->MaximumSphere(domain);
        vi += device->contribk(domain,WP);
        domain.WalkOnSphere();
        vbare += device->passage(domain);
      }
      vi = WP->Q[i] * ( 0.5 * vi + vbare );
      vrun += vi;
    }                                            /// end run/particle loop
    vrun *= 0.5; /// since we used the antithetic pair.
    pe += vrun;
    dpe += vrun * vrun;
  }
  pe *= h_norm;
  dpe *= h_norm;
  dpe = sqrt(h_norm * std::abs( dpe - pe * pe ));
  //  std::cout << "antithetic: " << pe << '\t' << dpe << std::endl;
  //  exit(-1);
  return pe;
}


/// correlated sampling
WOSPotential::ValueType
WOSPotential::method2(ParticleSet& P)
{
  double V0 = 0;
  /// intialise the particles in WP;
  WP->setP(P);
  Domain domain;   /// create domain;
  double pe = 0.0;
  double dpe = 0.0;
  for(int irun = 0; irun < m_runs; irun++)
  {
    domain.runner = WP->R0;                /// initialise runner
    domain.in_device = true;               /// runner is inside device
    device->MaximumSphere(domain);         /// calc d_{0}
    domain.WalkOnSphere();
    double vD0 = device->OC_contrib0(domain.radius,domain.runner,WP);
    double vbare = device->OC_passage(V0,domain,WP);
    WP->calcwt();
    double vol = 0.0;
    while(domain.in_device)
    {
      device->MaximumSphere(domain);
      vol += device->contribk(domain,WP);
      domain.WalkOnSphere();
      vbare += device->OC_passage(V0,domain,WP);
    }
    vol *= WP->qwt;   /// the half has been included
    double vrun = vol + vbare + vD0;
    pe += vrun;
    dpe += vrun * vrun;
  }
  pe *= m_norm;
  dpe *= m_norm;
  dpe = ( dpe - pe * pe )/static_cast<double>(m_runs-1);
  /// CHANGE FOR DMC, WARNING Tau is zero for VMC WOS
  pe += dpe * Tau; // sigma^2 * Tau
  //  std::cout << "VWOS: "<< pe << '\t' << Tau << std::endl;
  P.Properties(WOSVAR) = -dpe;
  //P.Properties(WOSVAR) = dpe;
  return pe;
}



/// correlated sampling
WOSPotential::ValueType
WOSPotential::method6(ParticleSet& P)
{
  double V0 = 0;
  /// intialise the particles in WP;
  WP->setP(P);
  Domain domain;   /// create domain;
  double pe = 0.0;
  double dpe = 0.0;
  for(int irun = 0; irun < m_runs; irun++)
  {
    domain.runner = WP->R0;                /// initialise runner
    domain.in_device = true;               /// runner is inside device
    device->MaximumSphere(domain);         /// calc d_{0}
    domain.WalkOnSphere();
    double vD0 = device->OC_contrib0(domain.radius,domain.runner,WP);
    double vbare = device->OC_passage(V0,domain,WP);
    WP->calcwt();
    double vol = 0.0;
    while(domain.in_device)
    {
      device->MaximumSphere(domain);
      vol += device->contribk(domain,WP);
      domain.WalkOnSphere();
      vbare += device->OC_passage(V0,domain,WP);
    }
    vol *= WP->qwt;   /// the half has been included
    double vrun = vol + vbare + vD0;
    pe += vrun;
    dpe += vrun * vrun;
  }
  pe *= m_norm;
  dpe *= m_norm;
  dpe = ( dpe - pe * pe )/static_cast<double>(m_runs-1);
  double stsq = dpe * Tau * Tau;
  double correction = 1.0 + stsq * (0.25 + stsq /
                                    static_cast<double>(3*(m_runs+3))
                                   ) / static_cast<double>(2*(m_runs+1));
  /// CHANGE FOR DMC, WARNING Tau is zero for VMC WOS
  //  pe += dpe * Tau; // sigma^2 * Tau
  double pe1 = 0.0;
  double dpe1 = 0.0;
  for(int irun = 0; irun < 1; irun++)
  {
    domain.runner = WP->R0;                /// initialise runner
    domain.in_device = true;               /// runner is inside device
    device->MaximumSphere(domain);         /// calc d_{0}
    domain.WalkOnSphere();
    double vD0 = device->OC_contrib0(domain.radius,domain.runner,WP);
    double vbare = device->OC_passage(V0,domain,WP);
    WP->calcwt();
    double vol = 0.0;
    while(domain.in_device)
    {
      device->MaximumSphere(domain);
      vol += device->contribk(domain,WP);
      domain.WalkOnSphere();
      vbare += device->OC_passage(V0,domain,WP);
    }
    vol *= WP->qwt;   /// the half has been included
    double vrun = vol + vbare + vD0;
    pe1 += vrun;
    dpe1 += vrun * vrun;
  }
  //pe1 *= m_norm;
  P.Properties(WOSVAR) = 2*(pe - pe1)/Tau + dpe;// * correction;
  return pe1;
}

/// importance sampling
WOSPotential::ValueType
WOSPotential::method3(ParticleSet& P)
{
  WP->setP(P);
  Domain domain;
  double pe = 0.0;
  double dpe = 0.0;
  for(int irun = 0; irun < m_runs; irun++)
  {
    double vrun = 0.0;
    for(int i = 0; i < WP->m_qpts; i++)
      /// start run
    {
      domain.runner = WP->R[i];
      domain.in_device = true;
      device->MaximumSphere(domain);
      double vi = device->contrib0(i,domain,WP);
      double Gwt = domain.ImportanceWalk();
      double vbare = device->passage(domain);
      double vol = 0.0;
      while(domain.in_device)
      {
        device->MaximumSphere(domain);
        vol += device->contribk(domain,WP);
        domain.WalkOnSphere();
        vbare += device->passage(domain);
      }                                  /// end walk: runner on boundary
      //      vi = 0.0; vol = 0.0;
      vi = WP->Q[i] * (0.25*( 0.5 * vol + vbare )/Gwt + 0.5*vi);
      vrun += vi;
    }                                       /// end run/particle loop
    pe += vrun;
    dpe += vrun * vrun;         /// collect statistics
  }                                         /// end runners-loop
  /// compute estimates
  pe *= m_norm;
  dpe *= m_norm;
  dpe = sqrt(m_norm * fabs ( dpe - pe * pe ));
  double var = dpe - pe * pe;
  std::cout << "importance: " << pe << '\t' << dpe << '\t' << var << std::endl;
  exit(-1);
  return pe;
}

WOSPotential::ValueType
WOSPotential::method4(ParticleSet& P)
{
  double V0 = 0.0;
  WP->setP(P);
  Domain domain;
  double pe = 0.0;
  double dpe = 0.0;
  int hfruns = m_runs/2;
  for(int irun = 0; irun < hfruns; irun++)
  {
    domain.runner = WP->R0;                /// initialise runner
    domain.in_device = true;               /// runner is inside device
    device->MaximumSphere(domain);         /// calc d_{0}
    domain.WalkOnSphere();
    posvec_t r1(-domain.runner[0],-domain.runner[1],-domain.runner[2]);
    double vD0 = device->OC_contrib0(domain.radius,domain.runner,WP);
    double vbare = device->OC_passage(V0,domain,WP);
    WP->calcwt();
    double vol = 0.0;
    while(domain.in_device)
    {
      device->MaximumSphere(domain);
      vol += device->contribk(domain,WP);
      domain.WalkOnSphere();
      vbare += device->OC_passage(V0,domain,WP);
    }
    vol *= WP->qwt;   /// the half has been included
    double vrun = vol + vbare + vD0;
    pe += vrun;
    dpe += vrun * vrun;
    domain.runner = WP->R0;                /// initialise runner
    domain.in_device = true;               /// runner is inside device
    device->MaximumSphere(domain);         /// calc d_{0}
    domain.runner = r1;
    //domain.WalkOnSphere();
    vD0 = device->OC_contrib0(domain.radius,domain.runner,WP);
    vbare = device->OC_passage(V0,domain,WP);
    WP->calcwt();
    vol = 0.0;
    while(domain.in_device)
    {
      device->MaximumSphere(domain);
      vol += device->contribk(domain,WP);
      domain.WalkOnSphere();
      vbare += device->OC_passage(V0,domain,WP);
    }
    vol *= WP->qwt;   /// the half has been included
    vrun = vol + vbare + vD0;
    pe += vrun;
    dpe += vrun * vrun;
  }
  pe *= m_norm;
  dpe *= m_norm;
  dpe = sqrt(m_norm * fabs ( dpe - pe * pe ));
  std::cout << pe << '\t' << dpe << std::endl;
  exit(-1);
  return pe;
}

WOSPotential::ValueType
WOSPotential::method5(ParticleSet& P)
{
  static const double V0 = 0;
  /// intialise the particles in WP;
  WP->setP(P);
  Domain domain;   /// create domain;
  double branch = 0.0;
  double tau = 0.01;
  double v0 = -8.709;
  double pe = 0.0;
  double dpe = 0.0;
  for(int irun = 0; irun < m_runs; irun++)
  {
    domain.runner = WP->R0;                /// initialise runner
    domain.in_device = true;               /// runner is inside device
    device->MaximumSphere(domain);         /// calc d_{0}
    //    domain.WalkOnSphere();
    double Gwt = domain.ImportanceWalk();
    double vD0 = device->OC_contrib0(domain.radius,domain.runner,WP);
    double vbare = device->OC_passage(V0,domain,WP);
    WP->calcwt();
    double vol = 0.0;
    while(domain.in_device)
    {
      device->MaximumSphere(domain);
      vol += device->contribk(domain,WP);
      domain.WalkOnSphere();
      vbare += device->OC_passage(V0,domain,WP);
    }
    vol *= WP->qwt;   /// the half has been included
    //    double vrun = vol + vbare + vD0;
    double vrun = 0.25*(vol + vbare)/Gwt + vD0;
    pe += vrun;
    dpe += vrun * vrun;
    //    std::cout << vrun << '\t' << v0 << '\t' << branch << std::endl;
  }
  pe *= m_norm;
  branch *= m_norm;
  dpe *= m_norm;
  dpe = sqrt(m_norm * fabs ( dpe - pe * pe ));
  double var = dpe - pe * pe;
  std::cout << "correlated: " << pe << '\t' << dpe << '\t' << var << '\t'
       << branch << '\t' << branch*exp(-0.5*var*tau*tau) << std::endl;
  exit(-1);
  return pe;
}

