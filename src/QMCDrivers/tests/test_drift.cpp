//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Particle/ParticleSet.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCDrivers/DriftOperators.h"

namespace qmcplusplus
{

TEST_CASE("drift pbyp and node correction real", "[drivers][drift]")
{
  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  MCWalkerConfiguration elec;

  elec.setName("elec");
  elec.setBoundBox(false);
  std::vector<int> agroup(1);
  agroup[0] = 1;
  elec.create(agroup);

  ParticleSet::RealType tau = 0.5;
  ParticleSet::RealType mass= 0.85;
  std::vector<ParticleSet::RealType> massinv(1,1./mass);
  ParticleSet::ParticlePos_t drift(1);

  // check from -xtot/2 to xtot/2 in step size of dx i.e. np.arange(-xtot/2,xtot/2,dx) 
  double xtot  = 10.;
  int    nx    = 100;
  double gradx = -xtot/2.;
  double dx    = xtot/nx;

  //app_log() << " begin printing" << std::endl;
  for (int ix=0;ix<nx;ix++)
  {
    elec.G[0][0] = gradx;
    setScaledDriftPbyPandNodeCorr(tau,massinv,elec.G,drift);
    double dval = drift[0][0]; 

    double scale_factor = (-1.+std::sqrt(1.+2.*gradx*gradx*tau/mass))/(gradx*gradx*tau/mass);
    REQUIRE( dval == Approx(scale_factor*gradx*tau/mass) );

    //app_log() << gradx << " " << dval << std::endl;
    gradx += dx;
  }
  //app_log() << " end printing." << std::endl;
}

#ifdef QMC_COMPLEX
TEST_CASE("drift pbyp and node correction complex", "[drivers][drift]")
{ // basically copy and pasted from real test, except "myi"
  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  MCWalkerConfiguration elec;

  elec.setName("elec");
  elec.setBoundBox(false);
  std::vector<int> agroup(1);
  agroup[0] = 1;
  elec.create(agroup);

  ParticleSet::RealType tau = 0.5;
  ParticleSet::RealType mass= 0.85;
  std::vector<ParticleSet::RealType> massinv(1,1./mass);
  ParticleSet::ParticlePos_t drift(1);

  // check from -xtot/2 to xtot/2 in step size of dx i.e. np.arange(-xtot/2,xtot/2,dx) 
  double xtot  = 10.;
  int    nx    = 100;
  double gradx = -xtot/2.;
  double dx    = xtot/nx;

  // imaginary component of wf gradient should NOT affect drift
  std::complex<double> myi(0,1.9);
  for (int ix=0;ix<nx;ix++)
  {
    elec.G[0][0] = gradx+myi;
    setScaledDriftPbyPandNodeCorr(tau,massinv,elec.G,drift);
    double dval = drift[0][0]; 

    double scale_factor = (-1.+std::sqrt(1.+2.*gradx*gradx*tau/mass))/(gradx*gradx*tau/mass);
    REQUIRE( dval == Approx(scale_factor*gradx*tau/mass) );

    gradx += dx;
  }
}
#endif

TEST_CASE("get scaled drift real", "[drivers][drift]")
{
  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  MCWalkerConfiguration elec;

  elec.setName("elec");
  elec.setBoundBox(false);
  std::vector<int> agroup(1);
  agroup[0] = 1;
  elec.create(agroup);

  ParticleSet::RealType tau = 0.5;
  ParticleSet::RealType mass= 0.85;
  std::vector<ParticleSet::RealType> massinv(1,1./mass);
  ParticleSet::PosType drift;

  // check from -xtot/2 to xtot/2 in step size of dx i.e. np.arange(-xtot/2,xtot/2,dx) 
  double xtot  = 10.;
  int    nx    = 100;
  double gradx = -xtot/2.;
  double dx    = xtot/nx;

  for (int ix=0;ix<nx;ix++)
  {
    elec.G[0][0] = gradx;
    getScaledDrift(tau/mass,elec.G[0],drift);
    double dval = drift[0]; 

    double scale_factor = (-1.+std::sqrt(1.+2.*gradx*gradx*tau/mass))/(gradx*gradx*tau/mass);
    REQUIRE( dval == Approx(scale_factor*gradx*tau/mass) );

    gradx += dx;
  }
}

#ifdef QMC_COMPLEX
TEST_CASE("get scaled drift complex", "[drivers][drift]")
{
  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  MCWalkerConfiguration elec;

  elec.setName("elec");
  elec.setBoundBox(false);
  std::vector<int> agroup(1);
  agroup[0] = 1;
  elec.create(agroup);

  ParticleSet::RealType tau = 0.5;
  ParticleSet::RealType mass= 0.85;
  std::vector<ParticleSet::RealType> massinv(1,1./mass);
  ParticleSet::PosType drift;

  // check from -xtot/2 to xtot/2 in step size of dx i.e. np.arange(-xtot/2,xtot/2,dx) 
  double xtot  = 10.;
  int    nx    = 100;
  double gradx = -xtot/2.;
  double dx    = xtot/nx;

  // imaginary component of wf gradient should NOT affect drift
  std::complex<double> myi(0,1.9);
  for (int ix=0;ix<nx;ix++)
  {
    elec.G[0][0] = gradx+myi;
    getScaledDrift(tau/mass,elec.G[0],drift);
    double dval = drift[0]; 

    double scale_factor = (-1.+std::sqrt(1.+2.*gradx*gradx*tau/mass))/(gradx*gradx*tau/mass);
    REQUIRE( dval == Approx(scale_factor*gradx*tau/mass) );

    gradx += dx;
  }
}
#endif

}

