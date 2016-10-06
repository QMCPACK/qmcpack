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
    
    



#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "OhmmsPETE/OhmmsVector.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimCubicSpline.h"

#define FUNCTION(r) cos(r)
#define DFUNCTION(r) sin(r)
//#define USE_PBC

int main(int argc, char** argv)
{
  double ri = 0.0;
  double rf = 1.0;
  int npts = 101;
  const double nk=2;
  const double twopi = 8.0*atan(1.0)*nk;
  LinearGrid<double> agrid;
  agrid.set(ri,rf,npts);
#if defined(USE_PBC)
  typedef OneDimCubicSplinePBC<double> OrbitalType;
#else
  typedef OneDimCubicSplineFirst<double> OrbitalType;
#endif
  OrbitalType aorb(&agrid);
  aorb.resize(agrid.size());
  for(int i=0; i<agrid.size(); i++)
  {
    aorb(i)=FUNCTION(twopi*agrid(i));
  }
  aorb.spline(0,twopi*DFUNCTION(twopi*agrid.rmin()),agrid.size()-1, twopi*DFUNCTION(twopi*agrid.rmax()));
  std::string fname("testpbc.dat");
  std::ofstream dfile(fname.c_str());
  dfile.setf(std::ios::scientific, std::ios::floatfield);
  dfile.setf(std::ios::left,std::ios::adjustfield);
  dfile.precision(15);
  const double c1=1.0/twopi;
  const double c2=c1*c1;
  double du,d2u,_r,_rinv,y;
#if defined(USE_PBC)
  for(int ig=agrid.size()/2; ig<agrid.size() ; ig++)
  {
    _r = agrid(ig)+0.5*agrid.dr(ig)-agrid.rmax();
    _rinv = 1.0/_r;
    //aorb.setgrid(_r);
    y=aorb.evaluate(_r,_rinv,du,d2u);
    dfile << std::setw(30) << _r << std::setw(30) << FUNCTION(twopi*_r) << std::setw(30) << y << std::setw(30) << du*c1 << std::setw(3) << d2u*c2 << std::endl;
  }
#endif
  for(int ig=0; ig<agrid.size()-1; ig++)
  {
    _r = agrid(ig)+0.5*agrid.dr(ig);
    _rinv = 1.0/_r;
    //aorb.setgrid(_r);
    y=aorb.evaluate(_r,_rinv,du,d2u);
    dfile << std::setw(30) << _r << std::setw(30) << FUNCTION(twopi*_r) << std::setw(30) << y
          << std::setw(30) << du*c1 << std::setw(3) << d2u*c2 << std::endl;
  }
#if defined(USE_PBC)
  for(int ig=0; ig<agrid.size()/2; ig++)
  {
    _r = agrid(ig)+0.5*agrid.dr(ig)+agrid.rmax();
    _rinv = 1.0/_r;
    //aorb.setgrid(_r);
    y=aorb.evaluate(_r,_rinv,du,d2u);
    dfile << std::setw(30) << _r << std::setw(30) << FUNCTION(twopi*_r) << std::setw(30) << y << std::setw(30) << du*c1 << std::setw(3) << d2u*c2 << std::endl;
  }
#endif
  dfile << std::endl;
  return 0;
}
