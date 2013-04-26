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
  string fname("testpbc.dat");
  std::ofstream dfile(fname.c_str());
  dfile.setf(ios::scientific, ios::floatfield);
  dfile.setf(ios::left,ios::adjustfield);
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
    dfile << setw(30) << _r << setw(30) << FUNCTION(twopi*_r) << setw(30) << y << setw(30) << du*c1 << setw(3) << d2u*c2 << endl;
  }
#endif
  for(int ig=0; ig<agrid.size()-1; ig++)
  {
    _r = agrid(ig)+0.5*agrid.dr(ig);
    _rinv = 1.0/_r;
    //aorb.setgrid(_r);
    y=aorb.evaluate(_r,_rinv,du,d2u);
    dfile << setw(30) << _r << setw(30) << FUNCTION(twopi*_r) << setw(30) << y
          << setw(30) << du*c1 << setw(3) << d2u*c2 << endl;
  }
#if defined(USE_PBC)
  for(int ig=0; ig<agrid.size()/2; ig++)
  {
    _r = agrid(ig)+0.5*agrid.dr(ig)+agrid.rmax();
    _rinv = 1.0/_r;
    //aorb.setgrid(_r);
    y=aorb.evaluate(_r,_rinv,du,d2u);
    dfile << setw(30) << _r << setw(30) << FUNCTION(twopi*_r) << setw(30) << y << setw(30) << du*c1 << setw(3) << d2u*c2 << endl;
  }
#endif
  dfile << endl;
  return 0;
}
