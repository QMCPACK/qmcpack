#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "Numerics/TriCubicSplineT.h"

struct TestFunc {

  double k0,k1,k2;
  double d2factor;

  TestFunc(int nk0=1, int nk1=1, int nk2=1) {
    const double twopi = 8.0*atan(1.0);
    k0=twopi*static_cast<double>(nk0);
    k1=twopi*static_cast<double>(nk1);
    k2=twopi*static_cast<double>(nk2);
    d2factor = -(k0*k0+k1*k1+k2*k2);
  }

  //inline double operator()(const TinyVector<double,3>& pos) {
  //  return sin(k0*pos[0])*sin(k1*pos[1])*sin(k2*pos[2]);
  //}
  //inline double operator()(double x, double y, double z) {
  //  return sin(k0*x)*sin(k1*y)*sin(k2*z);
  //}
  //
  inline double f(const TinyVector<double,3>& pos) {
    return cos(k0*pos[0])*cos(k1*pos[1])*cos(k2*pos[2]);
  }

  inline double f(double x, double y, double z) {
    return cos(k0*x)*cos(k1*y)*cos(k2*z);
  }

  inline double d2f(const TinyVector<double,3>& pos) {
    return d2factor*f(pos);
  }

  inline double d2f(double x, double y, double z) {
    return d2factor*f(x,y,z);
  }

};

struct ComboFunc {

  std::vector<double> C;
  std::vector<TestFunc*> F;

  ComboFunc() {}
  void push_back(double c, TestFunc* fn) { C.push_back(c); F.push_back(fn);}


  inline double f(double x, double y, double z) {
    double res=0;
    for(int i=0; i<C.size(); i++) res += C[i]*F[i]->f(x,y,z);
    return res;
  }

  inline double d2f(double x, double y, double z) {
    double res=0;
    for(int i=0; i<C.size(); i++) res += C[i]*F[i]->d2f(x,y,z);
    return res;
  }

  inline double f(const TinyVector<double,3>& pos) {
    return f(pos[0],pos[1],pos[2]);
  }

  inline double d2f(const TinyVector<double,3>& pos) {
    return d2f(pos[0],pos[1],pos[2]);
  }

};

int main(int argc, char** argv) {

  double ri = 0.0;
  double rf = 1.0;
  int nptX = 101, nptY = 51, nptZ = 81;

  double xcut=0.23;
  double ycut=0.67;

  const int nk0=1;
  const int nk1=1;
  const int nk2=1;

  //Create one-dimensional grids for three orthogonal directions
  typedef LinearGrid<double> GridType;
  GridType gridX, gridY, gridZ;
  gridX.set(ri,rf,nptX);
  gridY.set(ri,rf,nptY);
  gridZ.set(ri,rf,nptZ);

  //Create XYZCubicGrid
  XYZCubicGrid<double> grid3(&gridX,&gridY,&gridZ);

  //Create a TriCubicSpline with PBC: have to think more about fixed-boundary conditions
  TriCubicSplineT<double> aorb(&grid3);

  //Create an analytic function for assignment
  ComboFunc infunc;
  infunc.push_back(0.5,new TestFunc(1,1,1));
  infunc.push_back(0.3,new TestFunc(1,1,2));
  infunc.push_back(0.1,new TestFunc(1,2,1));
  infunc.push_back(-0.3,new TestFunc(7,2,3));

  //Assign the values
  for(int ix=0; ix<nptX; ix++) {
    double x(gridX(ix));
    for(int iy=0; iy<nptY; iy++) {
      double y(gridY(iy));
      for(int iz=0; iz<nptZ; iz++) {
         aorb(ix,iy,iz)=infunc.f(x,y,gridZ(iz));
      }
    }
  }

  //Reset the coefficients
  aorb.reset();
 
  string fname("spline3d.dat");
  std::ofstream dfile(fname.c_str());
  dfile.setf(ios::scientific, ios::floatfield);
  dfile.setf(ios::left,ios::adjustfield);
  dfile.precision(15);

  double lap;
  TinyVector<double,3> grad;
  for(int k=0; k<nptY-1; k++) {
    //TinyVector<double,3> p(xcut,ycut,gridZ(k)+0.11*gridZ.dr(k));
    TinyVector<double,3> p(xcut,gridY(k)+0.11*gridY.dr(k),ycut);
    aorb.setgrid(p);
    double y=aorb.evaluate(p,grad,lap);
    dfile << setw(30) << p[1] << setw(30) << infunc.f(p) << setw(30) << y << setw(30) << infunc.d2f(p) << setw(30) << lap << endl;
  }

  return 0;
}
