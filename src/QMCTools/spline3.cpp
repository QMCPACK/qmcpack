#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include "OhmmsData/OhmmsElementBase.h"
#include "Numerics/HDFSTLAttrib.h"
#include "Numerics/TriCubicSplineT.h"
#include "Numerics/HDFTriCubicSpline.h"
#include "Utilities/Clock.h"

struct TestFunc
{

  double k0,k1,k2;
  double d2factor;

  TestFunc(int nk0=1, int nk1=1, int nk2=1)
  {
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
  inline double f(const TinyVector<double,3>& pos)
  {
    return sin(k0*pos[0])*sin(k1*pos[1])*sin(k2*pos[2]);
  }

  inline double f(double x, double y, double z)
  {
    return sin(k0*x)*sin(k1*y)*sin(k2*z);
  }

  inline double d2f(const TinyVector<double,3>& pos)
  {
    return d2factor*f(pos);
  }

  inline double d2f(double x, double y, double z)
  {
    return d2factor*f(x,y,z);
  }

};

struct ComboFunc
{

  std::vector<double> C;
  std::vector<TestFunc*> F;

  ComboFunc() {}
  void push_back(double c, TestFunc* fn)
  {
    C.push_back(c);
    F.push_back(fn);
  }


  inline double f(double x, double y, double z)
  {
    double res=0;
    for(int i=0; i<C.size(); i++)
      res += C[i]*F[i]->f(x,y,z);
    return res;
  }

  inline double d2f(double x, double y, double z)
  {
    double res=0;
    for(int i=0; i<C.size(); i++)
      res += C[i]*F[i]->d2f(x,y,z);
    return res;
  }

  inline double f(const TinyVector<double,3>& pos)
  {
    return f(pos[0],pos[1],pos[2]);
  }

  inline double d2f(const TinyVector<double,3>& pos)
  {
    return d2f(pos[0],pos[1],pos[2]);
  }

};

int main(int argc, char** argv)
{
  double ri = 0.0;
  double rf = 1.0;
  std::vector<int> npts(3);
  npts[0]=51;
  npts[1]=51;
  npts[2]=51;
  double xcut=0.23;
  double ycut=0.67;
  const int nk0=1;
  const int nk1=1;
  const int nk2=1;
  //Create one-dimensional grids for three orthogonal directions
  typedef LinearGrid<double> GridType;
  GridType gridX, gridY, gridZ;
  gridX.set(ri,rf,npts[0]);
  gridY.set(ri,rf,npts[1]);
  gridZ.set(ri,rf,npts[2]);
  //Create an analytic function for assignment
  ComboFunc infunc;
  infunc.push_back(0.5,new TestFunc(1,1,1));
  //infunc.push_back(0.3,new TestFunc(1,1,2));
  //infunc.push_back(0.1,new TestFunc(1,2,1));
  //infunc.push_back(0.01,new TestFunc(2,1,1));
  //infunc.push_back(0.01,new TestFunc(2,2,1));
  //infunc.push_back(0.001,new TestFunc(2,1,2));
  //infunc.push_back(0.001,new TestFunc(2,2,2));
  //infunc.push_back(0.001,new TestFunc(5,5,5));
  //infunc.push_back(-0.3,new TestFunc(7,2,3));
  //infunc.push_back(0.01,new TestFunc(7,7,7));
  //infunc.push_back(0.001,new TestFunc(5,5,5));
  //Write to an array
  std::vector<double> inData(npts[0]*npts[1]*npts[2]);
  std::vector<double>::iterator it(inData.begin());
  Pooma::Clock timer;
  timer.start();
  //Assign the values
  for(int ix=0; ix<npts[0]; ix++)
  {
    double x(gridX(ix));
    for(int iy=0; iy<npts[1]; iy++)
    {
      double y(gridY(iy));
      for(int iz=0; iz<npts[2]; iz++)
      {
        (*it)=infunc.f(x,y,gridZ(iz));
        ++it;
      }
    }
  }
  timer.stop();
  cout << "Time to evaluate " << timer.cpu_time() << endl;
  //Test TriCubicSplineT function
  //Create XYZCubicGrid
  XYZCubicGrid<double> grid3(&gridX,&gridY,&gridZ);
  //Create a TriCubicSpline with PBC: have to think more about fixed-boundary conditions
  TriCubicSplineT<double> aorb(&grid3);
  //Reset the coefficients
  aorb.reset(inData.begin(), inData.end());
  double lap,val;
  TinyVector<double,3> grad;
  //aorb.reset();
  //Write for vtk ImageData
  string fname("spline3d.vti");
  std::ofstream dfile(fname.c_str());
  dfile.setf(ios::scientific, ios::floatfield);
  dfile.setf(ios::left,ios::adjustfield);
  dfile.precision(10);
  dfile << "<?xml version=\"1.0\"?>" << endl;
  dfile << "<VTKFile type=\"ImageData\" version=\"0.1\">" << endl;
  dfile << "  <ImageData WholeExtent=\"0 " << npts[0]-2 << " 0 " << npts[1]-2 << " 0 " << npts[2]-2
        << "\" Origin=\"0 0 0\" Spacing=\"1 1 1\">"<< endl;
  dfile << "    <Piece Extent=\"0 " << npts[0]-2 << " 0 " << npts[1]-2 << " 0 " << npts[2]-2 << "\">" << endl;
  dfile << "       <PointData Scalars=\"wfs\">" << endl;
  dfile << "          <DataArray type=\"Float32\" Name=\"wfs\">" << endl;
  timer.start();
  int ng=0;
  for(int ix=0; ix<npts[0]-1; ix++)
  {
    double x(gridX(ix));
    for(int iy=0; iy<npts[1]-1; iy++)
    {
      double y(gridY(iy));
      for(int iz=0; iz<npts[2]-1; iz++, ng++)
      {
        TinyVector<double,3> p(x,y,gridZ(iz));
        //aorb.setgrid(p);
        //Timing with the ofstream is not correct.
        //Uncomment the line below and comment out the next two line.
        //double t=aorb.evaluate(p,grad,lap);
        dfile << setw(20) << aorb.evaluate(p,grad,lap);
        if(ng%5 == 4)
          dfile << endl;
      }
    }
  }
  timer.stop();
  cout << "Time to evaluate with spline " << timer.cpu_time() << endl;
  dfile << "          </DataArray>" << endl;
  dfile << "       </PointData>" << endl;
  dfile << "    </Piece>" << endl;
  dfile << "  </ImageData>" << endl;
  dfile << "</VTKFile>" << endl;
  hid_t h_file = H5Fcreate("spline3d.h5",H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  HDFAttribIO<std::vector<double> > dump(inData,npts);
  dump.write(h_file,"orb0000");
  HDFAttribIO<TriCubicSplineT<double> > dump1(aorb);
  dump1.write(h_file,"spline0000");
  H5Fclose(h_file);
  //double lap;
  //TinyVector<double,3> grad;
  //for(int k=0; k<nptY-1; k++) {
  //  //TinyVector<double,3> p(xcut,ycut,gridZ(k)+0.11*gridZ.dr(k));
  //  TinyVector<double,3> p(xcut,gridY(k)+0.11*gridY.dr(k),ycut);
  //  aorb.setgrid(p);
  //  double y=aorb.evaluate(p,grad,lap);
  //  dfile << setw(30) << p[1] << setw(30) << infunc.f(p) << setw(30) << y << setw(30) << infunc.d2f(p) << setw(30) << lap << endl;
  //}
  return 0;
}
