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
#include <sstream>
#include "OhmmsData/OhmmsElementBase.h"
#include "Numerics/HDFSTLAttrib.h"
#include "Numerics/TriCubicSplineT.h"
#include "Numerics/HDFTriCubicSpline.h"
#include "Utilities/Clock.h"


int main(int argc, char** argv)
{
  Pooma::Clock timer;
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
  //Write to an array
  std::vector<double> inData(npts[0]*npts[1]*npts[2]);
  timer.start();
  hid_t h_file = H5Fopen("spline3d.h5",H5F_ACC_RDWR,H5P_DEFAULT);
  HDFAttribIO<std::vector<double> > dummy(inData,npts);
  dummy.read(h_file,"orb0000");
  H5Fclose(h_file);
  timer.stop();
  std::cout << "Time to read image data " << timer.cpu_time() << std::endl;
  //Create XYZCubicGrid
  XYZCubicGrid<double> grid3(&gridX,&gridY,&gridZ);
  //Create a TriCubicSpline with PBC: have to think more about fixed-boundary conditions
  TriCubicSplineT<double> aorb(&grid3);
  //Reset the coefficients
  timer.start();
  aorb.reset(inData.begin(), inData.end());
  timer.stop();
  std::cout << "Time to set up spline coefficients " << timer.cpu_time() << std::endl;
  double lap,val;
  TinyVector<double,3> grad;
  //vector<double>::iterator it(inData.begin());
  timer.start();
  for(int ix=0; ix<npts[0]-1; ix++)
  {
    double x(gridX(ix));
    for(int iy=0; iy<npts[1]-1; iy++)
    {
      double y(gridY(iy));
      int ng = npts[2]*(iy+ix*npts[1]);
      for(int iz=0; iz<npts[2]-1; iz++)
      {
        TinyVector<double,3> p(x,y,gridZ(iz));
        //aorb.setgrid(p);
        inData[ng++]=aorb.evaluate(p,grad,lap);
        //(*it) = aorb.evaluate(p,grad,lap); ++it;
      }
    }
  }
  timer.stop();
  std::cout << "Time to evaluate the values  " << timer.cpu_time() << std::endl;
  timer.start();
  h_file = H5Fcreate("spline3d_writeback.h5",H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  HDFAttribIO<std::vector<double> > dump(inData,npts);
  dump.write(h_file,"orb0000");
  HDFAttribIO<TriCubicSplineT<double> > dump1(aorb);
  dump1.write(h_file,"spline0000");
  H5Fclose(h_file);
  timer.stop();
  std::cout << "Time to write to hdf5 " << timer.cpu_time() << std::endl;
  return 0;
}
