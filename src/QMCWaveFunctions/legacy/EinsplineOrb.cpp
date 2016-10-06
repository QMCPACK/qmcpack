//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "EinsplineOrb.h"
#include "Numerics/HDFNumericAttrib.h"

namespace qmcplusplus
{

void
EinsplineOrb<double,2>::read (hid_t h5file, std::string groupPath)
{
  std::cerr << "2D orbital reads not yet implemented.\n";
  abort();
}

void
EinsplineOrb<double,3>::read (hid_t h5file, std::string groupPath)
{
  std::string centerName = groupPath + "center";
  std::string vectorName = groupPath + "eigenvector";
  std::string  valueName = groupPath + "eigenvalue";
  std::string radiusName = groupPath + "radius";
  std::string   uminName = groupPath + "umin";
  std::string   umaxName = groupPath + "umax";
  HDFAttribIO<PosType> h_Center(Center);
  HDFAttribIO<PosType> h_uMin(uMin);
  HDFAttribIO<PosType> h_uMax(uMax);
  HDFAttribIO<double>  h_Radius(Radius);
  HDFAttribIO<double>  h_Energy(Energy);
  h_Center.read(h5file, centerName.c_str());
  h_Radius.read(h5file, radiusName.c_str());
  h_Energy.read(h5file,  valueName.c_str());
  Localized = Radius > 0.0;
  Array<std::complex<double>,3> rawData;
  Array<double,3> realData;
  HDFAttribIO<Array<std::complex<double>,3> > h_rawData(rawData);
  h_rawData.read(h5file, vectorName.c_str());
  int nx, ny, nz;
  nx = rawData.size(0);
  ny=rawData.size(1);
  nz=rawData.size(2);
  if (!Localized)
  {
    nx--;
    ny--;
    nz--;
  }
  realData.resize(nx,ny,nz);
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++)
        realData(ix,iy,iz) = rawData(ix,iy,iz).real();
  Ugrid x_grid, y_grid, z_grid;
  BCtype_d xBC, yBC, zBC;
  if (Localized)
  {
    xBC.lCode = NATURAL;
    xBC.rCode = NATURAL;
    yBC.lCode = NATURAL;
    yBC.rCode = NATURAL;
    zBC.lCode = NATURAL;
    zBC.rCode = NATURAL;
    x_grid.start = uMin[0];
    x_grid.end = uMax[0];
    x_grid.num = nx;
    y_grid.start = uMin[1];
    y_grid.end = uMax[1];
    y_grid.num = ny;
    z_grid.start = uMin[2];
    z_grid.end = uMax[2];
    z_grid.num = nz;
  }
  else
  {
    xBC.lCode = PERIODIC;
    xBC.rCode = PERIODIC;
    yBC.lCode = PERIODIC;
    yBC.rCode = PERIODIC;
    zBC.lCode = PERIODIC;
    zBC.rCode = PERIODIC;
    x_grid.start = 0.0;
    x_grid.end = 1.0;
    x_grid.num = nx;
    y_grid.start = 0.0;
    y_grid.end = 1.0;
    y_grid.num = ny;
    z_grid.start = 0.0;
    z_grid.end = 1.0;
    z_grid.num = nz;
  }
  if (Localized)
    fprintf (stderr, "  Center = (%8.5f, %8.5f %8.5f)   Radius = %8.5f  Mesh = %dx%dx%d\n",
             Center[0], Center[1], Center[2], Radius, nx, ny, nz);
  Spline = create_UBspline_3d_d (x_grid, y_grid, z_grid,
                                 xBC, yBC, zBC, &realData(0,0,0));
}

void
EinsplineOrb<std::complex<double>,2>::read (hid_t h5file, std::string groupPath)
{
  std::cerr << "2D orbital reads not yet implemented.\n";
  abort();
}

void
EinsplineOrb<std::complex<double>,3>::read (hid_t h5file, std::string groupPath)
{
  uMin   = PosType(0.0, 0.0, 0.0);
  uMax   = PosType(1.0, 1.0, 1.0);
  Center = PosType(0.5, 0.5, 0.5);
  std::string centerName = groupPath + "center";
  std::string vectorName = groupPath + "eigenvector";
  std::string  valueName = groupPath + "eigenvalue";
  std::string radiusName = groupPath + "radius";
  std::string   uminName = groupPath + "umin";
  std::string   umaxName = groupPath + "umax";
  HDFAttribIO<PosType> h_Center(Center);
  HDFAttribIO<PosType> h_uMin(uMin);
  HDFAttribIO<PosType> h_uMax(uMax);
  HDFAttribIO<double>  h_Radius(Radius);
  HDFAttribIO<double>  h_Energy(Energy);
  h_Center.read(h5file, centerName.c_str());
  uCenter = Lattice.toUnit (Center);
  h_Radius.read(h5file, radiusName.c_str());
  h_Energy.read(h5file,  valueName.c_str());
  h_uMax.read(h5file, umaxName.c_str());
  h_Radius.read(h5file, radiusName.c_str());
  Localized = Radius > 0.0;
  Array<std::complex<double>,3> rawData;
  HDFAttribIO<Array<std::complex<double>,3> > h_rawData(rawData);
  h_rawData.read(h5file, vectorName.c_str());
  int nx, ny, nz;
  nx = rawData.size(0);
  ny=rawData.size(1);
  nz=rawData.size(2);
  Ugrid x_grid, y_grid, z_grid;
  BCtype_z xBC, yBC, zBC;
  if (Localized)
  {
    xBC.lCode = NATURAL;
    xBC.rCode = NATURAL;
    yBC.lCode = NATURAL;
    yBC.rCode = NATURAL;
    zBC.lCode = NATURAL;
    zBC.rCode = NATURAL;
    x_grid.start = uMin[0];
    x_grid.end = uMax[0];
    x_grid.num = nx;
    y_grid.start = uMin[1];
    y_grid.end = uMax[1];
    y_grid.num = ny;
    z_grid.start = uMin[2];
    z_grid.end = uMax[2];
    z_grid.num = nz;
    Spline = create_UBspline_3d_z (x_grid, y_grid, z_grid, xBC, yBC, zBC, rawData.data());
  }
  else
  {
    Array<std::complex<double>,3> splineData;
    splineData.resize(nx-1,ny-1,nz-1);
    for (int ix=0; ix<nx-1; ix++)
      for (int iy=0; iy<ny-1; iy++)
        for (int iz=0; iz<nz-1; iz++)
          splineData(ix,iy,iz) = rawData(ix,iy,iz);
    xBC.lCode = PERIODIC;
    xBC.rCode = PERIODIC;
    yBC.lCode = PERIODIC;
    yBC.rCode = PERIODIC;
    zBC.lCode = PERIODIC;
    zBC.rCode = PERIODIC;
    x_grid.start = 0.0;
    x_grid.end = 1.0;
    x_grid.num = nx-1;
    y_grid.start = 0.0;
    y_grid.end = 1.0;
    y_grid.num = ny-1;
    z_grid.start = 0.0;
    z_grid.end = 1.0;
    z_grid.num = nz-1;
    Spline = create_UBspline_3d_z (x_grid, y_grid, z_grid, xBC, yBC, zBC, splineData.data());
  }
}
}

