//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#ifndef QMCPLUSPLUS_EINSPLINE_ORB_H
#define QMCPLUSPLUS_EINSPLINE_ORB_H

#include <Configuration.h>
#if defined(__xlC__)
#include <type_traits/scalar_traits.h>
#else
#include "QMCWaveFunctions/EinsplineWrapper.h"
#endif
#include "Numerics/HDFNumericAttrib.h"
#include "Lattice/CrystalLattice.h"
#include <cmath>

namespace qmcplusplus
{

#if defined(__xlC__)
template<typename T, unsigned D>
class EinsplineOrb
{
public:
  typedef typename scalar_traits<T>::real_type real_type;
  typedef TinyVector<real_type,D> PosType;
  PosType kVec;
  void evaluate (const PosType& u, T &psi) {}
  void evaluate (const PosType& u, T &psi, TinyVector<T,D> &grad, Tensor<T,D> &hess) {}
  void read (hid_t h5file, const std::string& groupPath)
  {
    APP_ABORT("Using xlC compiler. Cannot use EinsplineOrb");
  }
};


#else
template<typename T, int N>
class EinsplineOrb //: public QMCTraits
{
};




#if OHMMS_DIM==2
template<>
class EinsplineOrb<double,2> //: public QMCTraits
{
public:
  CrystalLattice<double,OHMMS_DIM> Lattice;

  typedef TinyVector<double,2> PosType;
  PosType Center, uCenter, uMin, uMax;
  // Reflection controls whether or not to reflect the orbital
  // across each axis.
  PosType Reflection;
  std::vector<PosType> uCenters, Reflections;
  double Radius, Energy;
  bool Localized;
  BsplineClass_2d_d *Bspline;
  // UBspline_2d_d *Spline;
  PosType kVec;

  inline void evaluate (const PosType& r, double &psi)
  {
    psi = (*Bspline)(r);
    // eval_UBspline_2d_d (Spline, r[0], r[1], &psi);
  }
  inline void evaluate (const PosType& r, double &psi,
                        TinyVector<double,2> &grad,
                        double &lapl)
  {
    Bspline->evaluate (r, psi, grad, lapl);
    // eval_UBspline_2d_d_vgl (Spline, r[0], r[1], &psi, &(grad[0]), &lapl);
  }
  inline void evaluate (const PosType& r, double &psi,
                        TinyVector<double,2> &grad,
                        Tensor<double,2> &hess)
  {
    Bspline->evaluate (r, psi, grad, hess);
//       eval_UBspline_2d_d_vgh (Spline, r[0], r[1], &psi,
// 			      &(grad[0]), &(hess(0,0)));
  }
  void read (hid_t h5file, std::string baseName)
  {
  }
  EinsplineOrb() : Center(PosType()), Radius(0.0), Energy(0.0),
    Localized(false), Reflection(1.0, 1.0)
  {
  }
};

template<>
class EinsplineOrb<std::complex<double>,2> //: public QMCTraits
{
public:
  CrystalLattice<double,OHMMS_DIM> Lattice;

  typedef TinyVector<double,2> PosType;
  PosType Center, uCenter, uMin, uMax;
  // Reflection controls whether or not to reflect the orbital
  // across each axis.
  PosType Reflection;
  std::vector<PosType> uCenters, Reflections;
  double Radius, Energy;
  bool Localized;
  BsplineClass_2d_z *Bspline;
  // UBspline_2d_z *Spline;
  PosType kVec;

  inline void evaluate (const PosType& r, std::complex<double> &psi)
  {
    psi = (*Bspline)(r);
//       eval_UBspline_2d_z (Spline, r[0], r[1], &psi);
  }
  inline void evaluate (const PosType& r, std::complex<double> &psi,
                        TinyVector<std::complex<double>,2> &grad,
                        std::complex<double> &lapl)
  {
    Bspline->evaluate (r, psi, grad, lapl);
//       eval_UBspline_2d_z_vgl (Spline, r[0], r[1],
// 			      &psi, &(grad[0]), &lapl);
  }
  inline void evaluate (const PosType& r, std::complex<double> &psi,
                        TinyVector<std::complex<double>,2> &grad,
                        Tensor<std::complex<double>,2> &hess)
  {
    Bspline->evaluate (r, psi, grad, hess);
//       eval_UBspline_2d_z_vgh (Spline, r[0], r[1],
// 			      &psi, &(grad[0]), &(hess(0,0)));
  }
  void read (hid_t h5file, const std::string& baseName)
  {
  }

  EinsplineOrb() :
    Center(PosType()), Radius(0.0), Energy(0.0), Localized(false),
    Reflection(1.0, 1.0)
  {
  }

};



#elif OHMMS_DIM==3
template<>
class EinsplineOrb<double,3> //: public QMCTraits
{
public:
  CrystalLattice<double,OHMMS_DIM> Lattice;

  typedef TinyVector<double,3> PosType;
  PosType Center, uCenter, uMin, uMax;
  // Reflection controls whether or not to reflect the orbital
  // across each axis.
  PosType Reflection;
  std::vector<PosType> uCenters, Reflections;
  double Radius, Energy;
  bool Localized;
  BsplineClass_3d_d *Bspline;
  // UBspline_3d_d *Spline;
  PosType kVec;

  inline void evaluate (const PosType& u, double &psi)
  {
    if (Localized)
    {
      PosType udiff = u - uCenter;
      udiff[0] -= round (udiff[0]);
      udiff[1] -= round (udiff[1]);
      udiff[2] -= round (udiff[2]);
      PosType rdiff = Lattice.toCart (udiff);
      if (dot (rdiff,rdiff) < Radius*Radius)
      {
        //for (int i=0; i<3; i++) {
        //  assert (udiff[i] > -0.5);
        //  assert (udiff[i] < 0.5);
        //}
        for (int i=0; i<3; i++)
          udiff[i] *= Reflection[i];
        udiff[0]+=0.5;
        udiff[1]+=0.5;
        udiff[2]+=0.5;
        psi = (*Bspline)(udiff);
        // eval_UBspline_3d_d (Spline, udiff[0], udiff[1], udiff[2], &psi);
      }
      else
        psi = 0.0;
    }
    else
      psi = (*Bspline)(u);
    // eval_UBspline_3d_d (Spline, u[0], u[1], u[2], &psi);
  }
  inline void evaluate (const PosType& u, double &psi, TinyVector<double,3> &grad,
                        Tensor<double,3> &hess)
  {
    if (Localized)
    {
      PosType udiff = u - uCenter;
      PosType urad  = udiff;
      udiff[0] -= round (udiff[0]);
      udiff[1] -= round (udiff[1]);
      udiff[2] -= round (udiff[2]);
      PosType rdiff = Lattice.toCart (udiff);
      if (dot (rdiff,rdiff) < Radius*Radius)
      {
        PosType uBox = uMax - uMin;
        for (int i=0; i<3; i++)
          udiff[i] *= Reflection[i];
        udiff[0]+=0.5;
        udiff[1]+=0.5;
        udiff[2]+=0.5;
        Bspline->evaluate (udiff, psi, grad, hess);
// 	  eval_UBspline_3d_d_vgh (Spline, udiff[0], udiff[1], udiff[2],
// 				  &psi, &(grad[0]), &(hess(0,0)));
        for (int i=0; i<3; i++)
        {
          grad[i] *= Reflection[i];
          for (int j=0; j<3; j++)
            hess(i,j) *= Reflection[i]*Reflection[j];
        }
      }
      else
      {
        psi = grad[0] = grad[1] = grad[2] = 0.0;
        hess(0,0) = hess(0,1) = hess(0,2) = 0.0;
        hess(1,0) = hess(1,1) = hess(1,2) = 0.0;
        hess(2,0) = hess(2,1) = hess(2,2) = 0.0;
      }
    }
    else
      Bspline->evaluate (u, psi, grad, hess);
// 	eval_UBspline_3d_d_vgh (Spline, u[0], u[1], u[2],
// 				&psi, &(grad[0]), &(hess(0,0)));
  }

  void read (hid_t h5file, const std::string& groupPath)
  {
    uMin   = PosType(0.0, 0.0, 0.0);
    uMax   = PosType(1.0, 1.0, 1.0);
    Center = PosType(0.5, 0.5, 0.5);
    // This controls the relative grid spacing at the edges versus
    // the center for nonuniform grids
    PosType clusterfactor(0.0, 0.0, 0.0);
    std::vector<PosType> centers;
    std::string centerName = groupPath + "center";
    std::string centersName = groupPath + "centers";
    std::string reflectName = groupPath + "reflections";
    std::string vectorName = groupPath + "eigenvector";
    std::string  valueName = groupPath + "eigenvalue";
    std::string radiusName = groupPath + "radius";
    std::string   uminName = groupPath + "umin";
    std::string   umaxName = groupPath + "umax";
    std::string xfactorName = groupPath + "xgrid/clusterfactor";
    std::string yfactorName = groupPath + "ygrid/clusterfactor";
    std::string zfactorName = groupPath + "zgrid/clusterfactor";
    HDFAttribIO<PosType> h_Center(Center);
    HDFAttribIO<PosType> h_uMin(uMin);
    HDFAttribIO<PosType> h_uMax(uMax);
    HDFAttribIO<double>  h_Radius(Radius);
    HDFAttribIO<double>  h_Energy(Energy);
    HDFAttribIO<std::vector<PosType> > h_Centers(centers);
    HDFAttribIO<std::vector<PosType> > h_Reflections (Reflections);
    HDFAttribIO<double>  h_xfactor (clusterfactor[0]);
    HDFAttribIO<double>  h_yfactor (clusterfactor[1]);
    HDFAttribIO<double>  h_zfactor (clusterfactor[2]);
    h_Center.read (h5file, centerName.c_str());
    h_Centers.read(h5file, centersName.c_str());
    uCenter = Lattice.toUnit (Center);
    h_Reflections.read (h5file, reflectName.c_str());
    h_Center.read  (h5file,  centerName.c_str());
    h_Centers.read (h5file, centersName.c_str());
    h_Radius.read  (h5file,  radiusName.c_str());
    h_Energy.read  (h5file,   valueName.c_str());
    h_uMin.read    (h5file,    uminName.c_str());
    h_uMax.read    (h5file,    umaxName.c_str());
    h_Radius.read  (h5file,  radiusName.c_str());
    h_xfactor.read (h5file, xfactorName.c_str());
    h_yfactor.read (h5file, yfactorName.c_str());
    h_zfactor.read (h5file, zfactorName.c_str());
    if (Reflections.size() > 0)
      Reflection = Reflections[0];
    else
      Reflection = PosType (1.0, 1.0, 1.0);
    // Localized = Radius > 0.0;
    Localized = (uMin[0] > 0.0);
    bool nonuniform = ((clusterfactor[0] > 0.0) &&
                       (clusterfactor[1] > 0.0) &&
                       (clusterfactor[2] > 0.0));
    uCenter = Lattice.toUnit (Center);
    uCenters.resize(centers.size());
    for (int i=0; i<centers.size(); i++)
      uCenters[i] = Lattice.toUnit (centers[i]);
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
    double dx, dy, dz;
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
      dx = (uMax[0] - uMin[0]) / (double)(nx-1);
      dy = (uMax[1] - uMin[1]) / (double)(ny-1);
      dz = (uMax[2] - uMin[2]) / (double)(nz-1);
    }
    else
    {
      uMin[0] = uMin[1] = uMin[2] = 0.0;
      uMax[0] = uMax[1] = uMax[2] = 1.0;
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
      dx = 1.0 / (double)nx;
      dy = 1.0 / (double)ny;
      dz = 1.0 / (double)nz;
    }
    if (Localized)
      fprintf(stderr, "  Center = (%8.5f, %8.5f %8.5f)   Radius = %8.5f  Mesh = %dx%dx%d\n",
              Center[0], Center[1], Center[2], Radius, nx, ny, nz);
    if (nonuniform)
      Bspline = new NUBsplineClass_3d_d (uMin, uMax, clusterfactor, xBC, yBC, zBC, realData);
    else
      Bspline = new UBsplineClass_3d_d (uMin, uMax, xBC, yBC, zBC, realData);
//       Spline = create_UBspline_3d_d (x_grid, y_grid, z_grid,
// 				     xBC, yBC, zBC, realData.data());
    // Now test spline to make sure it interpolates the data
    // for (int ix=0; ix<nx-1; ix++)
    //   for (int iy=0; iy<ny-1; iy++)
    //     for (int iz=0; iz<nz-1; iz++) {
    //       double ux = uMin[0] + dx*(double)ix;
    //       double uy = uMin[1] + dy*(double)iy;
    //       double uz = uMin[2] + dz*(double)iz;
    //       double val;
    //       eval_UBspline_3d_d (Spline, ux, uy, uz, &val);
    //       if (std::abs((val - realData(ix,iy,iz))) > 1.0e-12) {
    //         std::cerr << "Error in spline interpolation at ix=" << ix
    //	        << " iy=" << iy << " iz=" << iz << std::endl;
    //       }
    //    }
  }

  EinsplineOrb() :
    Center(PosType()), Radius(0.0), Energy(0.0), Localized(false),
    Reflection(1.0, 1.0, 1.0)
  {
  }

  EinsplineOrb(EinsplineOrb &orb)
  {
    Lattice   = orb.Lattice;
    Center    = orb.Center;
    uCenter   = orb.uCenter;
    uMin      = orb.uMin;
    uMax      = orb.uMax;
    uCenters.resize(orb.uCenters.size());
    uCenters  = orb.uCenters;
    Radius    = orb.Radius;
    Energy    = orb.Energy;
    Localized = orb.Localized;
    // Spline    = orb.Spline;
    Bspline      = orb.Bspline;
    kVec      = orb.kVec;
    Reflection = orb.Reflection;
    Reflections.resize(orb.Reflections.size());
    Reflections = orb.Reflections;
  }
};

template<>
class EinsplineOrb<std::complex<double>,3> //: public QMCTraits
{
public:
  CrystalLattice<double,OHMMS_DIM> Lattice;
  typedef TinyVector<double,3> PosType;
  PosType Center, uCenter, uMin, uMax;
  // Reflection controls whether or not to reflect the orbital
  // across each axis.
  PosType Reflection;
  std::vector<PosType> uCenters, Reflections;
  double Radius, Energy;
  bool Localized;
  BsplineClass_3d_z *Bspline;
  // UBspline_3d_z *Spline;
  PosType kVec;

  inline void evaluate (const PosType& u, std::complex<double> &psi)
  {
    if (Localized)
    {
      PosType udiff = u - uCenter;
      udiff[0] -= round (udiff[0]);
      udiff[1] -= round (udiff[1]);
      udiff[2] -= round (udiff[2]);
      PosType rdiff = Lattice.toCart (udiff);
      if (dot (rdiff,rdiff) <= Radius*Radius)
      {
        //for (int i=0; i<3; i++) {
        //  assert (udiff[i] > -0.5);
        //  assert (udiff[i] < 0.5);
        //}
        udiff = Reflection * udiff;
        udiff[0]+=0.5;
        udiff[1]+=0.5;
        udiff[2]+=0.5;
        psi = (*Bspline)(udiff);
        // eval_UBspline_3d_z (Spline, udiff[0], udiff[1], udiff[2], &psi);
      }
      else
        //psi = 1.0e-10;
        psi = std::complex<double>();
    }
    else
      psi = (*Bspline)(u);
// 	eval_UBspline_3d_z (Spline, u[0], u[1], u[2], &psi);
  }
  inline void evaluate (const PosType& u, std::complex<double> &psi,
                        TinyVector<std::complex<double>,3> &grad,
                        Tensor<std::complex<double>,3> &hess)
  {
    if (Localized)
    {
      PosType udiff = u - uCenter;
      udiff[0] -= round (udiff[0]);
      udiff[1] -= round (udiff[1]);
      udiff[2] -= round (udiff[2]);
      PosType rdiff = Lattice.toCart (udiff);
      if (dot (rdiff,rdiff) <= Radius*Radius)
      {
        // PosType uBox = uMax - uMin;
        // udiff[0]+=0.5;  udiff[1]+=0.5;  udiff[2]+=0.5;
        // 	  Bspline->evaluate (udiff, psi, grad, hess);
        // eval_UBspline_3d_z_vgh (Spline, udiff[0], udiff[1], udiff[2],
        //         		     &psi, &(grad[0]), &(hess(0,0)));
        PosType uBox = uMax - uMin;
        for (int i=0; i<3; i++)
          udiff[i] *= Reflection[i];
        udiff[0]+=0.5;
        udiff[1]+=0.5;
        udiff[2]+=0.5;
        Bspline->evaluate (udiff, psi, grad, hess);
// 	  eval_UBspline_3d_d_vgh (Spline, udiff[0], udiff[1], udiff[2],
// 				  &psi, &(grad[0]), &(hess(0,0)));
        for (int i=0; i<3; i++)
        {
          grad[i] *= Reflection[i];
          for (int j=0; j<3; j++)
            hess(i,j) *= Reflection[i]*Reflection[j];
        }
      }
      else
      {
        //psi = 1.0e-10;
        psi       = std::complex<double>();
        grad[0]   = grad[1]   = grad[2]   = std::complex<double>();
        hess(0,0) = hess(0,1) = hess(0,2) = std::complex<double>();
        hess(1,0) = hess(1,1) = hess(1,2) = std::complex<double>();
        hess(2,0) = hess(2,1) = hess(2,2) = std::complex<double>();
      }
    }
    else
      Bspline->evaluate (u, psi, grad, hess);
// 	eval_UBspline_3d_z_vgh (Spline, u[0], u[1], u[2],
// 				&psi, &(grad[0]), &(hess(0,0)));
  }
  void read (hid_t h5file, const std::string& groupPath)
  {
    uMin   = PosType(0.0, 0.0, 0.0);
    uMax   = PosType(1.0, 1.0, 1.0);
    Center = PosType(0.5, 0.5, 0.5);
    std::vector<PosType> centers;
    PosType clusterfactor (0.0, 0.0, 0.0);
    std::string centerName  = groupPath + "center";
    std::string centersName = groupPath + "centers";
    std::string reflectName = groupPath + "reflections";
    std::string vectorName  = groupPath + "eigenvector";
    std::string  valueName  = groupPath + "eigenvalue";
    std::string radiusName  = groupPath + "radius";
    std::string   uminName  = groupPath + "umin";
    std::string   umaxName  = groupPath + "umax";
    std::string xfactorName = groupPath + "xgrid/clusterfactor";
    std::string yfactorName = groupPath + "ygrid/clusterfactor";
    std::string zfactorName = groupPath + "zgrid/clusterfactor";
    HDFAttribIO<PosType> h_Center(Center);
    HDFAttribIO<PosType> h_uMin(uMin);
    HDFAttribIO<PosType> h_uMax(uMax);
    HDFAttribIO<double>  h_Radius(Radius);
    HDFAttribIO<double>  h_Energy(Energy);
    HDFAttribIO<std::vector<PosType> > h_Centers(centers);
    HDFAttribIO<std::vector<PosType> > h_Reflections(Reflections);
    HDFAttribIO<double>  h_xfactor (clusterfactor[0]);
    HDFAttribIO<double>  h_yfactor (clusterfactor[1]);
    HDFAttribIO<double>  h_zfactor (clusterfactor[2]);
    h_Center.read  (h5file,  centerName.c_str());
    h_Centers.read (h5file, centersName.c_str());
    h_Radius.read  (h5file,  radiusName.c_str());
    h_Energy.read  (h5file,   valueName.c_str());
    h_uMin.read    (h5file,    uminName.c_str());
    h_uMax.read    (h5file,    umaxName.c_str());
    h_Radius.read  (h5file,  radiusName.c_str());
    h_xfactor.read (h5file, xfactorName.c_str());
    h_yfactor.read (h5file, yfactorName.c_str());
    h_zfactor.read (h5file, zfactorName.c_str());
    // Localized = Radius > 0.0;
    Localized = (uMin[0] > 0.0);
    bool nonuniform = ((clusterfactor[0] > 0.0) &&
                       (clusterfactor[1] > 0.0) &&
                       (clusterfactor[2] > 0.0));
    uCenter = Lattice.toUnit (Center);
    uCenters.resize(centers.size());
    for (int i=0; i<centers.size(); i++)
      uCenters[i] = Lattice.toUnit (centers[i]);
    h_Reflections.read (h5file, reflectName.c_str());
    if (Reflections.size() > 0)
      Reflection = Reflections[0];
    else
      Reflection = PosType (1.0, 1.0, 1.0);
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
      if (nonuniform)
        Bspline = new NUBsplineClass_3d_z (uMin, uMax, clusterfactor, xBC, yBC, zBC, rawData);
      else
        Bspline = new UBsplineClass_3d_z (uMin, uMax, xBC, yBC, zBC, rawData);
// 	x_grid.start = uMin[0];  x_grid.end = uMax[0];  x_grid.num = nx;
// 	y_grid.start = uMin[1];  y_grid.end = uMax[1];  y_grid.num = ny;
// 	z_grid.start = uMin[2];  z_grid.end = uMax[2];  z_grid.num = nz;
// 	Spline = create_UBspline_3d_z (x_grid, y_grid, z_grid,
// 				       xBC, yBC, zBC, rawData.data());
    }
    else
    {
      Array<std::complex<double>,3> splineData(nx-1,ny-1,nz-1);
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
      TinyVector<double,3> start (0.0, 0.0, 0.0), end (1.0, 1.0, 1.0);
      Bspline = new UBsplineClass_3d_z (start, end, xBC, yBC, zBC,
                                        splineData);
// 	x_grid.start = 0.0;  x_grid.end = 1.0;  x_grid.num = nx-1;
// 	y_grid.start = 0.0;  y_grid.end = 1.0;  y_grid.num = ny-1;
// 	z_grid.start = 0.0;  z_grid.end = 1.0;  z_grid.num = nz-1;
// 	Spline = create_UBspline_3d_z (x_grid, y_grid, z_grid,
// 				       xBC, yBC, zBC, splineData.data());
    }
    if (Localized)
      fprintf(stderr, "  uCenter = (%8.5f, %8.5f %8.5f)   Radius = %8.5f  Mesh = %dx%dx%d\n",
              uCenter[0], uCenter[1], uCenter[2], Radius, nx, ny, nz);
  }

  EinsplineOrb() : Center(PosType()), Radius(0.0), Energy(0.0), Localized(false)
  {
  }
  EinsplineOrb(EinsplineOrb &orb)
  {
    Lattice     = orb.Lattice;
    Center      = orb.Center;
    uCenter     = orb.uCenter;
    uMin        = orb.uMin;
    uMax        = orb.uMax;
    uCenters.resize(orb.uCenters.size());
    uCenters    = orb.uCenters;
    Radius      = orb.Radius;
    Energy      = orb.Energy;
    Localized   = orb.Localized;
    // Spline    = orb.Spline;
    Bspline     = orb.Bspline;
    kVec        = orb.kVec;
    Reflection  = orb.Reflection;
    Reflections.resize(orb.Reflections.size());
    Reflections = orb.Reflections;
  }
};
#endif
#endif
}
#endif
