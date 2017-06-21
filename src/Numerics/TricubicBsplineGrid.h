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
    
    



#ifndef TRICUBIC_B_SPLINE_GRID_H
#define TRICUBIC_B_SPLINE_GRID_H

#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/Tensor.h"
#include "OhmmsPETE/OhmmsArray.h"
#include "Lattice/CrystalLattice.h"
#include "Numerics/GridBConds.h"
#include "type_traits/scalar_traits.h"
//#include <blitz/array.h>
//#include <blitz/tinymat.h>
//using namespace blitz;

namespace qmcplusplus
{

template<typename T, int BC0=NO_GBC, int BC1=NO_GBC, int BC2=NO_GBC>
struct TricubicBsplineGrid
{

  typedef typename scalar_traits<T>::real_type real_type;
  typedef typename scalar_traits<T>::value_type value_type;

  bool Interpolating, Periodic;
  // The grid sizes
  int Nx, Ny, Nz;
  //// The starting and ending values for the uniform grids
  //real_type xStart, xEnd, yStart, yEnd, zStart, zEnd;
  //// The box dimensions and their inverses
  //real_type Lx, LxInv, Ly, LyInv, Lz, LzInv;
  // The grid spacing and inverse
  real_type dx, dxInv, dy, dyInv, dz, dzInv;
  int ix0,ix1,ix2,ix3,iy0,iy1,iy2,iy3,iz0,iz1,iz2,iz3;

  GridBCond<real_type,BC0> gridX;
  GridBCond<real_type,BC1> gridY;
  GridBCond<real_type,BC2> gridZ;

  Tensor<real_type,4> A, dA, d2A, d3A;
  TinyVector<real_type,4> px, py, pz;
  TinyVector<real_type,4> a,b,c;
  TinyVector<real_type,4> da,db,dc;
  TinyVector<real_type,4> d2a,d2b,d2c;

  TricubicBsplineGrid();
  TricubicBsplineGrid<T,BC0,BC1,BC2>& operator=(const TricubicBsplineGrid<T,BC0,BC1,BC2>& rhs);

  //void setGrid(real_type xi, real_type xf, real_type yi, real_type yf,
  //    real_type zi, real_type zf, int nx, int ny, int nz,
  //    bool interp=true, bool periodic=true,bool openend=true);
  void setGrid(real_type xi, real_type xf, real_type yi, real_type yf,
               real_type zi, real_type zf, int nx, int ny, int nz,
               bool periodicx, bool periodicy, bool pbcz, bool openend=true);

  void setGrid(const GridBCond<real_type,BC0>& g0,
               const GridBCond<real_type,BC1>& g1,
               const GridBCond<real_type,BC2>& g2);

  bool Find(real_type x, real_type y, real_type z);

  bool FindAll(real_type x, real_type y, real_type z);

  T evaluate(const Array<T,3>& P) const ;

  T evaluate(const Array<T,3>& P, TinyVector<T,3>& grad, T& laplacian) const ;

  T evaluate(const Array<T,3>& P, TinyVector<T,3>& grad, Tensor<T,3>& secDerivs) const ;

  void Init(const Array<T,3>& data, Array<T,3>& P);
  void SolvePeriodicInterp (const Array<T,3> &data, Array<T,3>& P);
  void SolveFirstDerivInterp (const Array<T,3> &data, Array<T,3>& P);
  void MakePeriodic(Array<T,3>& P);

  /* return the distance between the center with PBC */
  inline real_type getSep2(real_type x, real_type y, real_type z) const
  {
    gridX.applyBC(x);
    gridY.applyBC(y);
    gridZ.applyBC(z);
    ////x-=nearbyint(x*LxInv)*Lx;
    ////y-=nearbyint(y*LyInv)*Ly;
    ////z-=nearbyint(z*LzInv)*Lz;
    //real_type x1=std::fmod(x*LxInv,1.0);
    //real_type y1=std::fmod(y*LyInv,1.0);
    //real_type z1=std::fmod(z*LzInv,1.0);
    //x=Lx*(x1-static_cast<int>(2.0*x1));
    //y=Ly*(y1-static_cast<int>(2.0*y1));
    //z=Lz*(z1-static_cast<int>(2.0*z1));
    return x*x+y*y+z*z;
  }
};

#include "Numerics/TricubicBsplineGrid.cpp"

// template<typename T>
//   struct TricubicBsplineTraits: public OrbitalTraits<T>
//   {
//     typedef typename OrbitalTraits<T>::real_type real_type;
//     typedef typename OrbitalTraits<T>::value_type value_type;
//     typedef TinyVector<real_type,3> PosType;
//     typedef TricubicBsplineGrid<T> GridType;
//     typedef Array<T,3>             StorageType;
//     int ObjectID;
//     real_type Rcut2;
//     Tensor<real_type,3> GGt;
//     CrystalLattice<real_type,3> Lattice;

//     inline void setRcut(real_type rc)
//     {
//       Rcut2=rc*rc;
//     }
//     virtual ~TricubicBsplineTraits() {}

//     void setLattice(const CrystalLattice<real_type,3>& lat)
//     {
//       Lattice.set(lat);
//       Lattice.print(std::cout);
//       GGt=dot(Lattice.G,transpose(Lattice.G));
//     }
//   };
}
#endif
