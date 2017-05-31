//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/** @file TricubicBsplineGridBC.h
 * @brief Definition of TricubicBsplineGridBC
 */
#ifndef TRICUBIC_B_SPLINE_GRID_BCTEMPLATED_H
#define TRICUBIC_B_SPLINE_GRID_BCTEMPLATED_H

#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/Tensor.h"
#include "OhmmsPETE/OhmmsArray.h"
#include "Lattice/CrystalLattice.h"
#include "Numerics/GridBConds.h"
#include "QMCWaveFunctions/OrbitalTraits.h"

namespace qmcplusplus
{

template<typename T, bool PBC0=true, bool PBC1=true, bool PBC2=true>
struct TricubicBsplineGridBC: public OrbitalTraits<T>
{

  typedef typename OrbitalTraits<T>::real_type real_type;
  typedef typename OrbitalTraits<T>::value_type value_type;

  bool Interpolating, Periodic;
  // The grid sizes
  int Nx, Ny, Nz;
  // The grid spacing and inverse
  real_type dx, dxInv, dy, dyInv, dz, dzInv;
  // index
  int ix0,ix1,ix2,ix3,iy0,iy1,iy2,iy3,iz0,iz1,iz2,iz3;

  GridBCond<real_type,PBC0> BC0;
  GridBCond<real_type,PBC1> BC1;
  GridBCond<real_type,PBC2> BC2;

  Tensor<real_type,4> A, dA, d2A, d3A;
  TinyVector<real_type,4> px, py, pz;
  TinyVector<real_type,4> a,b,c;
  TinyVector<real_type,4> da,db,dc;
  TinyVector<real_type,4> d2a,d2b,d2c;

  TricubicBsplineGridBC();
  TricubicBsplineGridBC(const TricubicBsplineGridBC<T,PBC0,PBC1,PBC2>& rhs);
  TricubicBsplineGridBC& operator=(const TricubicBsplineGridBC<T,PBC0,PBC1,PBC2>& rhs);

  void setGrid(real_type xi, real_type xf, real_type yi, real_type yf,
               real_type zi, real_type zf, int nx, int ny, int nz,
               bool openend);

  bool Find(real_type x, real_type y, real_type z);

  bool FindAll(real_type x, real_type y, real_type z);

  T evaluate(const Array<T,3>& P) const ;

  T evaluate(const Array<T,3>& P, TinyVector<T,3>& grad, T& laplacian) const ;

  T evaluate(const Array<T,3>& P, TinyVector<T,3>& grad, Tensor<T,3>& secDerivs) const ;

  void Init(const Array<T,3>& data, Array<T,3>& P);
  void SolvePeriodicInterp (const Array<T,3> &data, Array<T,3>& P);
  void SolveFirstDerivInterp (const Array<T,3> &data, Array<T,3>& P);
  void MakePeriodic(Array<T,3>& P);

  ///* return the distance between the center with PBC */
  inline real_type getSep2(real_type x, real_type y, real_type z)
  {
    BC0.applyBC(x);
    BC1.applyBC(y);
    BC2.applyBC(z);
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

#include "Numerics/TricubicBsplineGridBC.cpp"

}
#endif
