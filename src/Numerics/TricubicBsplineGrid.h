//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  Kenneth Esler 
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Modified by Jeongnim Kim for qmcpack
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef TRICUBIC_B_SPLINE_GRID_H
#define TRICUBIC_B_SPLINE_GRID_H

#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/Tensor.h"
#include "OhmmsPETE/OhmmsArray.h"
#include "QMCWaveFunctions/OrbitalTraits.h"
//#include <blitz/array.h>
//#include <blitz/tinymat.h>
//using namespace blitz;

namespace qmcplusplus {

  template<typename T>
  struct TricubicBsplineGrid: public OrbitalTraits<T> {

    typedef typename OrbitalTraits<T>::real_type real_type;
    typedef typename OrbitalTraits<T>::value_type value_type;

    bool Interpolating, Periodic;
    // The grid sizes
    int Nx, Ny, Nz;
    // The starting and ending values for the uniform grids
    real_type xStart, xEnd, yStart, yEnd, zStart, zEnd;
    // The box dimensions and their inverses
    real_type Lx, LxInv, Ly, LyInv, Lz, LzInv;
    // The grid spacing and inverse
    real_type dx, dxInv, dy, dyInv, dz, dzInv;
    int ix0,ix1,ix2,ix3,iy0,iy1,iy2,iy3,iz0,iz1,iz2,iz3;

    Tensor<real_type,4> A, dA, d2A, d3A;
    TinyVector<real_type,4> px, py, pz;
    TinyVector<real_type,4> a,b,c;
    TinyVector<real_type,4> da,db,dc;
    TinyVector<real_type,4> d2a,d2b,d2c;

    TricubicBsplineGrid();
    TricubicBsplineGrid(const TricubicBsplineGrid<T>& rhs);
    TricubicBsplineGrid& operator=(const TricubicBsplineGrid<T>& rhs);

    void setGrid(real_type xi, real_type xf, real_type yi, real_type yf,
        real_type zi, real_type zf, int nx, int ny, int nz, 
        bool interp=true, bool periodic=true,bool openend=true);

    void Find(real_type x, real_type y, real_type z); 

    void FindAll(real_type x, real_type y, real_type z);

    T evaluate(const Array<T,3>& P) const ;

    T evaluate(const Array<T,3>& P, TinyVector<T,3>& grad, T& laplacian) const ;

    T evaluate(const Array<T,3>& P, TinyVector<T,3>& grad, Tensor<T,3>& secDerivs) const ;

    void Init(const Array<T,3>& data, Array<T,3>& P);
    void SolvePeriodicInterp (const Array<T,3> &data, Array<T,3>& P);
    void MakePeriodic(Array<T,3>& P);

    /* return the distance between the center with PBC */
    inline T getSep2(real_type x, real_type y, real_type z)
    {
      x-=nearbyint(x*LxInv)*Lx;
      y-=nearbyint(y*LyInv)*Ly;
      z-=nearbyint(z*LzInv)*Lz;
      return x*x+y*y+z*z;
    }
  };

#include "Numerics/TricubicBsplineGrid.cpp"

  template<typename T>
    struct TricubicBsplineTraits: public OrbitalTraits<T>
    {
      typedef typename OrbitalTraits<T>::real_type real_type;
      typedef typename OrbitalTraits<T>::value_type value_type;
      typedef TinyVector<real_type,3> PosType;
      typedef TricubicBsplineGrid<T> GridType;
      typedef Array<T,3>             StorageType;
      int ObjectID;
      virtual ~TricubicBsplineTraits() {}
    };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
