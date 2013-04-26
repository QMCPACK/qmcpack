/////////////////////////////////////////////////////////////////
// (c) Copyright 2007-  Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Modified by Jeongnim Kim for qmcpack
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file SplineGridHandler.h
 * @brief Define SplineGridHandler
 */
#ifndef QMCPLUSPLUS_SPLINE_GRID_HANDLER_H
#define QMCPLUSPLUS_SPLINE_GRID_HANDLER_H

#include "Lattice/CrystalLattice.h"

namespace qmcplusplus
{

/** handler of BsplineGrid */
template<typename T, unsigned D>
struct SplineGridHandler
{
  typedef CrystalLattice<T,D> unitcell_type;
  typedef TinyVector<T,D> pos_type;
  typedef Tensor<T,D>     tensor_type;
  ///boolean
  bool Orthorhombic;
  ///cutoff radius
  T Rcut;
  ///cutoff Rcut*Rcut
  T Rcut2;
  ///metric tensor to handle generic unitcell
  tensor_type GGt;
  /** Lattice for bspline grids
   */
  unitcell_type Lattice;
  ///unitcell for tiling
  unitcell_type UnitLattice;
  ///center in cartesian coordinate
  pos_type Center;
  ///lower corner in lattice
  pos_type Origin;

  SplineGridHandler(): Orthorhombic(true),
    Rcut(numeric_limits<T>::max()), Rcut2(numeric_limits<T>::max())
  {}

  ///destructor
  ~SplineGridHandler() {}

  /** set Rcut and Rcut2
   * @param rc cutoff radius
   */
  inline void setRcut(T rc)
  {
    Rcut=rc;
    Rcut2=rc*rc;
  }

  /** set the lattice of the spline sets
   * @param lat input lattice
   *
   * Both Lattice and UnitLattice is assigned to the same value
   */
  void setLattice(const unitcell_type& lat)
  {
    Lattice.set(lat);
    UnitLattice.set(lat);
    GGt=dot(lat.G,transpose(lat.G));
    Center=Lattice.toCart(pos_type(0.5));
    Origin=0.0;
    T offdiag=0.0;
    for(int idim=0; idim<D; idim++)
      for(int jdim=0; jdim<D; jdim++)
      {
        if(idim != jdim)
          offdiag+=abs(Lattice.R(idim,jdim));
      }
    Orthorhombic=(offdiag< numeric_limits<T>::epsilon());
  }
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2013 $   $Date: 2007-05-22 16:47:09 -0500 (Tue, 22 May 2007) $
 * $Id: TricubicBsplineSet.h 2013 2007-05-22 21:47:09Z jnkim $
 ***************************************************************************/
