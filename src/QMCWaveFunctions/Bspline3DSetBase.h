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
/** @file Bspline3DSet.h
 * @brief Define Bspline3DSetBase and its derived classes
 *
 * - Bspline3DSet_Ortho : orthorhombic unit cell
 * - Bspline3DSet_Gen : non-orthorhombic unit cell
 * - Bspline3DSet_Ortho_Trunc: orthorhombic unit cell with localized orbitals
 * - Bspline3DSet_Gen_Trunc : non-orthorhombic unit cell with localized orbitals
 */
#ifndef TRICUBIC_BSPLINE3D_SINGLEORBITALSET_BASE_H
#define TRICUBIC_BSPLINE3D_SINGLEORBITALSET_BASE_H

#include "Configuration.h"
#include "QMCWaveFunctions/SPOSetBase.h"
#include "Numerics/TricubicBsplineGrid.h"
#include "Optimize/VarList.h"

namespace qmcplusplus {

  /** base class for Bspline3DSet<bool ORTHO, bool TRUNC> and Bspline3DSet_XYZ
  */
#if defined(QMC_COMPLEX)
  struct Bspline3DSetBase: public OrbitalTraits<std::complex<OHMMS_PRECISION> >
#else
  struct Bspline3DSetBase: public OrbitalTraits<OHMMS_PRECISION>
#endif
  {

    typedef TinyVector<real_type,OHMMS_DIM>  PosType;
    typedef TinyVector<value_type,OHMMS_DIM> GradType;
    typedef TricubicBsplineGrid<value_type>  GridType;
    typedef Array<value_type,OHMMS_DIM>      StorageType;

    typedef SPOSetBase::IndexVector_t IndexVector_t;
    typedef SPOSetBase::ValueVector_t ValueVector_t;
    typedef SPOSetBase::ValueMatrix_t ValueMatrix_t;
    typedef SPOSetBase::GradVector_t  GradVector_t;
    typedef SPOSetBase::GradMatrix_t  GradMatrix_t;

    ///boolean 
    bool Orthorhombic;
    ///number of orbitals
    int NumOrbitals;
    ///index to keep track this object
    int ObjectID;
    ///cutoff radius
    real_type Rcut2;
    ///-|K|^2
    real_type mK2;
    ///TwistAngle in Cartesian Coordinate
    PosType TwistAngle;
    ///number of copies for each direction: only valid with localized orbitals
    TinyVector<int,OHMMS_DIM> Ncopy;
    ///metric tensor to handle generic unitcell
    Tensor<real_type,OHMMS_DIM> GGt;
    ///Lattice
    CrystalLattice<real_type,OHMMS_DIM> Lattice;
    ///grid
    GridType bKnots;
    ///centers
    vector<PosType> Centers;
    ///bspline data
    std::vector<const StorageType*> P;

    ///default constructor
    Bspline3DSetBase();

    ///virtual destructor
    virtual ~Bspline3DSetBase();

    inline void setRcut(real_type rc)
    {
      Rcut2=rc*rc;
    }

    /** set the lattice of the spline sets */
    void setLattice(const CrystalLattice<real_type,OHMMS_DIM>& lat);
    /** resize the internal storage P and Centers
     * @param norbs number of orbitals of this set
     */
    void resize(int norbs);
    /** set the twist angle 
     * @param tangle twist angle in Cartesian
     */
    void setTwistAngle(const PosType& tangle);
    /** set the grid
     * @param knots copy the grid
     */
    void setGrid(const GridType& knots);
    /** set the grid
     */
    void setGrid(real_type xi, real_type xf, 
        real_type yi, real_type yf, real_type zi, real_type zf, 
        int nx, int ny, int nz, 
        bool interp=true, bool periodic=true,bool openend=true);
    /** reset optimizable variables. 
     */
    void resetParameters(VarRegistry<real_type>& vlist);

    /** add a bspline orbital
     * @param i index of the orbital
     * @param data input data
     * @param curP interpolated data
     */
    void add(int i, const PosType& c, const StorageType& data, StorageType* curP);

    /** add a bspline orbital
     * @param i index of the orbital
     * @param curP interpolated data
     */
    void add(int i,const PosType& c,  StorageType* curP);

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2013 $   $Date: 2007-05-22 16:47:09 -0500 (Tue, 22 May 2007) $
 * $Id: TricubicBsplineSet.h 2013 2007-05-22 21:47:09Z jnkim $
 ***************************************************************************/
