//////////////////////////////////////////////////////////////////
// (c) Copyright 1998- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_CRYSTALLATTICE_H
#define OHMMS_CRYSTALLATTICE_H
#include <math.h>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/Tensor.h"
#include "Lattice/ParticleBConds.h"

/**@file CrystalLattice.h
 *@brief Declaration of CrystalLattice<T,D>
 */

#ifndef TWOPI
#ifndef M_PI
#define TWOPI 6.3661977236758134308E-1
#else
#define TWOPI 2*M_PI
#endif /* M_PI */
#endif /* TWOPI */

/** class to assist copy and unit conversion operations on position vectors
 */
struct PosUnit {

  /** enumeraton for the unit of position types.
   */
  enum {CartesianUnit=0,/*!< indicates that the values are in Cartesian units*/
        LatticeUnit/*!< indicates that the values are in Lattice units*/
       };
};

/** a class that defines a supercell in D-dimensional Euclean space.
 *
 *CrystalLattice handles the physical properties of a supercell, such
 *as lattice vectors, reciprocal vectors and metric tensors and provides
 *interfaces to access the lattice properties and convert units of
 *position vectors or a single-particle position from Cartesian to
 *Lattice Unit vice versa.
 *
 *The indices for R, G and D are chosen to perform
 *expression template operations with variable-cell algorithms.
 *
 */
template<class T, unsigned D>
struct CrystalLattice{

  ///enumeration for the dimension of the lattice
  enum {DIM = D};
  //@{
  ///the type of scalar
  typedef T                            Scalar_t;
  ///the type of a D-dimensional position vector 
  typedef TinyVector<T,D>              SingleParticlePos_t;
  ///the type of a D-dimensional index vector 
  typedef TinyVector<int,D>            SingleParticleIndex_t;
  ///the type of a D-dimensional Tensor
  typedef Tensor<T,D>                  Tensor_t;
  //@}

  //@{ 
  /**@brief Physcial properties of a supercell*/
  /// Volume of a supercell
  T  Volume; 
  ///Real-space unit vectors. R(i,j) i=vector and j=x,y,z
  Tensor_t R;
  ///Reciprocal unit vectors. G(j,i) i=vector and j=x,y,z
  Tensor_t G;
  ///Metric tensor
  Tensor_t M;
  ///Metric tensor for G vectors
  Tensor_t Mg;
  /**@brief Real-space unit vectors. 
   *
   *Introduced to efficiently return one vector at a time.
   *Rv[i] is D-dim vector of the ith direction. 
   */
  TinyVector<SingleParticlePos_t,D> Rv;
  /**@brief Reciprocal unit vectors. 
   *
   *Introduced to efficiently return one vector at a time.
   *Gv[i] is D-dim vector of the ith direction.
   */
  TinyVector<SingleParticlePos_t,D> Gv;
  //@}

  //@{
  /**@brief Parameters defining boundary Conditions */
  ///Functors that apply boundary conditions on the position vectors
  ParticleBConds<T,D> BConds;
  ///The boundary condition in each direction.
  TinyVector<int,D> BoxBConds;
  //@}

  ///default constructor, assign a huge supercell
  CrystalLattice();
  /** copy constructor
      @param rhs An existing SC object is copied to this SC.
  */
  CrystalLattice(const CrystalLattice<T,D>& rhs);
  
  ///destructor
  virtual ~CrystalLattice(){ }

  /**@param i the index of the directional vector, \f$i\in [0,D)\f$
   *@return The lattice vector of the ith direction
   *@brief Provide interfaces familiar to fotran users
   */
  inline SingleParticlePos_t a(int i) const {
    return Rv[i];
  }

  /**@param i the index of the directional vector, \f$i\in [0,D)\f$
   *@return The reciprocal vector of the ith direction
   *@brief Provide interfaces familiar to fotran users
   */
  inline SingleParticlePos_t b(int i) const {
    return Gv[i];
  }

  /** Convert a cartesian vector to a unit vector.
   * Boundary conditions are not applied.
   */
  inline SingleParticlePos_t toUnit(const SingleParticlePos_t &r) const {
    return dot(r,G);
  }

  /** Convert a unit vector to a cartesian vector.
   * Boundary conditions are not applied.
   */
  inline SingleParticlePos_t toCart(const SingleParticlePos_t &c) const {
    return dot(c,R);
  }

  /** evaluate the cartesian distance
   *@param ra a vector in the supercell unit
   *@param rb a vector in the supercell unit
   *@return Cartesian distance with two vectors in SC unit
   *
   @note The distance between two cartesian vectors are handled 
   *by dot function defined in OhmmsPETE/TinyVector.h
   */
  inline T Dot(const SingleParticlePos_t &ra, 
	       const SingleParticlePos_t &rb) const {
    return dot(ra,dot(M,rb));
  }

  /** conversion of a reciprocal-vector 
   *@param kin an input reciprocal vector in the Reciprocal-vector unit
   *@return k(reciprocal vector) in cartesian unit
  */
  inline SingleParticlePos_t k_cart(const SingleParticlePos_t& kin) const {
    return TWOPI*dot(G,kin);
  }

  /** evaluate \f$k^2\f$
   *
   *@param kin an input reciprocal vector in reciprocal-vector unit
   *@return \f$k_{in}^2\f$ 
   */
  inline T ksq(const SingleParticlePos_t& kin) const {
    return dot(kin,dot(Mg,kin));
  }

  ///assignment operator
  CrystalLattice<T,D>& operator=(const CrystalLattice<T,D>& rhs);

  /** assignment operator
   *@param rhs a tensor representing a unit cell
   */
  CrystalLattice<T,D>& operator=(const Tensor<T,D>& rhs);

  /** scale the lattice vectors by sc. All the internal data are reset.
   *@param sc the scaling value
   *@return a new CrystalLattice
   */
  CrystalLattice<T,D>& operator*=(T sc);
  
  /** set the lattice vector from the command-line options
   *@param argc the number of arguments
   *@param argv the argument lists
   *
   *This function is to provide a simple interface for testing.
   */
  void set(int argc, char **argv);
  
  /** set the lattice vector from the command-line options stored in a vector
   *@param argv the argument lists
   *
   *This function is to provide a simple interface for testing.
   */
  void set(vector<string>& argv);

  /** set the lattice vector by an array containing DxD T
   *@param sc a scalar to scale the input lattice parameters
   *@param lat the starting address of DxD T-elements representing a supercell
   */
  void set(T sc, T* lat= NULL);

  /** set the lattice vector by a CrystalLattice and expand it by integers
   *@param oldlat An input supercell to be copied.
   *@param uc An array to expand a supercell.
   */
  void set(const CrystalLattice<T,D>& oldlat, int* uc= NULL);

  /** set the lattice vector from the command-line options
   *@param lat a tensor representing a supercell
   */
  void set(const Tensor<T,D>& lat);

  /** Evaluate the reciprocal vectors, volume and metric tensor
   */
  void reset();

  //@{
  /* Copy functions with unit conversion*/
  template<class PA> void convert(const PA& pin, PA& pout) const;
  template<class PA> void convert2Unit(const PA& pin, PA& pout) const;
  template<class PA> void convert2Cart(const PA& pin, PA& pout) const;
  template<class PA> void convert2Unit(PA& pout) const;
  template<class PA> void convert2Cart(PA& pout) const;
  //@}

  template<class PA> void applyBC(const PA& pin, PA& pout, T del=1.e-6) const;
  template<class PA> void applyBC(PA& pos, T del=1.e-6) const;

  //! Print out CrystalLattice Data
  void print(ostream& , int level=2) const;
};

//including the definitions of the member functions
#include "Lattice/CrystalLattice.cpp"

/** Copy operation pout = pin with the unit checking of the in/out vectors
 *
 *@param pin an input position array
 *@param pout an output position array
 *
 *@note The units of in/out vectors are conserved.
 */
template<class T, unsigned D>
template<class PA>
inline void CrystalLattice<T,D>::convert(const PA& pin, PA& pout) const
{
  if(pin.getUnit() == pout.getUnit())   {
    pout = pin;
    return;
  }
  if(pin.getUnit() == PosUnit::LatticeUnit) {
    //convert to CartesianUnit
    for(int i=0; i<pin.size(); i++) pout[i] = dot(pin[i],R);
  } else {
    //convert to LatticeUnit
    for(int i=0; i<pin.size(); i++) pout[i] = dot(pin[i],G);
  }
}

/** Copy operation pout = pin with the unit checking of the in/out vectors
 *
 *@param pin an input position array
 *@param pout an output position array
 *
 *@note The unit of the out vector is forced to be Cartesian.
 */
template<class T, unsigned D>
template<class PA>
inline void CrystalLattice<T,D>::convert2Cart(const PA& pin, PA& pout) const
{
  pout.setUnit(PosUnit::CartesianUnit);
  if(pin.getUnit() == PosUnit::CartesianUnit) 
    pout = pin;
  else //need to convert to CartesianUnit
    for(int i=0; i<pin.size(); i++)  pout[i] = dot(pin[i],R);
}

/** Copy operation pout = pin with the unit checking of the in/out vectors
 *
 *@param pin an input position array
 *@param pout an output position array
 *
 *@note The unit of the out vector is forced to be Lattice unit.
 */
template<class T, unsigned D>
template<class PA>
inline void CrystalLattice<T,D>::convert2Unit(const PA& pin, PA& pout) const
{
  pout.setUnit(PosUnit::LatticeUnit);
  if(pin.getUnit() == PosUnit::LatticeUnit) 
    pout = pin;
  else //need to convert to LatticeUnit
    for(int i=0; i<pin.size(); i++)  pout[i] = dot(pin[i],G);
}

/** Conversion of a position array to Cartesian unit.
 *
 *@param pos in/out position array
 *
 *@note The unit of the out vector is forced to be Cartesian. If the unit
 *of pos is in Cartesian, no action is taken.
 */
template<class T, unsigned D>
template<class PA>
inline void CrystalLattice<T,D>::convert2Cart(PA& pos) const
{
  if(pos.getUnit() == PosUnit::LatticeUnit) {
    for(int i=0; i<pos.size(); i++)  pos[i] = dot(pos[i],R);
    pos.setUnit(PosUnit::CartesianUnit);
  }
}

/** Conversion of a position array to Lattice unit.
 *
 *@param pos in/out position array
 *
 *@note The unit of the out vector is forced to be in Lattice unit. 
 *If the unit of pos is in Lattice unit, no action is taken.
 */
template<class T, unsigned D>
template<class PA>
inline void CrystalLattice<T,D>::convert2Unit(PA& pos) const
{
  if(pos.getUnit() == PosUnit::CartesianUnit) {
    for(int i=0; i<pos.size(); i++)  pos[i] = dot(pos[i],G);
    pos.setUnit(PosUnit::LatticeUnit);
  }
}

/** Copy an in vector to an out vector after applying boundary conditions.
 *
 *@param pin an input position array
 *@param pout an output position array
 *@param del a tolerance for wrapping the values within [0,1)
 *
 *@note The values of pout are all within a bounding box.
 *The units of the in/out vectors are preserved. Numerical problems
 *can appear by using del=0 due to rounding errors. For instance,
 *when the box is partioned by 3x3x3, a terrible thing can happen.
 */
template<class T, unsigned D>
template<class PA>
void 
CrystalLattice<T,D>::applyBC(const PA& pin, PA& pout, T del) const
{
  int mode = pin.InUnit*2+pout.InUnit;
  switch(mode) {
  case(0):
    for(int i=0; i<pin.size(); i++)
      pout[i] = dot(BConds.wrap(dot(pin[i],G)+del),R);
    break;
  case(1):
    for(int i=0; i<pin.size(); i++)
      pout[i] = BConds.wrap(dot(pin[i],G)+del);
    break;
  case(2):
    for(int i=0; i<pin.size(); i++) pout[i] = dot(BConds.wrap(pin[i]+del),R);
    break;
  case(3):
    for(int i=0; i<pin.size(); i++) pout[i] = BConds.wrap(pin[i]+del);
    break;
  }

//   switch(mode) {
//   case(0):
//     for(int i=0; i<pin.size(); i++)
//       pout[i] = dot(BConds.wrap(dot(pin[i],G)),R);
//     break;
//   case(1):
//     for(int i=0; i<pin.size(); i++)
//       pout[i] = BConds.wrap(dot(pin[i],G));
//     break;
//   case(2):
//     for(int i=0; i<pin.size(); i++) pout[i] = dot(BConds.wrap(pin[i]),R);
//     break;
//   case(3):
//     for(int i=0; i<pin.size(); i++) pout[i] = BConds.wrap(pin[i]);
//     break;
//   }

}

/** Copy an in vector to an out vector after applying boundary conditions.
 *
 *@param pos an in/out position array
 *@param del a tolerance for wrapping the values within [0,1)
 *
 *@note The values of pos are all within a bounding box.
 *The unit is preserved. Numerical problems can appear 
 *by using del=0 due to rounding errors.
 */
template<class T, unsigned D>
template<class PA>
void CrystalLattice<T,D>::applyBC(PA& pos, T del) const
{
  if(pos.InUnit) {
    for(int i=0; i<pos.size(); i++) pos[i] = BConds.wrap(pos[i]);
  } else {
    for(int i=0; i<pos.size(); i++) pos[i] = BConds.wrap(dot(pos[i],G)+del);
    pos.InUnit = true;
  }

//   if(pos.InUnit) {
//     for(int i=0; i<pos.size(); i++) pos[i] = BConds.wrap(pos[i]);
//   } else {
//     for(int i=0; i<pos.size(); i++) pos[i] = BConds.wrap(dot(pos[i],G));
//     pos.InUnit = true;
//   }

}

#endif
  


/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
