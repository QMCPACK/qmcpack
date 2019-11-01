//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file CrystalLattice.h
 *@brief Declaration of CrystalLattice<T,D>
 */
#ifndef OHMMS_CRYSTALLATTICE_H
#define OHMMS_CRYSTALLATTICE_H
#include <limits>
#include <iostream>
#include <config/stdlib/math.hpp>
#include <OhmmsPETE/TinyVector.h>
#include <OhmmsPETE/Tensor.h>
#include <Lattice/LRBreakupParameters.h>

namespace qmcplusplus
{
/** enumeration to classify a CrystalLattice
 *
 * Use std::bitset<3> for all the dimension
 */
enum
{
  SUPERCELL_OPEN = 0, /*!< nnn */
  SUPERCELL_WIRE = 1, /*!< nnp */
  SUPERCELL_SLAB = 3, /*!< npp */
  SUPERCELL_BULK = 7, /*!< ppp */
  SOA_OFFSET     = 32 /*!< const to differentiate AoS and SoA */
};

/** enum class to assist copy and unit conversion operations on position vectors
*/
enum class PosUnit
{
  Cartesian = 0, /*!< indicates that the values are in Cartesian units*/
  Lattice        /*!< indicates that the values are in Lattice units*/
};
/** write unit type in human readable format
 *
 *  This could break tools if they rely on parsing log.
 */
std::ostream& operator<<(std::ostream& o_stream, PosUnit pos_unit);
/** Read unit type recorded in int
 *
 *  This should really be human readable
 *  TODO: support both until going to string only.
 */
std::istream& operator>>(std::istream& i_stream, PosUnit& pos_unit);

/** a class that defines a supercell in D-dimensional Euclean space.
 *
 *CrystalLattice handles the physical properties of a supercell, such
 *as lattice vectors, reciprocal vectors and metric tensors and provides
 *interfaces to access the lattice properties and convert units of
 *position vectors or a single-particle position from Cartesian to
 *Lattice Unit vice versa.
 */
template<class T, unsigned D>
struct CrystalLattice : public LRBreakupParameters<T, D>
{
  /// alias to the base class
  using Base = LRBreakupParameters<T, D>;

  ///enumeration for the dimension of the lattice
  enum
  {
    DIM = D
  };
  //@{
  ///the type of scalar
  typedef T Scalar_t;
  ///the type of a D-dimensional position vector
  typedef TinyVector<T, D> SingleParticlePos_t;
  ///the type of a D-dimensional index vector
  typedef TinyVector<int, D> SingleParticleIndex_t;
  ///the type of a D-dimensional Tensor
  typedef Tensor<T, D> Tensor_t;
  //@}

  ///true, if off-diagonal elements are zero so that other classes can take advantage of this
  bool DiagonalOnly;
  ///supercell enumeration
  int SuperCellEnum;
  ///The boundary condition in each direction.
  TinyVector<int, D> BoxBConds;
  ///The scale factor for adding vacuum.
  T VacuumScale;
  //@{
  /**@brief Physical properties of a supercell*/
  /// Volume of a supercell
  Scalar_t Volume;
  /// Wigner-Seitz cell radius
  Scalar_t WignerSeitzRadius;
  /// simulation cell radii
  Scalar_t SimulationCellRadius;
  /// SimulationCellRadius*SimulationCellRadius
  Scalar_t CellRadiusSq;
  /// Wigner-Seitz cell radius in reciprocal space
  Scalar_t WignerSeitzRadius_G;
  ///Real-space unit vectors. R(i,j) i=vector and j=x,y,z
  Tensor_t R;
  ///Reciprocal unit vectors. G(j,i) i=vector and j=x,y,z
  Tensor_t G;
  ///Transpose of reciprocal unit vectors:
  Tensor_t Gt;
  ///Metric tensor
  Tensor_t M;
  ///Metric tensor for G vectors
  Tensor_t Mg;
  ///Length[idim] length of the idim-th lattice vector
  SingleParticlePos_t Length;
  ///OneOverLength[idim] 1/length of the idim-th lattice vector
  SingleParticlePos_t OneOverLength;
  ///Center of the cell sum(Rv[i])/2
  SingleParticlePos_t Center;
  /**@brief Real-space unit vectors.
   *
   *Introduced to efficiently return one vector at a time.
   *Rv[i] is D-dim vector of the ith direction.
   */
  TinyVector<SingleParticlePos_t, D> Rv;
  /**@brief Reciprocal unit vectors.
   *
   *Introduced to efficiently return one vector at a time.
   *Gv[i] is D-dim vector of the ith direction.
   */
  TinyVector<SingleParticlePos_t, D> Gv;
  //@}
  //angles between the two lattice vectors
  SingleParticlePos_t ABC;
  ///true, the lattice is defined by the input instead of an artificial default
  bool explicitly_defined;

  ///default constructor, assign a huge supercell
  CrystalLattice();

  ///copy constructor
  CrystalLattice(const CrystalLattice&) = default;

  ///destructor
  virtual ~CrystalLattice() {}

  /**@param i the index of the directional vector, \f$i\in [0,D)\f$
   *@return The lattice vector of the ith direction
   *@brief Provide interfaces familiar to fotran users
   */
  inline SingleParticlePos_t a(int i) const { return Rv[i]; }

  /**@param i the index of the directional vector, \f$i\in [0,D)\f$
   *@return The reciprocal vector of the ith direction
   *@brief Provide interfaces familiar to fotran users
   */
  inline SingleParticlePos_t b(int i) const { return Gv[i]; }

  /** Convert a cartesian vector to a unit vector.
   * Boundary conditions are not applied.
   */
  template<class T1>
  inline SingleParticlePos_t toUnit(const TinyVector<T1, D>& r) const
  {
    return dot(r, G);
  }

  template<class T1>
  inline SingleParticlePos_t toUnit_floor(const TinyVector<T1, D>& r) const
  {
    SingleParticlePos_t val_dot;
    val_dot = toUnit(r);
    for (int i = 0; i < D; i++)
      if (-std::numeric_limits<T1>::epsilon() < val_dot[i] && val_dot[i] < 0)
        val_dot[i] = T1(0.0);
      else
        val_dot[i] -= std::floor(val_dot[i]);
    return val_dot;
  }

  /** Convert a unit vector to a cartesian vector.
   * Boundary conditions are not applied.
   */
  template<class T1>
  inline SingleParticlePos_t toCart(const TinyVector<T1, D>& c) const
  {
    return dot(c, R);
  }

  /// return true if all the open direction of reduced coordinates u are in the range [0,1)
  inline bool isValid(const TinyVector<T, D>& u) const
  {
    bool inside = true;
    for (int dim = 0; dim < D; dim++)
      inside &= (BoxBConds[dim] || (u[dim] >= T(0) && u[dim] < T(1)));
    return inside;
  }

  /// return true if any direction of reduced coordinates u goes larger than 0.5
  inline bool outOfBound(const TinyVector<T, D>& u) const
  {
    for (int i = 0; i < D; ++i)
      if (std::abs(u[i]) > 0.5)
        return true;
    return false;
  }

  inline void applyMinimumImage(TinyVector<T, D>& c) const
  {
    if (SuperCellEnum)
    {
      TinyVector<T, D> u = dot(c, G);
      for (int i = 0; i < D; ++i)
        u[i] = u[i] - round(u[i]);
      c = dot(u, R);
    }
  }

  /** evaluate the cartesian distance
   *@param ra a vector in the supercell unit
   *@param rb a vector in the supercell unit
   *@return Cartesian distance with two vectors in SC unit
   *
   @note The distance between two cartesian vectors are handled
   *by dot function defined in OhmmsPETE/TinyVector.h
   */
  inline T Dot(const SingleParticlePos_t& ra, const SingleParticlePos_t& rb) const { return dot(ra, dot(M, rb)); }

  /** conversion of a reciprocal-vector
   *@param kin an input reciprocal vector in the Reciprocal-vector unit
   *@return k(reciprocal vector) in cartesian unit
  */
  inline SingleParticlePos_t k_cart(const SingleParticlePos_t& kin) const { return TWOPI * dot(G, kin); }

  /** conversion of a caresian reciprocal-vector to unit k-vector
   *@param kin an input reciprocal vector in cartesian form
   *@return k(reciprocal vector) as unit vector
  */
  inline SingleParticlePos_t k_unit(const SingleParticlePos_t& kin) const { return dot(R, kin) / TWOPI; }

  /** evaluate \f$k^2\f$
   *
   *@param kin an input reciprocal vector in reciprocal-vector unit
   *@return \f$k_{in}^2\f$
   */
  inline T ksq(const SingleParticlePos_t& kin) const { return dot(kin, dot(Mg, kin)); }

  ///assignment operator
  template<typename T1>
  CrystalLattice<T, D>& operator=(const CrystalLattice<T1, D>& rhs)
  {
    Base::LR_dim_cutoff = rhs.LR_dim_cutoff;
    Base::LR_kc         = rhs.LR_kc;
    Base::LR_rc         = rhs.LR_rc;

    explicitly_defined = rhs.explicitly_defined;
    BoxBConds          = rhs.BoxBConds;
    VacuumScale        = rhs.VacuumScale;
    R                  = rhs.R;
    reset();
    return *this;
  }

  /** scale the lattice vectors by sc. All the internal data are reset.
   *@param sc the scaling value
   *@return a new CrystalLattice
   */
  CrystalLattice<T, D>& operator*=(T sc);

  /** set the lattice vector from the command-line options
   *@param lat a tensor representing a supercell
   */
  template<class TT>
  void set(const Tensor<TT, D>& lat);

  /** Evaluate the reciprocal vectors, volume and metric tensor
   */
  void reset();

  //  //@{
  //  /* Copy functions with unit conversion*/
  //  template<class PA> void convert(const PA& pin, PA& pout) const;
  //  template<class PA> void convert2Unit(const PA& pin, PA& pout) const;
  //  template<class PA> void convert2Cart(const PA& pin, PA& pout) const;
  //  template<class PA> void convert2Unit(PA& pout) const;
  //  template<class PA> void convert2Cart(PA& pout) const;
  //  //@}
  //
  //  template<class PA> void applyBC(const PA& pin, PA& pout) const;
  //  template<class PA> void applyBC(PA& pos) const;
  //  template<class PA> void applyBC(const PA& pin, PA& pout, int first, int last) const;

  //! Print out CrystalLattice Data
  void print(std::ostream&, int level = 2) const;
};

} // namespace qmcplusplus
//including the definitions of the member functions
#include "Lattice/CrystalLattice.cpp"

#endif
