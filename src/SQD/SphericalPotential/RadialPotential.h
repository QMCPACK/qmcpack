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
    
    


#ifndef OHMMSHF_RADIALPOTENTIALBASE_H
#define OHMMSHF_RADIALPOTENTIALBASE_H

#include "SQD/HFConfiguration.h"
#include "OhmmsData/libxmldefs.h"

class Clebsch_Gordan;
/**
 *@defgroup RadialPotential
 *@brief Classes to define potentials on a radial grid.
 *
 *A radial potential is defined on a radial grid,
 *OneDimGridBase<T> and RadialPotentialSet are used to represent
 *the LHS of Radial Schrodinger Equation.
 */
namespace ohmmshf
{

/**
 *@ingroup RadialPotential
 *@class RadialPotentialBase
 *@brief An abstract base class for a Radial Potential.
 *
 *Inherited classes implement the member function evaluate
 *to calculate matrix elements for a Hamiltonian term.
 */
struct RadialPotentialBase
{

  typedef SphericalOrbitalTraits::BasisSetType       BasisSetType;
  typedef SphericalOrbitalTraits::value_type         value_type;
  typedef SphericalOrbitalTraits::RadialGrid_t       RadialGrid_t;
  typedef SphericalOrbitalTraits::RadialOrbital_t    RadialOrbital_t;
  typedef SphericalOrbitalTraits::RadialOrbitalSet_t RadialOrbitalSet_t;

  ///lower-bound for eigenvalues
  value_type MinEigenValue;

  ///upper-bound for eigenvalues
  value_type MaxEigenValue;

  ///the core charge at infinity
  value_type Qinfty;
  ///the cutoff radius
  value_type Rcut;

  ///storage for an external potential
  RadialOrbital_t* Vext;
  ///storage for a temporary integrand
  RadialOrbital_t* integrand;

  ///constructor
  RadialPotentialBase()
    : MinEigenValue(0.0), MaxEigenValue(0.0), Qinfty(0.0), Rcut(1.0)
    , Vext(0), integrand(0)
  {}

  ///destructor
  virtual ~RadialPotentialBase();

  /**
   *@param psi the wavefunction
   *@param V the potential
   *@param norb the number of orbitals
   *@return The sum of the matrix elements of a Radial Potential \f$V(r)\f$
   *@brief Calculates and assigns the values of a Randial Potential for
   *each orbital.
   *
   \f[
   \sum_{k=0}^{N_{orb}} \langle k|V(r)|k \rangle = \sum_{k=0}^{N_{orb}}
   \int_0^{\infty} dr \: \psi_k^*(r)V(r)\psi_k(r)
   \f]
  */
  virtual
  value_type evaluate(const BasisSetType& psi,
                      RadialOrbitalSet_t& V,
                      int norb) = 0;

  /**
   *@param n the principal quantum number
   *@param l the angular quantum number
   *@return the number of radial nodes
   */
  virtual int getNumOfNodes(int n, int l)=0;

  /**
   *@param RootFileName the name of file root
   *@brief Output the internal storage of the potential.
   *@note Only applies for the Hartree and Exchange potentials.
   */
  virtual void getStorage(const BasisSetType& psi,
                          const std::string& RootFileName);

  /**
   *@param ig the grid index
   *@return the value of the external potential
  */
  inline
  value_type getV(int ig) const
  {
    return (*Vext)(ig);
  }


  /**@return the pointer of the data **/
  inline value_type* data() const
  {
    if(Vext)
      return Vext->data();
    else
      return NULL;
  }

  /**
   *@param cur the current xml node to process
   *@return true if the input is valid
   */
  virtual bool put(xmlNodePtr cur)
  {
    return true;
  }

};

/**
 *@ingroup RadialPotential
 *@class HartreePotential
 *@brief Implements the Hartree potential.
*/
struct HartreePotential: public RadialPotentialBase
{
  ///the Clebsch Gordan coefficient matrix
  Clebsch_Gordan *CG_coeff;
  ///aux functors
  RadialOrbital_t *Ykii_r;
  RadialOrbital_t *Ykjj_r;

  /**store the matrix elements
     \f$\langle \psi_i \psi_j |V| \psi_i \psi_j \rangle \f$ */
  Vector<value_type> storage;
  HartreePotential(Clebsch_Gordan* cg, int norb);
  ~HartreePotential();
  value_type evaluate(const BasisSetType& psi,
                      RadialOrbitalSet_t& V, int norb);
  void getStorage(const BasisSetType& psi,
                  const std::string& RootFileName);
  int getNumOfNodes(int n, int l)
  {
    return 0;
  }
};

/**
 *@ingroup RadialPotential
 *@class ExchangePotential
 *@brief Implements the exchange potential
 */
struct ExchangePotential: public RadialPotentialBase
{
  ///the Clebsch Gordan coefficient matrix
  Clebsch_Gordan *CG_coeff;
  RadialOrbital_t *Ykij_r;
  /**store the matrix elements
     \f$\langle \psi_i \psi_j |V| \psi_j \psi_i \rangle \f$ */
  Vector<value_type> storage;
  ExchangePotential(Clebsch_Gordan* cg, int norb);
  ~ExchangePotential();
  value_type evaluate(const BasisSetType& psi,
                      RadialOrbitalSet_t& V, int norb);
  void getStorage(const BasisSetType& psi,
                  const std::string& RootFileName);
  int getNumOfNodes(int n, int l)
  {
    return 0;
  }
};
}
#endif


