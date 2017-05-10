//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef OHMMS_YLMRNLSET_ONGRID_H
#define OHMMS_YLMRNLSET_ONGRID_H
#include <map>
#include <set>
#include <vector>
#include <stdlib.h>
#include "OhmmsData/OhmmsElementBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/OneDimIntegration.h"
#include "OhmmsPETE/TinyVector.h"

using namespace qmcplusplus;
/**@defgroup QuantumNumber
   @brief quantum numbers n, l, m and s
   *
   Several binary functions to compare qauntum numbers are used to sort
   them and map between the quantum numbers and the spherical orbitals.
*/

/**@ingroup QuantumNumber
   @brief binary function to compare the indices representing quantum
   numbers \f$nlm\f$
*/
struct ltnlm
{
  bool operator()(const TinyVector<int,3>& a,
                  const TinyVector<int,3>& b) const
  {
    if(a[0] > b[0])
      return false;
    if(a[0] < b[0])
      return true;
    if(a[1]> b[1])
      return false;
    if(a[1]< b[1])
      return true;
    return a[2]<b[2];
  }
};

/**@ingroup QuantumNumber
   @brief binary function to compare the indices representing quantum
   numbers \f$nl\f$
 */
struct ltnl
{
  bool operator()(const TinyVector<int,2>& a,
                  const TinyVector<int,2>& b) const
  {
    if(a[0] > b[0])
      return false;
    if(a[0] < b[0])
      return true;
    return a[1]<b[1];
  }
};

/**@ingroup QuantumNumber
   @brief binary function to compare the indices representing quantum numbers
   \f$nlm\f$ and spin indices
 */
template<unsigned D>
struct equal_nlms
{

  bool operator()(const TinyVector<int,D>& a,
                  const TinyVector<int,D>& b) const
  {
    int i=0;
    while(i<D)
    {
      if(a[i]!= b[i])
        return true;
      i++;
    }
    return false;
  }
};

/**
 @brief Manage a set of spherical orbitals.
 *
 *Only unique spherical orbitals are stored and std::map is used to map
 * a quantum number set \f$(n,l,m,s)\f$ to a spherical orbital.
 * Each radial orbital is stored in OneDimGridBase<T> and they share a
 * single grid.
 */
template<class T>
struct YlmRnlSet
{

  typedef T                                        value_type;
  typedef OneDimGridBase<value_type>               RadialGrid_t;
  typedef OneDimCubicSpline<value_type,value_type> RadialOrbital_t;
  typedef std::vector<RadialOrbital_t>                  RadialOrbitalSet_t;

  typedef TinyVector<int,2> NLIndex;
  typedef TinyVector<int,3> NLMIndex;
  typedef TinyVector<int,4> NLMSIndex;

  typedef std::map<NLIndex,int,ltnl >  NL_Map_t;
  typedef std::map<NLMIndex,int,ltnlm > NLM_Map_t;
  typedef std::map<NLMSIndex,int,equal_nlms<4> >  NLMS_Map_t;

  ///constructor
  YlmRnlSet(): m_grid(NULL), NumUniqueOrb(0), Restriction("none"),
    CuspParam(0.0), MinEigenValue(0.0), MaxEigenValue(0.0),
    Nup(0), Ndown(0) { }

  ~YlmRnlSet() { }


  bool add(int n, int l, int m, int s, value_type occ);

  //void initialize(const RadialOrbital_t& orb)
  //{
  //  m_grid=orb.m_grid;
  //}

  void applyRestriction(int norb);

  ///normailize all the orbitals
  void normalize(int norb)
  {
    for(int i=0; i < norb; i++)
      normalize_RK2(*psi[i]);
  }

  ///assigns orbital \f$ \psi_{iorb}(r) \f$
  inline RadialOrbital_t& operator()(int iorb)
  {
    return *psi[iorb];
  }

  ///returns a reference to orbital \f$ \psi_{iorb}(r) \f$
  inline const RadialOrbital_t& operator()(int iorb) const
  {
    return *psi[iorb];
  }

  ///return the number of orbitals
  inline int size() const
  {
    return psi.size();
  }

  ///return \f$ \psi_{iorb}(r_j) \f$
  inline value_type operator()(int iorb, int j) const
  {
    return psi[iorb]->operator()(j);
  }

  ///assigns orbital \f$ \psi_{iorb}(r_j) \f$
  inline value_type& operator()(int iorb, int j)
  {
    return psi[iorb]->operator()(j);
  }

  ///returns the derivative at the first grid point
  inline value_type first_deriv(int iorb) const
  {
    return psi[iorb]->yprime0;
  }

  ///returns the derivative at the last grid point
  inline value_type last_deriv(int iorb) const
  {
    return psi[iorb]->yprimeN;
  }

  ///reset the values of orbitals
  void reset()
  {
    for(int i=0; i < psi.size(); i++)
      psi[i]->m_Y = 0.0;
  }

  ///restriction type;
  std::string Restriction;

  ///number of unique orbitals
  int NumUniqueOrb;

  ///a common grid
  RadialGrid_t* m_grid;

  ///principal quantum number for each orbital
  std::vector<int> N;
  ///angular momentum quantum number for each orbital
  std::vector<int> L;
  ///z-angular momentum quantum number for each orbital
  std::vector<int> M;
  ///spin for each orbital (up = 1, down = -1)
  std::vector<int> S;
  ///coefficient for each orbtial (1 if occupied, 0 if not)
  std::vector<value_type> Occ;

  ///number of spin-up orbitals
  int Nup;
  ///number of spin-down orbitals
  int Ndown;
  ///charge of the nuclei
  int Charge;
  ///parameter for calculating orbital cusp conditions
  value_type CuspParam;

  ///lower-bound for the eigenvalues
  value_type MinEigenValue;

  ///upper-bound for the eigenvalues
  value_type MaxEigenValue;

  ///location of the radial function with \f$(n,l)\f$ orbital(s)
  NL_Map_t NL;

  ///location of the radial function with \f$(n,l,m)\f$ orbital(s)
  NLM_Map_t NLM;

  ///map for occupation number
  NLMS_Map_t OccNo;

  ///assign an id for each orbital
  std::vector<int> ID;

  ///keeps track of the number of orbitals with the same id
  std::vector<int> IDcount;

  ///the radial grid orbitals
  std::vector<RadialOrbital_t*> psi;

  bool put(xmlNodePtr cur);
  bool get(std::ostream& os);

  bool print_HDF5(const std::string&, const std::string&,
                  const std::vector<value_type>&);

  bool print_basis(const std::string&, const std::string&,
                   const std::string&);

};
#include "SQD/YlmRnlSet.cpp"
#include "SQD/YlmRnlSet.IO.h"
#endif

