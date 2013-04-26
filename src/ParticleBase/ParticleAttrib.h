//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//
//   Jeongnim Kim
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
/**@file ParticleAttrib.h
 *@brief Declaraton of ParticleAttrib<T>
 */

/*! \class ParticleAttrib
 *  \brief A one-dimensional vector class based on PETE.
 *
 *  Closely related to PETE STL vector example and written to substitute
 *  Poomar1::ParticleInteractAttrib.
 *
 *  Equivalent to blitz::Array<T,1>, pooma::Array<1,T>.
 *
 *  class C is a container class. Default is vector<T>
 *
 *  \todo Implement openMP compatible container class or evaluate function.
 *  \todo Implement get/put member functions for MPI-like parallelism
 */
#ifndef OHMMS_PARTICLEATTRIB_PEPE_H
#define OHMMS_PARTICLEATTRIB_PEPE_H

#include "PETE/PETE.h"
#include <vector>
#include <string>
#include <iostream>
#include "Utilities/OhmmsObject.h"

//#ifndef DEBUGMSG
//#define DEBUGMSG(msg)
//#endif
namespace qmcplusplus
{
template<class T>
class ParticleAttrib: public OhmmsObject
{
public:

  typedef T           Type_t;
  typedef vector<T>   Container_t;//! Using stl::vector as a container
  typedef typename Container_t::iterator iterator;
  typedef typename Container_t::const_iterator const_iterator;
  typedef ParticleAttrib<T> This_t;

  /// The unit type
  int InUnit;
  /// The number of local particle attributes
  int nLocal;
  /// The number of ghost particle attributes
  int nGhosts;

  inline ParticleAttrib(): InUnit(0),nLocal(0),nGhosts(0) {}

  /** @defgroup PtclAttribConst Constructors of ParticleAttrib
   * The ghost elements are NEVER copied.
   * @{
   */
  /**@brief nothing is created but the type and object name */
  inline ParticleAttrib(const std::string& tname, const std::string& oname):
    OhmmsObject(tname, oname), InUnit(0), nLocal(0), nGhosts(0)
  {
    ElementByteSize = sizeof(T);
  }

  /**@brief n-element is created but the type and object names */
  inline ParticleAttrib(const std::string& tname, const std::string& oname, int n):
    OhmmsObject(tname,oname), InUnit(0), nLocal(0), nGhosts(0)
  {
    ElementByteSize = sizeof(T);
    resize(n);
    assign(*this, T());
  }

  /**@brief n-element without type and object names **/
  explicit inline ParticleAttrib(unsigned n):InUnit(0), nLocal(0), nGhosts(0)
  {
    resize(n);
    for(int i=0; i<n; i++)
      X[i]=T();
    //assign(*this, T());
  }

  //! Destructor
  inline ~ParticleAttrib() { }

  //! return the current size
  inline unsigned size() const
  {
    return nLocal; //return X.size()-nGhosts;
  }

  //! resize the container (probably, should be removed)
  void resize(unsigned n);

  //! add n elements
  void create(unsigned n);

  //@{handling ghost elements
  //! return the number of ghost attributes
  inline int getNumGhosts() const
  {
    return nGhosts;
  }

  inline void clear()
  {
    X.clear();
    nLocal = 0;
    nGhosts = 0;
  }

  //! remove ghost elements
  inline void clearGhosts()
  {
    if(nGhosts)
    {
      X.erase(X.end()-nGhosts, X.end());
      nGhosts=0;
    }
  }

  //! add ghost elements
  inline void addGhosts(int n)
  {
    if(n)
    {
      X.insert(X.end(), n, T());
      nGhosts+=n;
    }
  }
  //@}

  // Assignment Operators
  inline This_t& operator=(const ParticleAttrib<T> &rhs)
  {
    return assign(*this,rhs);
  }

  inline const This_t& operator=(const ParticleAttrib<T> &rhs) const
  {
    return assign(*this, rhs);
  }

  template<class RHS>
  inline This_t& operator=(const RHS& rhs)
  {
    return assign(*this,rhs);
  }

  // Get and Set Operations
  inline Type_t& operator[](unsigned int i)
  {
    return X[i];
  }

  inline Type_t operator[](unsigned int i) const
  {
    return X[i];
  }

  inline Type_t& operator()(unsigned int i)
  {
    return X[i];
  }

  inline Type_t operator()( unsigned int i) const
  {
    return X[i];
  }

  ///write to a ostream
  bool get(std::ostream& ) const;

  ///read from istream
  bool put(std::istream& );

  ///read from an xmlNode
  bool put(xmlNodePtr cur);

  ///reset member data
  void reset() { }

  //@{set/set the unit
  inline void setUnit(int i)
  {
    InUnit = i;
  }
  inline int getUnit() const
  {
    return InUnit;
  }
  //@}

  //@{iterators consistent with stl::iterator
  inline const_iterator begin() const
  {
    return X.begin();
  }
  inline const_iterator end() const
  {
    return X.end()-nGhosts;
  }

  inline iterator begin()
  {
    return X.begin();
  }
  inline iterator end()
  {
    return X.end()-nGhosts;
  }

  inline Type_t* first_address()
  {
    return &(X[0]);
  }
  inline const Type_t* first_address() const
  {
    return &(X[0]);
  }

  inline Type_t* last_address()
  {
    return &(X[0])+nLocal;
  }
  inline const Type_t* last_address() const
  {
    return &(X[0])+nLocal;
  }
  //@}

  inline void begin_node(std::ostream& os) const
  {
    os << "<attrib name=\"" <<  objName()
       << "\" datatype=\""<< typeName()
       << "\" size=\""<< size()
       << "\" condition=\"" << InUnit
       << "\">" << endl;
  }

  inline void end_node(std::ostream& os) const
  {
    os << "</attrib>" << std::endl;
  }


  //----------------------------------------------------------------------
  // parallel communication
  //Message& putMessage(Message& m) const {
  //  m.setCopy(true);
  //  ::putMessage(m, X.begin(), X + D);
  //    return m;
  //}

  //Message& getMessage(Message& m) {
  //  ::getMessage(m, X, X + D);
  //  return m;
  //}


private:
  Container_t X;
};

///////////////////////////////////////////////////////////
// member functions
///////////////////////////////////////////////////////////
/** resize the attribute
 *
 *  Although resize exists, there is no reason to use it.
 * Use, ParticleAttrib<T>::create which add a number of particle attributes.
 * \warning For reasons I don't understand, KCC/SGI-CC on Origin cannot
 * resize stl::vector with X.resize(n) with optimization.
 * KCC on other platforms work perfectly fine.
 */
template<class T>
void ParticleAttrib<T>::resize(unsigned n)
{
  if(n)
  {
    X = Container_t(n);
    nLocal = n;
    nGhosts = 0;
  }
}

/** Adding n elements for all the particle attributes added by addAttrubte
 * @param n number of elements to be appended.
 */
template<class T>
void ParticleAttrib<T>::create(unsigned n)
{
  Container_t y(X);
  X = Container_t(n+y.size());
  for(int i=0; i<y.size(); i++)
    X[i] = y[i];
  nLocal+=n;
  //if(n) X.insert(X.end(), n, T());
}


/** Specialized to write the unit
 *\return true, if the attribute is relative to a unit
 *
 * Ad hoc get function to tell users if this attribute is absolute or relative value.
 * E.g., is Cartesian unit vs Lattice unit.
 */
template<class T>
bool ParticleAttrib<T>::get(std::ostream& os) const
{
  os << InUnit;
  return true;
}

/** Specialized to read the unit
 */
template<class T>
bool ParticleAttrib<T>::put(std::istream& is)
{
  is >> InUnit;
  return true;
}

/*@warning not fully implemented.*/
template<class T>
bool ParticleAttrib<T>::put(xmlNodePtr cur)
{
  return true;
}

}
#include "ParticleBase/ParticleAttrib.cpp"
#endif // OHMMS_PARTICLEATTRIB_PEPE_H


/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
