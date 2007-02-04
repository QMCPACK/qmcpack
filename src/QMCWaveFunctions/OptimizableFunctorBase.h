//////////////////////////////////////////////////////////////////
// (c) Copyright 2005-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
/** @file OptimizableFunctorBase.h
 * @brief Define a base class for one-dimensional functions that can be optimized.
 */
#ifndef QMCPLUSPLUS_OPTIMIZABLEFUNCTORBASE_H
#define QMCPLUSPLUS_OPTIMIZABLEFUNCTORBASE_H
#include <complex>
#include "OhmmsData/libxmldefs.h"
#include "Optimize/VarList.h"

template<class T> struct NumericTraits {};

template<>
struct NumericTraits<double> {
  typedef double          real_type;
  typedef double          value_type;
  typedef std::complex<double> complex_type;
};

template<>
struct NumericTraits<std::complex<double> > {
  typedef double          real_type;
  typedef std::complex<double> value_type;
  typedef std::complex<double> complex_type;
};

/** Base class for any functor used as a source for NumericalJastrow
*/
template<class T>
struct OptimizableFunctorBase: public NumericTraits<T> {

  typedef typename NumericTraits<T>::real_type real_type;
  typedef typename NumericTraits<T>::value_type value_type;

  ///reset the Jastrow Function
  virtual void reset()=0;

  /** evaluate the value at r
   * @param r distance
   *
   * virtual function necessary for a transformation to a numerical functor
   */
  virtual real_type f(real_type r)=0;

  /** evaluate the first derivate 
   * @param r distance
   *
   * virtual function necessary for a transformation to a numerical functor
   */
  virtual real_type df(real_type r)=0;

  /** process xmlnode and registers variables to optimize
   * @param cur xmlNode for a functor
   */
  virtual bool put(xmlNodePtr cur) = 0;

  /** add variables to be optimized
   * @param vlist VarRegistery<T>
   */
  virtual void addOptimizables(VarRegistry<real_type>& vlist) =0;

  /** empty virtual function to help builder classes
  */
  virtual void setDensity(real_type n) { }
};

/** Implements a linear combination of any functor
*/
template<class T>
struct ComboFunctor: public OptimizableFunctorBase<T> {

  typedef OptimizableFunctorBase<T> ComponentType;
  typedef typename NumericTraits<T>::real_type real_type;
  typedef typename NumericTraits<T>::value_type value_type;

  std::vector<real_type> C;
  std::vector<ComponentType*> Phi;
  std::vector<std::string> ID;

  ComboFunctor() { 
    C.reserve(8);
    Phi.reserve(8);
    ID.reserve(8);
  }

  int size() const { return Phi.size();}

  void add(ComponentType* func, real_type c,  const std::string& id) {
    C.push_back(c);
    Phi.push_back(func);
    ID.push_back(id);
  }

  inline void reset() {
    for(int i=0; i<Phi.size(); i++) Phi[i]->reset();
  }

  inline real_type f(real_type r) {
    real_type res=0;
    for(int i=0; i<Phi.size(); i++) { res += C[i]*Phi[i]->f(r);}
    return res;
  }

  inline real_type df(real_type r) {
    real_type res(0);
    for(int i=0; i<Phi.size(); i++) { res += C[i]*Phi[i]->df(r);}
    return res;
  }

  bool put(xmlNodePtr cur) 
  {
    return true;
  }

  void addOptimizables(VarRegistry<real_type>& vlist) {
    for(int i=0; i<C.size(); i++) {
      vlist.add(ID[i],&(C[i]),1);
    }
    for(int i=0; i<Phi.size(); i++) {
      Phi[i]->addOptimizables(vlist);
    }
  }
};


#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

