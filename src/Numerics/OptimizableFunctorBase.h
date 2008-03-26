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
 * @brief Define a base class for one-dimensional functions with optimizable variables
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

  ///define the real type
  typedef typename NumericTraits<T>::real_type  real_type;
  ///define the value type
  typedef typename NumericTraits<T>::value_type value_type;
  ///define the type of Optimizable Sets
  typedef VarRegistry<real_type>                OptimizableSetType;

  ///index of the first optimizable variable
  int FirstIndex;
  ///index of the last optimizable variable
  int LastIndex;
  ///default constructor
  OptimizableFunctorBase():FirstIndex(0),LastIndex(1) {}
  ///virtual destrutor
  virtual ~OptimizableFunctorBase(){}

  ///return the total number of variables to optimize
  inline int getNumOfVariables() const {return LastIndex-FirstIndex;}

  /** set the index bounds of the variables to optimize
   *
   * This is to utilize the vectorized container of variables.
   */
  inline void setBounds(int first, int last=-1)
  {
    FirstIndex=first;
    LastIndex=(last>first)?last:first+1;
  }

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
   * @param vlist list to which  derived classes add optimizable variables
   */
  virtual void addOptimizables(OptimizableSetType& vlist) =0;

  /** reset the optimizable variables
   *
   * @param optVariables list of active optimizable variables
   */
  virtual void resetParameters(OptimizableSetType& optVariables)=0;

  /** empty virtual function to help builder classes
  */
  virtual void setDensity(real_type n) { }
};

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

