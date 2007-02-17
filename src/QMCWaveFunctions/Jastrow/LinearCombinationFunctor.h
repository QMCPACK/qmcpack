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
/** @file LinearConbinationFunctor.h
 * @brief Define a generic functor which consists of multiple OptimizableFunctorBase<T>
 */
#ifndef QMCPLUSPLUS_COMBOFUNCTORWITHCONSTRAINTS_H
#define QMCPLUSPLUS_COMBOFUNCTORWITHCONSTRAINTS_H
#include "Numerics/OptimizableFunctorBase.h"

/** Implements a linear combination of any functors
 *
 * \f$ \phi(r) = \sum_i C_i f_i(r) \f$ where  \f$ f_i(r)\f$ is an one-dimensional
 * functor. Each one-dmensional functor is represented by a derived class from
 * OptimizableFunctorBase<T>. 
 */
template<class T>
struct LinearCombinationFunctor: public OptimizableFunctorBase<T> 
{

  typedef OptimizableFunctorBase<T> ComponentType;
  typedef typename NumericTraits<T>::real_type real_type;
  typedef typename NumericTraits<T>::value_type value_type;
  typedef typename ComponentType::OptimizableSetType OptimizableSetType;

  ///number of ComponentType*
  int NumComponents;
  ///list of linear coefficients
  std::vector<real_type> C;
  ///list of component functors
  std::vector<ComponentType*> Phi;
  ///list of C names
  std::vector<std::string> ID;

  LinearCombinationFunctor(): NumComponents(0)
  { 
    C.reserve(8);
    Phi.reserve(8);
    ID.reserve(8);
  }

  int size() const { return NumComponents;}

  void addComponent(ComponentType* func, real_type c,  const std::string& id, xmlNodePtr cur) 
  {
    C.push_back(c);
    Phi.push_back(func);
    ID.push_back(id);
    NumComponents++;
  }

  inline real_type f(real_type r) {
    real_type res=0;
    for(int i=0; i<NumComponents; i++) res += C[i]*Phi[i]->f(r);
    return res;
  }

  inline real_type df(real_type r) {
    real_type res(0);
    for(int i=0; i<NumComponents; i++) res += C[i]*Phi[i]->df(r);
    return res;
  }

  bool put(xmlNodePtr cur) 
  {
    return true;
  }

  void addOptimizables(OptimizableSetType& vlist) 
  {
    for(int i=0; i<NumComponents; i++) 
      vlist[ID[i]]=C[i];
    for(int i=0; i<NumComponents; i++) 
      Phi[i]->addOptimizables(vlist);
  }

  /** reset the coefficients
   * @param optVariables modified variables
   *
   * - update C[i] if optVariables contains the ID[i]
   * - call resetParameters of the component functors
   */
  void resetParameters(OptimizableSetType& optVariables) 
  {
    for(int i=0; i<NumComponents; i++) 
    {
      typename OptimizableSetType::iterator it(optVariables.find(ID[i]));
      if(it != optVariables.end())
      {
        C[i]=(*it).second;
      }
    }
    for(int i=0; i<NumComponents; i++) 
      Phi[i]->resetParameters(optVariables);
  }

};

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1697 $   $Date: 2007-02-04 15:12:32 -0600 (Sun, 04 Feb 2007) $
 * $Id: OptimizableFunctorBase.h 1697 2007-02-04 21:12:32Z jnkim $ 
 ***************************************************************************/

