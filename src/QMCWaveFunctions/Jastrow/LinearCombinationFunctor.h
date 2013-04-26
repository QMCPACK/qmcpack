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
struct LinearCombinationFunctor: public OptimizableFunctorBase
{

  typedef OptimizableFunctorBase ComponentType;

  ///number of ComponentType*
  int NumComponents;
  ///list of bool
  std::vector<bool> CanNotChange;
  ///list of linear coefficients
  std::vector<real_type> C;
  ///list of component functors
  std::vector<ComponentType*> Phi;

  LinearCombinationFunctor(): NumComponents(0)
  {
    CanNotChange.reserve(8);
    C.reserve(8);
    Phi.reserve(8);
  }

  OptimizableFunctorBase* makeClone() const
  {
    LinearCombinationFunctor<T>* myclone=new LinearCombinationFunctor<T>(*this);
    for(int i=0; i<NumComponents; ++i)
      myclone->Phi[i]=Phi[i]->makeClone();
    return myclone;
  }

  int size() const
  {
    return NumComponents;
  }

  void addComponent(ComponentType* func, real_type c,  std::string& id, bool fixit=false)
  {
    CanNotChange.push_back(fixit);
    C.push_back(c);
    Phi.push_back(func);
    int loc=myVars.size();
    myVars.insert(id,c);
    if(fixit)
      myVars.Index[loc]=-1;//freeze this
    NumComponents++;
  }

  inline void reset()
  {
    for(int i=0; i<NumComponents; ++i)
      Phi[i]->reset();
  }


  inline real_type f(real_type r)
  {
    real_type res=0;
    for(int i=0; i<NumComponents; i++)
      res += C[i]*Phi[i]->f(r);
    return res;
  }

  inline real_type df(real_type r)
  {
    real_type res(0);
    for(int i=0; i<NumComponents; i++)
      res += C[i]*Phi[i]->df(r);
    return res;
  }

  bool put(xmlNodePtr cur)
  {
    return true;
  }


  //disable optimization of exponents
  void checkInVariables(opt_variables_type& active)
  {
    active.insertFrom(myVars);
    //for(int i=0; i<NumComponents; i++) Phi[i]->checkInVariables(active);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active);
    //for(int i=0; i<NumComponents; i++) Phi[i]->checkOutVariables(active);
  }


  /** reset the coefficients
   * @param optVariables modified variables
   *
   * - update C[i] if optVariables contains the ID[i]
   * - call resetParameters of the component functors
   */
  void resetParameters(const opt_variables_type& active)
  {
    for(int i=0; i<NumComponents; ++i)
    {
      if(CanNotChange[i])
        continue;
      int loc=myVars.where(i);
      if(loc>=0)
        C[i]=myVars[i]=active[loc];
    }
    //for(int i=0; i<NumComponents; i++) Phi[i]->resetParameters(active);
  }

};

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1697 $   $Date: 2007-02-04 15:12:32 -0600 (Sun, 04 Feb 2007) $
 * $Id: OptimizableFunctorBase.h 1697 2007-02-04 21:12:32Z jnkim $
 ***************************************************************************/

