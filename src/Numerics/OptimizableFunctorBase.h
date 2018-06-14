//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/** @file OptimizableFunctorBase.h
 * @brief Define a base class for one-dimensional functions with optimizable variables
 */
#ifndef QMCPLUSPLUS_OPTIMIZABLEFUNCTORBASE_H
#define QMCPLUSPLUS_OPTIMIZABLEFUNCTORBASE_H
#include "Optimize/VariableSet.h"
#include "OhmmsData/OhmmsElementBase.h"
#include "OhmmsPETE/TinyVector.h"

/** Base class for any functor with optimizable parameters
 *
 * Derived classes from OptimizableFunctorBase are called "functor"s and
 * can be used as a template signature for  Jastrow functions.
 * - OneBodyJastroOrbital<FUNC>
 * - TwoBodyJastroOrbital<FUNC>
 * Functor in qmcpack denotes any function which returns a value at a point, e.g.,
 * GTO, STO, one-dimensional splines etc. OptimizableFunctorBase is introduced for
 * optimizations. The virtual functions are intended for non-critical operations that
 * are executed infrequently during optimizations.
 *
 * This class handles myVars of opt_variables_type (Optimize/VariableSet.h). A derived class
 * can insert any number of variables it handles during optimizations, by calling
 * myVars.insert(name,value);
 * Unlike VarList which uses map, VariableSet is serialized in that the internal order is according
 * to insert calls.
 */
struct OptimizableFunctorBase
{
  ///typedef for real values
  typedef optimize::VariableSet::real_type real_type;
  ///typedef for variableset: this is going to be replaced
  typedef optimize::VariableSet opt_variables_type;
  ///typedef for name-value lists
  typedef optimize::VariableSet::variable_map_type variable_map_type;
  ///maximum cutoff
  real_type cutoff_radius;
  ///set of variables to be optimized
  opt_variables_type myVars;
  ///default constructor
  inline OptimizableFunctorBase() {}
  ///virtual destrutor
  virtual ~OptimizableFunctorBase() {}

  inline void getIndex(const opt_variables_type& active)
  {
    myVars.getIndex(active);
  }

  virtual void checkInVariables(opt_variables_type& active)=0;

  virtual void checkOutVariables(const opt_variables_type& active)=0;

  /** reset the optimizable variables
   * @param active list of active optimizable variables
   */
  virtual void resetParameters(const opt_variables_type& active)=0;
  /** create a clone of this object
   */
  virtual OptimizableFunctorBase* makeClone() const =0;

  /** reset function
   */
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

  /** empty virtual function to help builder classes
  */
  virtual void setDensity(real_type n) { }

  virtual inline bool evaluateDerivatives (real_type r, std::vector<qmcplusplus::TinyVector<real_type,3> >& derivs)
  {
    return false;
  }

// mmorales: don't know how to solve a template problem for cusp correction,
//           so for now I do this
  virtual void setGridManager(bool willmanage) { }

};

#endif

