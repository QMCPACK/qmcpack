//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_OPTIMIZATIONFUNCION_BASE_H
#define QMCPLUSPLUS_OPTIMIZATIONFUNCION_BASE_H

#include "OhmmsData/OhmmsElementBase.h"
#include "Optimize/LeastSquaredFit.h"

/** Base class for any cost function
 *
 * Based on K. Esler's MinimizeFunction.
 */
template<class T=double>
class CostFunctionBase
{

public:

  typedef T Return_t;

  /** boolean to indicate if the cost function is valid.
   *
   * Can be used by optimizers to stop optimization.
   */
  bool IsValid;

  CostFunctionBase():IsValid(true) { }

  virtual ~CostFunctionBase() {}

  virtual int NumParams() = 0;

  virtual Return_t& Params(int i) = 0;

  virtual Return_t Params(int i) const = 0;

  virtual Return_t Cost(bool needGrad=true) = 0;

  virtual void GradCost(std::vector<Return_t>& PGradient, const std::vector<Return_t>& PM, Return_t FiniteDiff=0) = 0;

  virtual void Report() = 0;

  /** find dl that minimize an object function in cg direction
   * @param x0 current parameters
   * @param cg direction
   * @param dl displacement
   * @param val optimal value
   * @param lambda_max maximum displacement
   * @return true, if lineoptimization is done by the derived class
   *
   * Default implementation returns false to perform a simple method
   */
  virtual bool lineoptimization(const std::vector<T>& x0, const std::vector<T>& gr,
                                Return_t val0,
                                Return_t& dl, Return_t& vopt, Return_t& lambda_max)
  {
    return false;
  }

  /** evaluate gradient, \f$ gr(\alpha_i)= -\frac{\partial F}{\partial \alhpa_i}\f$
   * @param gradients
   * @return true, if a cost function evaluates the gradients
   *
   * Default implementation returns false to perform a finite-difference method
   * for gradients.
   */
  virtual bool evaluateGradients(std::vector<T>& gr)
  {
    return false;
  }
};

/** base class for optimization(minimization) classes
 *
 * Template parameter T is the numerical type
 */
template<class T=double>
struct MinimizerBase
{

  /** stream to write intermediate message
   */
  std::ostream* msg_stream;

  /** typedef of the object function to be optimized
   */
  typedef CostFunctionBase<T> ObjectFuncType;

  /** default constructor */
  MinimizerBase():msg_stream(0) {}

  /** virtual destructor */
  virtual ~MinimizerBase() {}

  /** set msg_stream
   * @param os_ptr pointer to std::ostream
   */
  void setOstream(std::ostream* os_ptr)
  {
    msg_stream = os_ptr;
  }

  /** optimize an object function
   * @param atarget ObjectFuncType to be optimized
   * @return true if converged
   */
  virtual bool optimize(ObjectFuncType* atarget) = 0;

  virtual bool put(xmlNodePtr cur)=0;
};

#endif
