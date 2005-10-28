//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim 
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_OPTIMIZATIONFUNCION_BASE_H
#define QMCPLUSPLUS_OPTIMIZATIONFUNCION_BASE_H

#include "OhmmsData/OhmmsElementBase.h"

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

  virtual Return_t Cost() = 0;

  virtual void Report() = 0;

};

/** base class for optimization(minimization) classes 
 *
 * Template parameter T is the numerical type
 */
template<class T=double>
struct MinimizerBase {

  /** stream to write intermediate message
   */
  ostream* msg_stream;

  /** typedef of the object function to be optimized
   */
  typedef CostFunctionBase<T> ObjectFuncType;

  /** default constructor */
  MinimizerBase():msg_stream(0) {}

  /** virtual destructor */
  virtual ~MinimizerBase(){}

  /** set msg_stream 
   * @param os_ptr pointer to ostream
   */
  void setOstream(ostream* os_ptr) { 
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
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
