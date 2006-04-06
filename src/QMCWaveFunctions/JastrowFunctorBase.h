//////////////////////////////////////////////////////////////////
// (c) Copyright 2005-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
// -*- C++ -*-
#ifndef QMCPLUSPLUS_JASTROWFUNCTORBASE_H
#define QMCPLUSPLUS_JASTROWFUNCTORBASE_H
#include "OhmmsData/libxmldefs.h"
#include "Optimize/VarList.h"
#include "QMCWaveFunctions/OrbitalTraits.h"

namespace qmcplusplus {
/** Base class for any functor used as a source for NumericalJastrow
 */
template<class RT>
struct JastrowFunctorBase {

  typedef QMCDataTrait<RT> data_type;
  ///typedef for the argument
  typedef typename data_type::real_type real_type;
  ///typedef for the return value
  typedef typename data_type::value_type value_type;

  ///reset the Jastrow Function
  virtual void reset()=0;

  /** evaluate the value at r
   * @param r distance
   *
   * virtual function necessary for a transformation to a numerical functor
   */
  virtual value_type f(real_type r)=0;

  /** evaluate the first derivate 
   * @param r distance
   *
   * virtual function necessary for a transformation to a numerical functor
   */
  virtual value_type df(real_type r)=0;

  /** process xmlnode and registers variables to optimize
   * @param cur xmlNode for a functor
   * @param vlist optimizable variables 
  */
  virtual void put(xmlNodePtr cur, VarRegistry<real_type>& vlist) = 0;
};

}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

