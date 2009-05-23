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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_LINEAR_NRC_H
#define QMCPLUSPLUS_LINEAR_NRC_H
 
#include "Optimize/NRCOptimization.h" 
#include "QMCDrivers/QMCCostFunctionBase.h"
namespace qmcplusplus {
  
  class QMCCostFunctionBase;

  template<class T>
  struct LinearLineMin: public NRCOptimization<T> {

    typedef T Return_t;
    
    ///target cost function to optimize
    QMCCostFunctionBase* optTarget;
    ///Parameters 
    vector<Return_t> Parms;
    ///Direction to minimize
    vector<Return_t> Grad;

    /** constructor
    */
    LinearLineMin() {}

    /** destructor */
    ~LinearLineMin() {}
    
    void setFunc(QMCCostFunctionBase* optPtr, vector<Return_t> ggevdir, vector<Return_t> curP)
    {
      optTarget=optPtr;
      Grad=ggevdir;
      Parms=curP;
    }
    
    Return_t Func(Return_t dl)
    {
      vector<Return_t> Pdl(Parms);
      for (int i=0;i<Parms.size();i++) Pdl[i] += dl*Grad[i];
      for (int i=0;i<Parms.size();i++) optTarget->Params(i) = Pdl[i];
      return optTarget->Cost();
    }
  
  };
}
/***************************************************************************
 * $RCSfile$   $Author: kesler $
 * $Revision: 3531 $   $Date: 2009-02-09 16:09:02 -0600 (Mon, 09 Feb 2009) $
 * $Id: LinearLineMin.h 3531 2009-02-09 22:09:02Z kesler $ 
 ***************************************************************************/
#endif
