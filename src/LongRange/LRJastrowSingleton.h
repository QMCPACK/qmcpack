//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnum Kim
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
/** @file LRJastrowSingleton.h
 * @brief Define a LRHandler with two template parameters
 */
#ifndef QMCPLUSPLUS_LRJASTROWSINGLETON_H
#define QMCPLUSPLUS_LRJASTROWSINGLETON_H

#include "LongRange/LRHandlerTemp.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"

namespace qmcplusplus {
  /** JastrowFunctor
   * 
   * A Func for LRHandlerTemp.  Four member functions have to be provided 
   * 
   * - reset(T volume) : reset the normalization factor
   * - operator() (T r, T rinv) : return a value of the original function e.g., 1.0/r
   * - Fk(T k, T rc)
   * - Xk(T k, T rc)
   * 
   */     

  template<class T=double>
  struct JastrowFunctor {
    T OneOverSqrtRs;
    T Rs;
    T NormFactor;
    inline JastrowFunctor(){}
    
    void reset(ParticleSet& ref) {
      NormFactor=4.0*M_PI/ref.Lattice.Volume;
      T Density=ref.getTotalNum()/ref.Lattice.Volume;
      Rs = std::pow(3.0/(4.0*M_PI*Density), 1.0/3.0);
      OneOverSqrtRs = 1.0 / std::sqrt(Rs);
    }
     
    void reset(ParticleSet& ref, T rs) {
       NormFactor=4.0*M_PI/ref.Lattice.Volume;
       Rs = rs;
       OneOverSqrtRs = 1.0 / std::sqrt(Rs);
    }
     
    inline T operator()(T r, T rinv) { 
       return Rs*rinv*(1.0 - std::exp(-r*OneOverSqrtRs));
    }
    
    inline T Fk(T k, T rc) {
       return -Xk(k,rc);
    }
    
    inline T Xk(T k, T rc) {
       T coskr = std::cos(k*rc);
       T sinkr = std::sin(k*rc);
       T oneOverK = 1.0/k;
       return -NormFactor * Rs * 
	  ( coskr*oneOverK*oneOverK - std::exp(-rc*OneOverSqrtRs) * 
	        (coskr - OneOverSqrtRs * sinkr * oneOverK) / (k*k + 1 / Rs) );
    }
  };


  struct LRJastrowSingleton {

    typedef LRHandlerTemp<JastrowFunctor<double>,LPQHIBasis> LRHandlerType;
    
    static LRHandlerType* JastrowHandler;

    static LRHandlerType* getHandler(ParticleSet& ref);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
