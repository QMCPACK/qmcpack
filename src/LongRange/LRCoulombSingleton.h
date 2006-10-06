//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
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
/** @file LRCoulombSingleton.h
 * @brief Define a LRHandler with two template parameters
 */
#ifndef QMCPLUSPLUS_LRCOULOMBSINGLETON_H
#define QMCPLUSPLUS_LRCOULOMBSINGLETON_H

#include "LongRange/LRHandlerTemp.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"

namespace qmcplusplus {

  /** CoulombFunctor
   *
   * An example for a Func for LRHandlerTemp. Four member functions have to be provided
   * - reset(T volume) : reset the normalization factor
   * - operator()(T r, T rinv): return a value of the original function, e.g., 1.0/r
   * - Fk(T k, T rc)
   * - Xk(T k, T rc)
   */
  template<class T=double>
  struct CoulombFunctor {
    T NormFactor;
    inline CoulombFunctor(){}
    void reset(T volume) {
      NormFactor=4.0*M_PI/volume;
    }
    inline T operator()(T r, T rinv) { return rinv;}
    inline T Fk(T k, T rc) {
      return NormFactor/(k*k)* std::cos(k*rc);
    }
    inline T Xk(T k, T rc) {
      return -NormFactor/(k*k)* std::cos(k*rc);
    }
    inline T df(T r) { return -1.0/r/r;}
  };

  template<class T=double>
  struct PseudoCoulombFunctor {
    typedef OneDimCubicSpline<T> RadialFunctorType;
    RadialFunctorType& radFunc;
    T NormFactor;
    inline PseudoCoulombFunctor(RadialFunctorType& rfunc):radFunc(rfunc){}
    void reset(T volume) {
      NormFactor=4.0*M_PI/volume;
    }
    inline T operator()(T r, T rinv) { return radFunc.splint(r);}
    inline T Fk(T k, T rc) {
      return NormFactor/(k*k)* std::cos(k*rc);
    }
    inline T Xk(T k, T rc) {
      return -NormFactor/(k*k)* std::cos(k*rc);
    }
    inline T df(T r) { return -1.0/r/r;}
  };


  struct LRCoulombSingleton {

    typedef LRHandlerTemp<CoulombFunctor<double>,LPQHIBasis> LRHandlerType;
    
    static LRHandlerType* CoulombHandler;

    static LRHandlerType* getHandler(ParticleSet& ref);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
