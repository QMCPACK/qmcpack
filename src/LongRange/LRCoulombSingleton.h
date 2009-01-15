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
#include "Numerics/OneDimLinearSpline.h"

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
    void reset(ParticleSet& ref) {
      NormFactor=4.0*M_PI/ref.Lattice.Volume;
    }
    void reset(ParticleSet& ref, T rs) {
      NormFactor=4.0*M_PI/ref.Lattice.Volume;
    }
    inline T operator()(T r, T rinv) { return rinv;}
    inline T df(T r) { return -1.0/(r*r); }
    inline T Fk(T k, T rc) {
      return NormFactor/(k*k)* std::cos(k*rc);
    }
    inline T Xk(T k, T rc) {
      return -NormFactor/(k*k)* std::cos(k*rc);
    }
  };

  template<class T=double>
  struct PseudoCoulombFunctor {
    //typedef OneDimCubicSpline<T> RadialFunctorType;
    typedef OneDimLinearSpline<T> RadialFunctorType;
    RadialFunctorType& radFunc;
    T NormFactor;
    inline PseudoCoulombFunctor(RadialFunctorType& rfunc):radFunc(rfunc){}
    void reset(ParticleSet& ref) {
      NormFactor=4.0*M_PI/ref.Lattice.Volume;
    }
    inline T operator()(T r, T rinv) { return radFunc.splint(r);}
    inline T df(T r) {
      T du, d2u;
      radFunc.splint(r, du, d2u);
      return du;
    }
    inline T Fk(T k, T rc) {
      return NormFactor/(k*k)* std::cos(k*rc);
    }
    inline T Xk(T k, T rc) {
      return -NormFactor/(k*k)* std::cos(k*rc);
    }
  };


  struct LRCoulombSingleton {

    typedef OHMMS_PRECISION                                    RealType;
    typedef LRHandlerTemp<CoulombFunctor<RealType>,LPQHIBasis> LRHandlerType;
    typedef LinearGrid<RealType>                               GridType;
    //    typedef OneDimLinearSpline<RealType>                       RadFunctorType;
    typedef OneDimCubicSpline<RealType>                       RadFunctorType;

    static LRHandlerType* CoulombHandler;

    static LRHandlerType* getHandler(ParticleSet& ref);

    /** create a linear spline function
     * @param aLR LRHandler
     * @param rcut cutoff radius
     * @param agrid pointer to a grid
     * @return a RadFunctorType
     *
     * The spline function is the short-range term after breaking up
     * \f$r V_{S} = r \times (V(r)-V_{L})\f$
     */
    static RadFunctorType* createSpline4RbyVs(LRHandlerType* aLR, RealType rcut,
        GridType* agrid=0);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
