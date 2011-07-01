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
/** @file LRHandlerTemp.h
 * @brief Define a LRHandler with two template parameters
 */
#include "LongRange/LRCoulombSingleton.h"
#include "LongRange/EwaldHandler.h"
#include <numeric>

namespace qmcplusplus {

  //initialization of the static data
  LRCoulombSingleton::LRHandlerType* LRCoulombSingleton::CoulombHandler=0;

  /** CoulombFunctor
   *
   * An example for a Func for LRHandlerTemp. Four member functions have to be provided
   * - reset(T volume) : reset the normalization factor
   * - operator()(T r, T rinv): return a value of the original function, e.g., 1.0/r
   * - Fk(T k, T rc)
   * - Xk(T k, T rc)
   */
#if OHMMS_DIM==3  
  template<class T=double>
  struct CoulombFunctor {
    T NormFactor;
    inline CoulombFunctor(){}
    void reset(ParticleSet& ref) 
    {
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

    inline T integrate_r2(T r) const 
    {
      return 0.5*r*r;
    }
  };
#elif OHMMS_DIM==2
    template<class T=double>
  struct CoulombFunctor {
    T NormFactor;
    inline CoulombFunctor(){}
    void reset(ParticleSet& ref) 
    {
      NormFactor=2.0*M_PI/ref.Lattice.Volume;
    }
    void reset(ParticleSet& ref, T rs) {
      NormFactor=2.0*M_PI/ref.Lattice.Volume;
    }
    inline T operator()(T r, T rinv) { return rinv;}
    inline T df(T r) { return -1.0/(r*r); }
    inline T Fk(T k, T rc) {
      return NormFactor/k* std::cos(k*rc);
    }
    inline T Xk(T k, T rc) {
      return -NormFactor/k* std::cos(k*rc);
    }

    inline T integrate_r2(T r) const 
    {
      return 0.5*r*r;
    }
  };
#endif  

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
    inline T integrate_r2(T r) const 
    {
      //fix this to return the integration
      return 0.5*r*r;
    }
  };


  LRCoulombSingleton::LRHandlerType*
    LRCoulombSingleton::getHandler(ParticleSet& ref) {
      if(CoulombHandler ==0) 
      {
#if OHMMS_DIM==3  
        if(ref.Lattice.SuperCellEnum == SUPERCELL_BULK)
        {
          app_log() << "\n  Creating CoulombHandler with the optimal breakup. " << endl;
          CoulombHandler= new LRHandlerTemp<CoulombFunctor<RealType>,LPQHIBasis>(ref);
        }
        else if(ref.Lattice.SuperCellEnum == SUPERCELL_SLAB)
        {
          app_log() << "\n   Creating CoulombHandler using quasi-2D Ewald method for the slab. " << endl;
          CoulombHandler= new EwaldHandler(ref);
        }
#elif OHMMS_DIM==2
          app_log() << "\n   Creating CoulombHandler using quasi-2D Ewald method for the slab. " << endl;
          CoulombHandler= new TwoDEwaldHandler(ref);
#endif

        CoulombHandler->initBreakup(ref);
        return CoulombHandler;
      }
      else
      {
        app_log() << "  Clone CoulombHandler. " << endl;
        return CoulombHandler->makeClone(ref);
      }
    }

  LRCoulombSingleton::RadFunctorType*
    LRCoulombSingleton::createSpline4RbyVs(LRHandlerType* aLR, RealType rcut,
        GridType* agrid)
    {
      if(agrid == 0) 
      {
        agrid = new GridType;
        agrid->set(0.0,rcut,1001);
      }

      int ng=agrid->size();
      vector<RealType> v(ng);
      RealType r=(*agrid)[0];

      //check if the first point is not zero
      v[0]=(r>numeric_limits<RealType>::epsilon())? r*aLR->evaluate(r,1.0/r):0.0; 
      for(int ig=1; ig<ng-1; ig++) {
        r=(*agrid)[ig];
        v[ig]=r*aLR->evaluate(r,1.0/r);
      }
      v[0] = 2.0*v[1] - v[2];

      v[ng-1]=0.0;
      RadFunctorType* V0=new RadFunctorType(agrid,v);
      RealType deriv=(v[1]-v[0])/((*agrid)[1]-(*agrid)[0]);
      V0->spline(0,deriv,ng-1,0.0);

      return V0;
    }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
