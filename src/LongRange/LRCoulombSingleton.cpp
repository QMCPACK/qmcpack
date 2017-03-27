//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/** @file LRHandlerTemp.h
 * @brief Define a LRHandler with two template parameters
 */
#include "LongRange/LRCoulombSingleton.h"
#if OHMMS_DIM==3
#include "LongRange/EwaldHandler.h"
#elif OHMMS_DIM==2
#include "LongRange/TwoDEwaldHandler.h"
#endif
#include <numeric>
namespace qmcplusplus
{

//initialization of the static data
LRCoulombSingleton::LRHandlerType* LRCoulombSingleton::CoulombHandler=0;
LRCoulombSingleton::LRHandlerType* LRCoulombSingleton::CoulombDerivHandler=0;
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
struct CoulombFunctor
{
  T NormFactor;
  inline CoulombFunctor() {}
  void reset(ParticleSet& ref)
  {
    NormFactor=4.0*M_PI/ref.LRBox.Volume;
  }
  void reset(ParticleSet& ref, T rs)
  {
    NormFactor=4.0*M_PI/ref.LRBox.Volume;
  }
  inline T operator()(T r, T rinv)
  {
    return rinv;
  }
  inline T df(T r)
  {
    return -1.0/(r*r);
  }
  inline T df2(T r)
  {
    return 2.0/(r*r*r);
  }
  inline T Vk(T k)
  {
    return NormFactor/(k*k);
  }
  
  inline T Xk_dk(T k){ return 0.0;}
  inline T Fk(T k, T rc)
  {
    return NormFactor/(k*k)* std::cos(k*rc);
  }
  inline T Xk(T k, T rc)
  {
    return -NormFactor/(k*k)* std::cos(k*rc);
  }
  
  inline T dVk_dk(T k)
  {
	  return -2*NormFactor/k/k/k;
  }
  inline T dFk_dk(T k, T rc)
  {
    return -NormFactor/k/k*(2.0/k*std::cos(k*rc)+rc*std::sin(k*rc));
  }
  
  inline T dXk_dk(T k, T rc)
  {
    return NormFactor/k/k*(2.0/k*std::cos(k*rc)+rc*std::sin(k*rc));
  }

  inline T integrate_r2(T r) const
  {
    return 0.5*r*r;
  }
};
#elif OHMMS_DIM==2
template<class T=double>
struct CoulombFunctor
{
  T NormFactor;
  inline CoulombFunctor() {}
  void reset(ParticleSet& ref)
  {
    NormFactor=2.0*M_PI/ref.LRBox.Volume;
  }
  void reset(ParticleSet& ref, T rs)
  {
    NormFactor=2.0*M_PI/ref.LRBox.Volume;
  }
  inline T operator()(T r, T rinv)
  {
    return rinv;
  }
  inline T df(T r)
  {
    return -1.0/(r*r);
  }
  inline T df2(T r)
  {
    return 2/(r*r*r);
  }
  inline T Fk(T k, T rc)
  {
    return NormFactor/k* std::cos(k*rc);
  }
  inline T Xk(T k, T rc)
  {
    return -NormFactor/k* std::cos(k*rc);
  }
  
  
 
  inline T integrate_r2(T r) const
  {
    return 0.5*r*r;
  }
};
#endif

template<class T=double>
struct PseudoCoulombFunctor
{
  //typedef OneDimCubicSpline<T> RadialFunctorType;
  typedef OneDimLinearSpline<T> RadialFunctorType;
  RadialFunctorType& radFunc;
  T NormFactor;
  inline PseudoCoulombFunctor(RadialFunctorType& rfunc):radFunc(rfunc) {}
  void reset(ParticleSet& ref)
  {
    NormFactor=4.0*M_PI/ref.LRBox.Volume;
  }
  inline T operator()(T r, T rinv)
  {
    return radFunc.splint(r);
  }
  inline T df(T r)
  {
    T du, d2u;
    radFunc.splint(r, du, d2u);
    return du;
  }
  inline T Fk(T k, T rc)
  {
    return NormFactor/(k*k)* std::cos(k*rc);
  }
  inline T Xk(T k, T rc)
  {
    return -NormFactor/(k*k)* std::cos(k*rc);
  }
  inline T integrate_r2(T r) const
  {
    //fix this to return the integration
    return 0.5*r*r;
  }
};


LRCoulombSingleton::LRHandlerType*
LRCoulombSingleton::getHandler(ParticleSet& ref)
{
  if(CoulombHandler ==0)
  {
#if OHMMS_DIM==3
    if(ref.SK->SuperCellEnum == SUPERCELL_SLAB)
    {
      app_log() << "\n   Creating CoulombHandler using quasi-2D Ewald method for the slab. " << std::endl;
      CoulombHandler= new EwaldHandler(ref);
    }
    else //if(ref.LRBox.SuperCellEnum == SUPERCELL_BULK)
    {
      app_log() << "\n  Creating CoulombHandler with the optimal breakup. " << std::endl;
      CoulombHandler= new LRHandlerTemp<CoulombFunctor<mRealType>,LPQHIBasis>(ref);
    //  CoulombHandler = new LRHandlerSRCoulomb<CoulombFunctor<mRealType>, LPQHISRCoulombBasis>(ref);
      
    }
//        else if(ref.LRBox.SuperCellEnum == SUPERCELL_SLAB)
//        {
//          app_log() << "\n   Creating CoulombHandler using quasi-2D Ewald method for the slab. " << std::endl;
//          CoulombHandler= new EwaldHandler(ref);
//        }
#elif OHMMS_DIM==2
    app_log() << "\n   Creating CoulombHandler using 2D Ewald method. " << std::endl;
    CoulombHandler= new TwoDEwaldHandler(ref);
#endif
    CoulombHandler->initBreakup(ref);
    return CoulombHandler;
  }
  else
  {
    app_log() << "  Clone CoulombHandler. " << std::endl;
    return CoulombHandler->makeClone(ref);
  }
}
LRCoulombSingleton::LRHandlerType*
LRCoulombSingleton::getDerivHandler(ParticleSet& ref)
{
  //APP_ABORT("SR Coulomb Basis Handler has cloning issues.  Stress also has some kinks");
  if(CoulombDerivHandler==0)
  {
    app_log() << "\n  Creating CoulombHandler with the optimal breakup of SR piece. " << std::endl;
    CoulombDerivHandler= new LRHandlerSRCoulomb<CoulombFunctor<mRealType>,LPQHISRCoulombBasis>(ref);
   // CoulombDerivHandler = new LRDerivHandler<CoulombFunctor<mRealType>, LPQHIBasis> (ref);
    //CoulombDerivHandler= new EwaldHandler(ref);
    CoulombDerivHandler->initBreakup(ref);
    //return CoulombDerivHandler;

    return CoulombDerivHandler;
  }
  else
  {
    app_log() << "  Clone CoulombDerivHandler. " << std::endl;
    return CoulombDerivHandler->makeClone(ref);
   // return CoulombDerivHandler;
  }
}

template<typename T>
OneDimCubicSpline<T>*  
createSpline4RbyVs_temp(LRHandlerBase* aLR, T rcut, LinearGrid<T>* agrid)
{
  typedef OneDimCubicSpline<T> func_type;
  if(agrid==0)
  {
    agrid = new LinearGrid<T>;
    agrid->set(0.0,rcut,1001);
  }
  int ng=agrid->size();
  std::vector<T> v(ng);
  T r=(*agrid)[0];
  //check if the first point is not zero
  v[0]=(r>std::numeric_limits<T>::epsilon())? r*aLR->evaluate(r,1.0/r):0.0;
  for(int ig=1; ig<ng-1; ig++)
  {
    r=(*agrid)[ig];
    v[ig]=r*aLR->evaluate(r,1.0/r);
  }
  v[0] = 2.0*v[1] - v[2];
  v[ng-1]=0.0;
  func_type* V0=new func_type(agrid,v);
  T deriv=(v[1]-v[0])/((*agrid)[1]-(*agrid)[0]);
  V0->spline(0,deriv,ng-1,0.0);
  return V0;
}

LRCoulombSingleton::RadFunctorType*
LRCoulombSingleton::createSpline4RbyVs(LRHandlerType* aLR, mRealType rcut,
                                       GridType* agrid)
{
  return createSpline4RbyVs_temp(aLR,static_cast<pRealType>(rcut),agrid);
}

}
