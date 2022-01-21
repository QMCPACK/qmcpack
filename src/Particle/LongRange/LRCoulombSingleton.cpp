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
#include "LRCoulombSingleton.h"
#if OHMMS_DIM == 3
#include "LongRange/EwaldHandler.h"
#include "LongRange/EwaldHandler3D.h"
#elif OHMMS_DIM == 2
#include "LongRange/TwoDEwaldHandler.h"
#endif
#include <numeric>
namespace qmcplusplus
{
//initialization of the static data
std::unique_ptr<LRCoulombSingleton::LRHandlerType> LRCoulombSingleton::CoulombHandler;
std::unique_ptr<LRCoulombSingleton::LRHandlerType> LRCoulombSingleton::CoulombDerivHandler;
LRCoulombSingleton::lr_type LRCoulombSingleton::this_lr_type = ESLER;
/** CoulombFunctor
 *
 * An example for a Func for LRHandlerTemp. Four member functions have to be provided
 * - reset(T volume) : reset the normalization factor
 * - operator()(T r, T rinv): return a value of the original function, e.g., 1.0/r
 * - Fk(T k, T rc)
 * - Xk(T k, T rc)
 */
#if OHMMS_DIM == 3
template<class T = double>
struct CoulombFunctor
{
  T NormFactor;
  inline CoulombFunctor() {}
  void reset(ParticleSet& ref) { NormFactor = 4.0 * M_PI / ref.getLRBox().Volume; }
  void reset(ParticleSet& ref, T rs) { NormFactor = 4.0 * M_PI / ref.getLRBox().Volume; }
  inline T operator()(T r, T rinv) const { return rinv; }
  inline T df(T r) const { return -1.0 / (r * r); }
  inline T Vk(T k) const { return NormFactor / (k * k); }

  inline T Fk(T k, T rc) const { return NormFactor / (k * k) * std::cos(k * rc); }
  inline T Xk(T k, T rc) const { return -NormFactor / (k * k) * std::cos(k * rc); }

  inline T dVk_dk(T k) const { return -2 * NormFactor / k / k / k; }

  inline T integrate_r2(T r) const { return 0.5 * r * r; }
};
#elif OHMMS_DIM == 2
template<class T = double>
struct CoulombFunctor
{
  T NormFactor;
  inline CoulombFunctor() {}
  void reset(ParticleSet& ref) { NormFactor = 2.0 * M_PI / ref.getLRBox().Volume; }
  void reset(ParticleSet& ref, T rs) { NormFactor = 2.0 * M_PI / ref.getLRBox().Volume; }
  inline T operator()(T r, T rinv) const { return rinv; }
  inline T df(T r) const { return -1.0 / (r * r); }
  inline T Fk(T k, T rc) const { return NormFactor / k * std::cos(k * rc); }
  inline T Xk(T k, T rc) const { return -NormFactor / k * std::cos(k * rc); }


  inline T integrate_r2(T r) const { return 0.5 * r * r; }
};
#endif


std::unique_ptr<LRCoulombSingleton::LRHandlerType> LRCoulombSingleton::getHandler(ParticleSet& ref)
{
  if (CoulombHandler == 0)
  {
#if OHMMS_DIM == 3
    if (ref.getSK().SuperCellEnum == SUPERCELL_SLAB)
    {
      app_log() << "\n   Creating CoulombHandler using quasi-2D Ewald method for the slab. " << std::endl;
      CoulombHandler = std::make_unique<EwaldHandler>(ref);
    }
    else //if(ref.getLRBox().SuperCellEnum == SUPERCELL_BULK)
    {
      if (this_lr_type == ESLER)
      {
        app_log() << "\n  Creating CoulombHandler with the Esler Optimized Breakup. " << std::endl;
        CoulombHandler = std::make_unique<LRHandlerTemp<CoulombFunctor<mRealType>, LPQHIBasis>>(ref);
      }
      else if (this_lr_type == EWALD)
      {
        app_log() << "\n  Creating CoulombHandler with the 3D Ewald Breakup. " << std::endl;
        CoulombHandler = std::make_unique<EwaldHandler3D>(ref);
      }
      else if (this_lr_type == NATOLI)
      {
        app_log() << "\n  Creating CoulombHandler with the Natoli Optimized Breakup. " << std::endl;
        CoulombHandler = std::make_unique<LRHandlerSRCoulomb<CoulombFunctor<mRealType>, LPQHISRCoulombBasis>>(ref);
      }
      else
      {
        APP_ABORT("\n  Long range breakup method not recognized.\n");
      }
    }
//        else if(ref.getLRBox().SuperCellEnum == SUPERCELL_SLAB)
//        {
//          app_log() << "\n   Creating CoulombHandler using quasi-2D Ewald method for the slab. " << std::endl;
//          CoulombHandler= new EwaldHandler(ref);
//        }
#elif OHMMS_DIM == 2
    app_log() << "\n   Creating CoulombHandler using 2D Ewald method. " << std::endl;
    CoulombHandler = std::make_unique<TwoDEwaldHandler>(ref);
#endif
    CoulombHandler->initBreakup(ref);
    return std::unique_ptr<LRHandlerType>(CoulombHandler->makeClone(ref));
  }
  else
  {
    app_log() << "  Clone CoulombHandler. " << std::endl;
    return std::unique_ptr<LRHandlerType>(CoulombHandler->makeClone(ref));
  }
}

std::unique_ptr<LRCoulombSingleton::LRHandlerType> LRCoulombSingleton::getDerivHandler(ParticleSet& ref)
{
#if OHMMS_DIM != 3
  APP_ABORT("energy derivative implemented for 3D only");
#endif
  //APP_ABORT("SR Coulomb Basis Handler has cloning issues.  Stress also has some kinks");
  if (CoulombDerivHandler == 0)
  {
    if (this_lr_type == EWALD)
    {
      app_log() << "\n  Creating CoulombDerivHandler with the 3D Ewald Breakup. " << std::endl;
      CoulombDerivHandler = std::make_unique<EwaldHandler3D>(ref);
    }
    else if (this_lr_type == NATOLI)
    {
      app_log() << "\n  Creating CoulombDerivHandler with the Natoli Optimized Breakup. " << std::endl;
      CoulombDerivHandler = std::make_unique<LRHandlerSRCoulomb<CoulombFunctor<mRealType>, LPQHISRCoulombBasis>>(ref);
    }
    else if (this_lr_type == ESLER)
    {
      APP_ABORT("\n  Derivatives are not supported with Esler Optimized Breakup.\n");
    }
    else
    {
      APP_ABORT("\n  Long range breakup method for derivatives not recognized.\n");
    }
    CoulombDerivHandler->initBreakup(ref);
    return std::unique_ptr<LRHandlerType>(CoulombDerivHandler->makeClone(ref));
  }
  else
  {
    app_log() << "  Clone CoulombDerivHandler. " << std::endl;
    return std::unique_ptr<LRHandlerType>(CoulombDerivHandler->makeClone(ref));
  }
}

template<typename T>
std::unique_ptr<OneDimCubicSpline<T>> createSpline4RbyVs_temp(LRHandlerBase* aLR, T rcut, const LinearGrid<T>* agrid)
{
  using func_type = OneDimCubicSpline<T>;
  std::unique_ptr<LinearGrid<T>> agrid_local;
  if (agrid == nullptr)
  {
    agrid_local = std::make_unique<LinearGrid<T>>();
    agrid_local->set(0.0, rcut, 1001);
    agrid = agrid_local.get();
  }
  const int ng = agrid->size();
  std::vector<T> v(ng);
  T r = (*agrid)[0];
  //check if the first point is not zero
  v[0] = (r > std::numeric_limits<T>::epsilon()) ? r * aLR->evaluate(r, 1.0 / r) : 0.0;
  for (int ig = 1; ig < ng - 1; ig++)
  {
    r     = (*agrid)[ig];
    v[ig] = r * aLR->evaluate(r, 1.0 / r);
  }
  v[0]      = 2.0 * v[1] - v[2];
  v[ng - 1] = 0.0;
  auto V0   = std::make_unique<func_type>(agrid->makeClone(), v);
  T deriv   = (v[1] - v[0]) / ((*agrid)[1] - (*agrid)[0]);
  V0->spline(0, deriv, ng - 1, 0.0);
  return V0;
}

template<typename T>
std::unique_ptr<OneDimCubicSpline<T>> createSpline4RbyVsDeriv_temp(LRHandlerBase* aLR,
                                                                   T rcut,
                                                                   const LinearGrid<T>* agrid)
{
  using func_type = OneDimCubicSpline<T>;
  std::unique_ptr<LinearGrid<T>> agrid_local;
  if (agrid == nullptr)
  {
    agrid_local = std::make_unique<LinearGrid<T>>();
    agrid_local->set(0.0, rcut, 1001);
    agrid = agrid_local.get();
  }
  int ng = agrid->size();
  std::vector<T> v(ng);
  T r = (*agrid)[0];
  //check if the first point is not zero
  v[0] = (r > std::numeric_limits<T>::epsilon()) ? r * aLR->evaluate(r, 1.0 / r) : 0.0;
  T v_val(0.0);
  T v_deriv(0.0);

  for (int ig = 1; ig < ng - 1; ig++)
  {
    r       = (*agrid)[ig];
    v_val   = aLR->evaluate(r, 1.0 / r);
    v_deriv = aLR->srDf(r, 1 / r);

    v[ig] = v_val + v_deriv * r;
  }
  v[0]      = 2.0 * v[1] - v[2];
  v[ng - 1] = 0.0;
  auto dV0  = std::make_unique<func_type>(agrid->makeClone(), v);
  T deriv   = (v[1] - v[0]) / ((*agrid)[1] - (*agrid)[0]);
  dV0->spline(0, deriv, ng - 1, 0.0);
  return dV0;
}


std::unique_ptr<LRCoulombSingleton::RadFunctorType> LRCoulombSingleton::createSpline4RbyVs(LRHandlerType* aLR,
                                                                                           mRealType rcut,
                                                                                           const GridType* agrid)
{
  return createSpline4RbyVs_temp(aLR, static_cast<pRealType>(rcut), agrid);
}

std::unique_ptr<LRCoulombSingleton::RadFunctorType> LRCoulombSingleton::createSpline4RbyVsDeriv(LRHandlerType* aLR,
                                                                                                mRealType rcut,
                                                                                                const GridType* agrid)
{
  return createSpline4RbyVsDeriv_temp(aLR, static_cast<pRealType>(rcut), agrid);
}


} // namespace qmcplusplus
