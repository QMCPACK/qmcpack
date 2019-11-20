//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file LRCoulombSingleton.h
 * @brief Define a LRHandler with two template parameters
 */
#ifndef QMCPLUSPLUS_LRCOULOMBSINGLETON_H
#define QMCPLUSPLUS_LRCOULOMBSINGLETON_H

#include <config.h>
#include "LongRange/LRHandlerTemp.h"
#include "LongRange/LRHandlerSRCoulomb.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/OneDimLinearSpline.h"

namespace qmcplusplus
{
struct LRCoulombSingleton
{
  DECLARE_COULOMB_TYPES

  typedef LRHandlerBase LRHandlerType;
  typedef LinearGrid<pRealType> GridType;
  typedef OneDimCubicSpline<pRealType> RadFunctorType;

  ///Stores the energ optimized LR handler.
  static LRHandlerType* CoulombHandler;
  ///Stores the force/stress optimized LR handler.
  static LRHandlerType* CoulombDerivHandler;
  ///This returns an energy optimized LR handler.  If non existent, it creates one.
  static LRHandlerType* getHandler(ParticleSet& ref);
  ///This returns a force/stress optimized LR handler.  If non existent, it creates one.
  static LRHandlerType* getDerivHandler(ParticleSet& ref, std::string method="ewald");

  //The following two helper functions are provided to spline the short-range component
  //of the coulomb potential and its derivative.  This is much faster than evaluating
  //via a sum over basis functions, which is typical in Ceperley/Esler style breakups.

  /** create a linear spline function
   * @param aLR LRHandler
   * @param rcut cutoff radius
   * @param agrid pointer to a grid
   * @return a RadFunctorType
   *
   * The spline function is the short-range term after breaking up
   * \f$r V_{S} = r \times (V(r)-V_{L})\f$
   */
  static RadFunctorType* createSpline4RbyVs(LRHandlerType* aLR, mRealType rcut, GridType* agrid = 0);
  /** create a linear spline of the derivative of short-range potential
   * @param aLR LRHandler
   * @param rcut cutoff radius
   * @param agrid pointer to a grid
   * @return a RadFunctorType
   *
   * The spline function is the short-range term after breaking up
   * \f$r \frac{d}{dr} V_{S} = \frac{d}{dr}\left(r \times (V(r)-V_{L})\right)\f$
   */
  static RadFunctorType* createSpline4RbyVsDeriv(LRHandlerType* aLR, mRealType rcut, GridType* agrid = 0);
};

} // namespace qmcplusplus
#endif
