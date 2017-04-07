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

  typedef LRHandlerBase                LRHandlerType;
  typedef LinearGrid<pRealType>        GridType;
  typedef OneDimCubicSpline<pRealType> RadFunctorType;

  static LRHandlerType* CoulombHandler;
  static LRHandlerType* CoulombDerivHandler;
  static LRHandlerType* getHandler(ParticleSet& ref);
  static LRHandlerType* getDerivHandler(ParticleSet& ref);
  /** create a linear spline function
   * @param aLR LRHandler
   * @param rcut cutoff radius
   * @param agrid pointer to a grid
   * @return a RadFunctorType
   *
   * The spline function is the short-range term after breaking up
   * \f$r V_{S} = r \times (V(r)-V_{L})\f$
   */
  static RadFunctorType* createSpline4RbyVs(LRHandlerType* aLR, mRealType rcut,
      GridType* agrid=0);
};

}
#endif
