//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
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

#include <config.h>
#include "LongRange/LRHandlerTemp.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/OneDimLinearSpline.h"

namespace qmcplusplus
{

struct LRCoulombSingleton
{

  typedef OHMMS_PRECISION                                    RealType;
  typedef LRHandlerBase                                      LRHandlerType;
  typedef LinearGrid<RealType>                               GridType;
  //    typedef OneDimLinearSpline<RealType>                 RadFunctorType;
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
