//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_ONEDIMGRIDFACTORY_H
#define QMCPLUSPLUS_ONEDIMGRIDFACTORY_H
#include "Configuration.h"
#include "Numerics/OneDimGridFunctor.h"

namespace qmcplusplus
{

/** Factory class using Singleton pattern
 */
struct OneDimGridFactory: public QMCTraits
{

  ///typedef of the one-dimensional grid
  typedef OneDimGridBase<RealType>   GridType;
  ///typedef of map(string,GridType*>
  typedef std::map<string,GridType*> GridObjectMapType;

  ///container of one-dimensional grids
  static GridObjectMapType GridObjects;

  /** return a GridType*
   * @param cur xmlnode for the grid definition
   */
  static GridType* createGrid(xmlNodePtr cur);

  static RealType setSmoothCutoff(GridType* agrid, xmlNodePtr cur);
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
