//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


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
  ///typedef of map( std::string,GridType*>
  typedef std::map<std::string,GridType*> GridObjectMapType;

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
