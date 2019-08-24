//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from ParticleSet.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MOVECONTEXT_H
#define QMCPLUSPLUS_MOVECONTEXT_H

#include <map>
#include <vector>
#include "OhmmsSoA/VectorSoaContainer.h"
#include "Configuration.h"
#include "Particle/Walker.h"

namespace qmcplusplus
{

class DistanceTableData;
/** Thread local context for moving walkers
 *
 *  created on the stack once per driver per crowd
 *  refactored out of ParticleSet, only the minimum
 */
class MoveContext
{
public:
  using ParticlePositions = PtclOnLatticeTraits::ParticlePos_t;
  using MCPWalker  = Walker<QMCTraits, PtclOnLatticeTraits>;
  using RealType = QMCTraits::RealType;

  MoveContext(int particles);
  
  void loadWalker(MCPWalker& awalker);

  ParticlePositions get_positions() const { return positions_; }
protected:
// only one of these should exist
  ///Positions
  ParticlePositions positions_;

  ///SoA copy of R
  VectorSoaContainer<RealType, OHMMS_DIM> positions_soa_;

  /** map to handle distance tables
   *
   * myDistTableMap[source-particle-tag]= locator in the distance table
   * myDistTableMap[ObjectTag] === 0
   */
  std::map<std::string, int> myDistTableMap;

  /// distance tables that need to be updated by moving this ParticleSet
  std::vector<DistanceTableData*> DistTables;

  /// Descriptions from distance table creation.  Same order as DistTables.
  std::vector<std::string> distTableDescriptions;

};

}
#endif
