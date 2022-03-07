//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "DynamicCoordinatesBuilder.h"
#include "Particle/RealSpacePositions.h"
#include "Particle/RealSpacePositionsOMPTarget.h"

namespace qmcplusplus
{
/** create DynamicCoordinates based on kind
 */
std::unique_ptr<DynamicCoordinates> createDynamicCoordinates(const DynamicCoordinateKind kind)
{
  if (kind == DynamicCoordinateKind::DC_POS)
    return std::make_unique<RealSpacePositions>();
  else if (kind == DynamicCoordinateKind::DC_POS_OFFLOAD)
    return std::make_unique<RealSpacePositionsOMPTarget>();
  // dummy return
  return std::unique_ptr<RealSpacePositions>();
}
} // namespace qmcplusplus
