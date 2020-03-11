#ifndef FS_UTILITIES_H
#define FS_UTILITIES_H

#include <vector>
#include "Configuration.h"

namespace qmcplusplus
{
typedef QMCTraits::PosType PosType;
typedef QMCTraits::RealType RealType;
typedef QMCTraits::IndexType IndexType;

void get_gridinfo_from_posgrid(const std::vector<PosType>& posgridlist, //list of grid points
                               const IndexType& axis,                   //the axis to get grid info.  0=x, 1=y,etc
                               RealType& lx,                            //the lower bound of grid (aka "left").
                               RealType& rx,                            //the upper bound (aka "right")
                               RealType& dx,                            // the grid spacing
                               IndexType& Nx);                          // the number of grid points

void getStats(const std::vector<RealType>& vals, RealType& avg, RealType& err);

} // namespace qmcplusplus
#endif
