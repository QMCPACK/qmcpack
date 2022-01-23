#ifndef FS_UTILITIES_H
#define FS_UTILITIES_H

#include <vector>
#include "Configuration.h"

namespace qmcplusplus
{
using PosType   = QMCTraits::PosType;
using RealType  = QMCTraits::RealType;
using IndexType = QMCTraits::IndexType;

void get_gridinfo_from_posgrid(const std::vector<PosType>& posgridlist, //list of grid points
                               const IndexType& axis,                   //the axis to get grid info.  0=x, 1=y,etc
                               RealType& lx,                            //the lower bound of grid (aka "left").
                               RealType& rx,                            //the upper bound (aka "right")
                               RealType& dx,                            // the grid spacing
                               IndexType& Nx);                          // the number of grid points

/** Simpleaverage and error estimate
 * 
 * input array of values to be averaged. returns avg and err.
 * Can pass a start index to average from "start" to the end of the array
 */
void getStats(const std::vector<RealType>& vals, RealType& avg, RealType& err, int start = 0);

/** estimate equilibration of block data
 *
 * Assumes vals is an array to be averaged. This will estimate which, if any, of the 
 * data should be thrown out as equilibration data. 
 * Returns the index of the array to start averaging from. 
 * frac is the fraction of the data to average initially to help determine the equilibration
 */
int estimateEquilibration(const std::vector<RealType>& vals, RealType frac = 0.75);

} // namespace qmcplusplus
#endif
