#ifndef GRID_UTILITIES_H
#define GRID_UTILITIES_H

#include <vector>
#include <algorithm>
#include "Configuration.h"

namespace qmcplusplus
{

typedef QMCTraits::PosType PosType;
typedef QMCTraits::RealType RealType;
typedef QMCTraits::IndexType IndexType;

void get_gridinfo_from_posgrid(const std::vector<PosType>& posgridlist,  //list of grid points
                                          const IndexType&        axis,  //the axis to get grid info.  0=x, 1=y,etc
                                                 RealType&          lx,  //the lower bound of grid (aka "left").  
                                                 RealType&          rx,  //the upper bound (aka "right")
                                                 RealType&          dx,  // the grid spacing
                                                IndexType&          Nx)  // the number of grid points
{
  vector<RealType> kx;
  kx.resize(posgridlist.size());

  for(IndexType i=0; i<posgridlist.size(); i++)
    kx[i]=posgridlist[i][axis];

  vector<RealType>::iterator it;

  std::sort(kx.begin(),kx.end());
  
  it=std::unique(kx.begin(),kx.end());

  lx=*(kx.begin());
  rx=*(it-1);
  Nx=it-kx.begin();
  dx=(rx-lx)/RealType(Nx-1);
}

}
#endif
