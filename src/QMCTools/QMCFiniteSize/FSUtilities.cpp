#include "FSUtilities.h"
#include <algorithm>

namespace qmcplusplus
{

void get_gridinfo_from_posgrid(const std::vector<PosType>& posgridlist, 
                               const IndexType& axis,                   
                               RealType& lx,                            
                               RealType& rx,                            
                               RealType& dx,                            
                               IndexType& Nx)                           
{
  std::vector<RealType> kx;
  kx.resize(posgridlist.size());

  for (IndexType i = 0; i < posgridlist.size(); i++)
    kx[i] = posgridlist[i][axis];

  std::vector<RealType>::iterator it;

  std::sort(kx.begin(), kx.end());

  it = std::unique(kx.begin(), kx.end());

  lx = *(kx.begin());
  rx = *(it - 1);
  Nx = it - kx.begin();
  dx = (rx - lx) / RealType(Nx - 1);
}

void getStats(const std::vector<RealType>& vals, RealType& avg, RealType& err)
{
  avg = 0.0;
  for (int i = 0; i < vals.size(); i++)
  {
    avg += vals[i];
  }
  avg /= vals.size();
  err = 0.0;
  for (int i = 0; i < vals.size(); i++)
  {
    err += (vals[i] - avg) * (vals[i] - avg);
  }
  err /= vals.size();
  err = std::sqrt(err);
}

} // namespace qmcplusplus
