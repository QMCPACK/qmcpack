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

void getStats(const std::vector<RealType>& vals, RealType& avg, RealType& err, int start)
{
  avg   = 0.0;
  int n = 0;
  for (int i = start; i < vals.size(); i++)
  {
    avg += vals[i];
    n++;
  }
  avg /= n;
  err = 0.0;
  for (int i = start; i < vals.size(); i++)
  {
    err += (vals[i] - avg) * (vals[i] - avg);
  }
  err /= n;
  err = std::sqrt(err);
}

int estimateEquilibration(const std::vector<RealType>& vals, RealType frac)
{
  int idx = int(frac * vals.size());
  RealType avg, err;
  getStats(vals, avg, err, idx);
  int c3, c2, c1;
  c3 = vals.size();
  c2 = vals.size();
  c1 = vals.size();
  for (int i = vals.size() - 2; i >= 0; i--)
  {
    if ((vals[i] - avg) * (vals[i + 1] - avg) < 0.0)
    {
      c3 = c2;
      c2 = c1;
      c1 = i;
    }
  }
  if (c3 > frac * vals.size())
    c3 = int(frac * vals.size());
  return c3;
}

} // namespace qmcplusplus
