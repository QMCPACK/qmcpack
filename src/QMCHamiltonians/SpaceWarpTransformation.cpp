#include "QMCHamiltonians/SpaceWarpTransformation.h"

namespace qmcplusplus
{

SpaceWarpTransformation::RealType SpaceWarpTransformation::f(RealType r)
{
  return 1.0/(r*r*r*r);
}

}
