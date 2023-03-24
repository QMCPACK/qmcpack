//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
//
// File created by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "TwoBodyJastrow.h"

namespace qmcplusplus
{

template<typename FT>
template<class T>
T TwoBodyJastrow<FT>::do_ratioT(ParticleSet& P, int iat)
{
  //only ratio, ready to compute it again
  UpdateMode = ORB_PBYP_RATIO;
  cur_Uat    = computeU(P, iat, P.getDistTableAA(my_table_ID_).getTempDists());
  return std::exp(static_cast<T>(Uat[iat] - cur_Uat));
}


} // namespace qmcplusplus