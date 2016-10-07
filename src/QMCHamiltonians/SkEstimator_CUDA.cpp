//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "QMCHamiltonians/SkEstimator_CUDA.h"
#include "Particle/MCWalkerConfiguration.h"

namespace qmcplusplus
{
void
SkEstimator_CUDA::addEnergy(MCWalkerConfiguration &W,
                            std::vector<RealType> &LocalEnergy)
{
  int nw = W.WalkerList.size();
  gpu::host_vector<CUDA_PRECISION_FULL> rhok_host;
  int stride = (int)(W.WalkerList[0]->get_rhok_ptr(1) -
                     W.WalkerList[0]->get_rhok_ptr(0));
  RealType OneOverNW = 1.0/(RealType)nw;
  std::vector<CUDA_PRECISION> rhok_total(2*NumK, 0.0);
  for (int iw=0; iw<nw; iw++)
  {
    Walker_t &walker = *(W.WalkerList[iw]);
    rhok_host = walker.Rhok_GPU;
    fill (rhok_total.begin(), rhok_total.end(), 0.0);
    for (int isp=0; isp<NumSpecies; isp++)
      for (int ik=0; ik<2*NumK; ik++)
        rhok_total[ik] += rhok_host[ik+isp*stride];
    if(hdf5_out)
    {
      for (int ik=0; ik<NumK; ik++)
        W.Collectables[myIndex+ik] += OneOverN*OneOverNW*
                                      (rhok_total[2*ik+0]*rhok_total[2*ik+0] +
                                       rhok_total[2*ik+1]*rhok_total[2*ik+1]);
    }
    else
    {
      for (int ik=0; ik<NumK; ik++)
        W.WalkerList[iw]->getPropertyBase()[NUMPROPERTIES+myIndex+ik] =
          OneOverN * (rhok_total[2*ik+0]*rhok_total[2*ik+0] +
                      rhok_total[2*ik+1]*rhok_total[2*ik+1]);
    }
  }
}
}
