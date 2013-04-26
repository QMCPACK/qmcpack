//////////////////////////////////////////////////////////////////
// (c) Copyright 2010-  by Ken Esler and Jeongnim Kim
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: esler@uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////

#include "QMCHamiltonians/SkEstimator_CUDA.h"
#include "Particle/MCWalkerConfiguration.h"

namespace qmcplusplus
{
void
SkEstimator_CUDA::addEnergy(MCWalkerConfiguration &W,
                            vector<RealType> &LocalEnergy)
{
  int nw = W.WalkerList.size();
  gpu::host_vector<CUDA_PRECISION> rhok_host;
  int stride = (int)(W.WalkerList[0]->get_rhok_ptr(1) -
                     W.WalkerList[0]->get_rhok_ptr(0));
  RealType OneOverNW = 1.0/(RealType)nw;
  vector<CUDA_PRECISION> rhok_total(2*NumK, 0.0);
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
