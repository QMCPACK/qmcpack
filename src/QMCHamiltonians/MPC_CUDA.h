//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_MPC_CUDA_H
#define QMCPLUSPLUS_MPC_CUDA_H

#include "QMCHamiltonians/MPC.h"
#include <einspline/bspline_create_cuda.h>

namespace qmcplusplus
{

class MPC_CUDA: public MPC
{
protected:
  //////////////////////////////////
  // Vectorized evaluation on GPU //
  //////////////////////////////////
  //// Short-range part
  UBspline_3d_s_cuda *CudaSpline;
  std::vector<int> IonFirst, IonLast;
  std::vector<ParticleSet*> myPtcl;
  // This is indexed by the ion species
  gpu::device_vector<CUDA_PRECISION>  SumGPU;
  gpu::host_vector<CUDA_PRECISION>  SumHost;
  gpu::device_vector<CUDA_PRECISION>  L, Linv;

  void initBreakup();

public:
  MPC_CUDA(ParticleSet& ref, double cutoff);

  QMCHamiltonianBase* makeClone(ParticleSet& qp,
                                TrialWaveFunction& psi);

  void addEnergy(MCWalkerConfiguration &W,
                 std::vector<RealType> &LocalEnergy);

};
}

#endif
