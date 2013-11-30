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
  vector<int> IonFirst, IonLast;
  vector<ParticleSet*> myPtcl;
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
                 vector<RealType> &LocalEnergy);

};
}

#endif
