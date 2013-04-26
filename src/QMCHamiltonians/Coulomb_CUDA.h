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
#ifndef QMCPLUSPLUS_COULOMBPBCAA_CUDA_H
#define QMCPLUSPLUS_COULOMBPBCAA_CUDA_H
#include "QMCHamiltonians/CoulombPotential.h"
#include "QMCHamiltonians/CudaCoulomb.h"

class MCWalkerConfiguration;

namespace qmcplusplus
{
struct CoulombAA_CUDA : public CoulombPotentialAA
{
  //////////////////////////////////
  // Vectorized evaluation on GPU //
  //////////////////////////////////
  CoulombAA_CUDA(ParticleSet& ref);

  ParticleSet &PtclRef;
  gpu::device_vector<CUDA_PRECISION>  SumGPU;
  gpu::host_vector<CUDA_PRECISION>  SumHost;
  void addEnergy(MCWalkerConfiguration &W,
                 vector<RealType> &LocalEnergy);

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);
};

struct CoulombAB_CUDA : public CoulombPotentialAB
{
  ParticleSet &PtclRef;
  //////////////////////////////////
  // Vectorized evaluation on GPU //
  //////////////////////////////////
  CoulombAB_CUDA(ParticleSet& ref, ParticleSet &ions);

  gpu::device_vector<CUDA_PRECISION>  SumGPU;
  gpu::host_vector<CUDA_PRECISION>  SumHost;
  gpu::device_vector<CUDA_PRECISION> ZionGPU;
  void addEnergy(MCWalkerConfiguration &W,
                 vector<RealType> &LocalEnergy);

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);
};

}
#endif
