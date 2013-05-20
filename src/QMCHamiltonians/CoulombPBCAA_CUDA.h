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
#include "QMCHamiltonians/CoulombPBCAA.h"
#include "QMCHamiltonians/CudaCoulomb.h"


namespace qmcplusplus
{
struct CoulombPBCAA_CUDA : public CoulombPBCAA
{
  //////////////////////////////////
  // Vectorized evaluation on GPU //
  //////////////////////////////////

  CoulombPBCAA_CUDA(ParticleSet& ref, bool active, bool cloning=false);

  ParticleSet &PtclRef;
  //// Short-range part
  TextureSpline *SRSpline;
  gpu::device_vector<CUDA_PRECISION>  SumGPU;
  gpu::host_vector<CUDA_PRECISION>  SumHost;
  gpu::device_vector<CUDA_PRECISION>  L, Linv;
  //// Long-range part
  int Numk;
  gpu::device_vector<CUDA_PRECISION> kpointsGPU;
  gpu::device_vector<int>            kshellGPU;
  // This has the same lengths as KshellGPU
  gpu::device_vector<CUDA_PRECISION> FkGPU;
  // The first vector index is the species number
  // Complex, stored as float2
  vector<gpu::device_vector<CUDA_PRECISION*> > RhoklistsGPU;
  vector<gpu::host_vector<CUDA_PRECISION*> > RhoklistsHost;
  gpu::device_vector<CUDA_PRECISION> RhokGPU;
  void setupLongRangeGPU(ParticleSet &P);
  void addEnergy(MCWalkerConfiguration &W,
                 vector<RealType> &LocalEnergy);

  void initBreakup(ParticleSet& P, bool cloning);
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

};
}
#endif
