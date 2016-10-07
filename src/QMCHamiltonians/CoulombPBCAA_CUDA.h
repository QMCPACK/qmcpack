//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
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
  gpu::device_vector<CUDA_PRECISION_FULL>  SumGPU;
  gpu::host_vector<CUDA_PRECISION_FULL>  SumHost;
  gpu::device_vector<CUDA_PRECISION_FULL>  L, Linv;
  //// Long-range part
  int Numk;
  gpu::device_vector<CUDA_PRECISION_FULL> kpointsGPU;
  gpu::device_vector<int>            kshellGPU;
  // This has the same lengths as KshellGPU
  gpu::device_vector<CUDA_PRECISION_FULL> FkGPU;
  // The first vector index is the species number
  // Complex, stored as float2
  std::vector<gpu::device_vector<CUDA_PRECISION_FULL*> > RhoklistsGPU;
  std::vector<gpu::host_vector<CUDA_PRECISION_FULL*> > RhoklistsHost;
  gpu::device_vector<CUDA_PRECISION_FULL> RhokGPU;
  void setupLongRangeGPU(ParticleSet &P);
  void addEnergy(MCWalkerConfiguration &W,
                 std::vector<RealType> &LocalEnergy);

  void initBreakup(ParticleSet& P, bool cloning);
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

};
}
#endif
