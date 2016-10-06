//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
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
                 std::vector<RealType> &LocalEnergy);

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
                 std::vector<RealType> &LocalEnergy);

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);
};

}
#endif
