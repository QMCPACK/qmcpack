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
    
#ifndef LOCAL_ECPOTENTIAL_CUDA_H
#define LOCAL_ECPOTENTIAL_CUDA_H

#include "QMCHamiltonians/LocalECPotential.h"
#include "QMCHamiltonians/CudaCoulomb.h"

class MCWalkerConfiguration;

namespace qmcplusplus
{
struct LocalECPotential_CUDA : public LocalECPotential
{
  //////////////////////////////////
  // Vectorized evaluation on GPU //
  //////////////////////////////////
  //// Short-range part
  int NumIons, NumElecs, NumElecGroups, NumIonSpecies;
  ParticleSet &ElecRef, &IonRef;
  std::vector<int> IonFirst, IonLast;
  // This is indexed by the ion species
  std::vector<TextureSpline*> SRSplines;
  TextureSpline *V0Spline;
  gpu::device_vector<CUDA_PRECISION>  SumGPU;
  gpu::host_vector<CUDA_PRECISION>    SumHost;
  gpu::device_vector<CUDA_PRECISION>  IGPU;
  gpu::device_vector<CUDA_PRECISION>  ZionGPU;

  std::vector<PosType> SortedIons;
  void add(int groupID, RadialPotentialType* ppot, RealType zion);

  void addEnergy(MCWalkerConfiguration &W,
                 std::vector<RealType> &LocalEnergy);

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  LocalECPotential_CUDA(ParticleSet& ions, ParticleSet& elns);

};
}

#endif
