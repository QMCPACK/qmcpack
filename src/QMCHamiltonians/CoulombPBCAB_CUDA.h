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


#ifndef QMCPLUSPLUS_COULOMBPBCAB_CUDA_H
#define QMCPLUSPLUS_COULOMBPBCAB_CUDA_H
#include "QMCHamiltonians/CoulombPBCAB.h"
#include "QMCHamiltonians/CudaCoulomb.h"


namespace qmcplusplus
{
struct CoulombPBCAB_CUDA : public CoulombPBCAB
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
  TextureSpline* V0Spline;
  gpu::device_vector<CUDA_PRECISION_FULL> SumGPU;
  gpu::host_vector<CUDA_PRECISION_FULL> SumHost;
  gpu::device_vector<CUDA_PRECISION> IGPU;
  gpu::device_vector<CUDA_PRECISION_FULL> L, Linv;
  //// Long-range part
  int Numk;
  gpu::device_vector<CUDA_PRECISION_FULL> kpointsGPU;
  gpu::device_vector<int> kshellGPU;
  // This has the same lengths as KshellGPU
  gpu::device_vector<CUDA_PRECISION_FULL> FkGPU;
  // The first vector index is the species number
  // Complex, stored as float2
  // This is for the electrons -- one per walker
  gpu::device_vector<CUDA_PRECISION_FULL*> RhoklistGPU;
  gpu::host_vector<CUDA_PRECISION_FULL*> RhoklistHost;
  // This stores rho_k for the electrons in one big array
  gpu::device_vector<CUDA_PRECISION_FULL> RhokElecGPU;

  std::vector<PosType> SortedIons;
  // This stores rho_k for the ions.  Index is species number
  std::vector<gpu::device_vector<CUDA_PRECISION_FULL>> RhokIonsGPU;
  void setupLongRangeGPU();

  void add(int groupID, RadFunctorType* ppot);

  void initBreakup(ParticleSet& P);

  void addEnergy(MCWalkerConfiguration& W, std::vector<RealType>& LocalEnergy);

  OperatorBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  CoulombPBCAB_CUDA(ParticleSet& ions, ParticleSet& elns, bool cloning = false);
};
} // namespace qmcplusplus

#endif
